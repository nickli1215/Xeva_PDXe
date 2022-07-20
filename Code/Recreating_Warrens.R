library(here)
library(magrittr)
library(tidyverse)
library(reticulate)
library(SummarizedExperiment)
conda_list()
use_condaenv(condaenv = 'Xeva_Analysis', required = TRUE)
library(Seurat)
source('analysis_helpers.R')
source('global_params.R')
source('Celligner_helpers.R')
source('Celligner_methods.R')
#plan("multiprocess", workers = 5)

#pset_folder<-"../Data"
pset_folder<-"G:/UHN/psets"

TCGA_tissue_data<-read.csv("../Data/TCGA_experiment_data.csv")
rownames(TCGA_tissue_data)<-TCGA_tissue_data$TCGA.studyid
studies<-na.omit(TCGA_tissue_data$TCGA.studyid)
studies<-c("UCS","SKCM","PRAD","LUAD","OV","BRCA")

exprs_mats<-list()
TCGA_tissues<-c()
#gene_order<-readRDS(paste0(pset_folder,"/TCGA/",studies[1],"/RNASeq.rds"))[,1]
gene_order<-readRDS(paste0(pset_folder,"/TCGA/",studies[1],".rds"))[,1]

for (study in studies){
    tissue_type<-study
    #if(file.exists(paste0(pset_folder,"/TCGA/",study,"/RNASeq.rds"))){
    if(file.exists(paste0(pset_folder,"/TCGA/",study,".rds"))){
    #exprs_mat<-readRDS(paste0(pset_folder,"/TCGA/",study,"/RNASeq.rds"))
    exprs_mat<-readRDS(paste0(pset_folder,"/TCGA/",study,".rds"))
    rownames(exprs_mat)<-exprs_mat[,1]
    exprs_mat<-exprs_mat[,-1]
    exprs_mat<-exprs_mat[gene_order,]
    
    #TCGA_clinical_data<-readRDS(paste0(pset_folder,"/TCGA/",study,"_clinicalMatrix.rds"))
    
    exprs_mats[[study]]<-exprs_mat
    print(study)
    TCGA_tissues<-c(TCGA_tissues,
                    rep(TCGA_tissue_data[study,"CCLE.tissueid"],
                                     ncol(exprs_mat)))
    }
  }

TCGA_mat<-do.call(cbind,exprs_mats)



CCLE.exprs<-readRDS('../Data/CCLE_expression.RDS')
gene_metadata<-CCLE.exprs@elementMetadata
rownames(gene_metadata)<-gene_metadata$gene_id
CCLE_mat<-assay(CCLE.exprs)
CCLE_mat<-CCLE_mat-log2(0.001)
gene_metadata<-gene_metadata[rownames(CCLE_mat),]
rownames(CCLE_mat)<-gene_metadata$gene_name



PDXe<-readRDS("../Data/PDXe_RNASeq_FPKM.rds")
PDXe.metadata<-pData(PDXe)
PDXe_mat<-exprs(PDXe)
PDXe.metadata<-PDXe.metadata[colnames(PDXe_mat),]
PDXe_tissue_data<-read.csv("../Data/PDXe_experiment_data.csv")
rownames(PDXe_tissue_data)<-PDXe_tissue_data$PDXe.tissueid
PDXe.metadata[["CCLE.tissueid"]]<-PDXe_tissue_data[PDXe.metadata$tissue,]$CCLE.tissueid


gene_list<-intersect(rownames(CCLE_mat),rownames(TCGA_mat))
gene_list<-intersect(gene_list,rownames(PDXe_mat))
gene_list<-intersect(gene_list,gene_metadata[gene_metadata$gene_type=="protein_coding",]$gene_name)
CCLE_mat<-CCLE_mat[gene_list,colSums(is.na(CCLE_mat))==0]
TCGA_mat<-TCGA_mat[gene_list,colSums(is.na(TCGA_mat))==0]
PDXe_mat<-PDXe_mat[gene_list,colSums(is.na(PDXe_mat))==0]

CCLE.celldata<-colData(CCLE.exprs)
CCLE.celldata<-CCLE.celldata[colnames(CCLE_mat),]


ann <- data.frame(sampleID = c(colnames(TCGA_mat), colnames(PDXe_mat),colnames(CCLE_mat)),
                  lineage = c(TCGA_tissues,PDXe.metadata$CCLE.tissueid,CCLE.celldata$tissueid),
                  subtype = NA,
                  type = c(rep('tumor', ncol(TCGA_mat)),rep('PDX',ncol(PDXe_mat)), rep('CL', ncol(CCLE_mat))))
ann$`Primary/Metastasis` <- NA

TCGA_ann <- dplyr::filter(ann, type=='tumor')
CCLE_ann <- dplyr::filter(ann, type=='CL')
PDXe_ann <- dplyr::filter(ann, type=='PDX')

CCLE_mat<-t(CCLE_mat)
TCGA_mat<-t(TCGA_mat)
PDXe_mat<-t(PDXe_mat)

gene_stats<-data.frame(
  Tumor_SD = apply(TCGA_mat, 2, sd, na.rm=T),
  CCLE_SD = apply(CCLE_mat, 2, sd, na.rm=T),
  PDXe_SD=apply(PDXe_mat,2,sd,na.rm=T),
  Tumor_mean = colMeans(TCGA_mat, na.rm=T),
  CCLE_mean = colMeans(CCLE_mat, na.rm=T),
  PDXe_mean = colMeans(PDXe_mat, na.rm=T),
  Gene = gene_list,
  stringsAsFactors = F)%>% 
  dplyr::mutate(max_SD = pmax(Tumor_SD, CCLE_SD,PDXe_SD,na.rm=T))

TCGA_obj <- create_Seurat_object(TCGA_mat, TCGA_ann, type='tumor')
CCLE_obj <- create_Seurat_object(CCLE_mat, CCLE_ann, type='CL')
PDXe_obj <- create_Seurat_object(PDXe_mat, PDXe_ann, type='PDX')

TCGA_obj <- cluster_data(TCGA_obj)
CCLE_obj <- cluster_data(CCLE_obj)
PDXe_obj <- cluster_data(PDXe_obj)

tumor_DE_genes <- find_differentially_expressed_genes(TCGA_obj)
CL_DE_genes <- find_differentially_expressed_genes(CCLE_obj)
PDX_DE_genes <- find_differentially_expressed_genes(PDXe_obj)
DE_genes <- full_join(tumor_DE_genes, PDX_DE_genes, by = 'Gene', suffix = c('_tumor', '_PDX'))
DE_genes<-full_join(DE_genes,CL_DE_genes, by = 'Gene', suffix = c('_tumor', '_CL'))
colnames(DE_genes)[4]<-"gene_stat_CL"
DE_genes<-DE_genes%>%
  mutate(
    tumor_rank = dplyr::dense_rank(-gene_stat_tumor),
    CL_rank = dplyr::dense_rank(-gene_stat_CL),
    PDX_rank = dplyr::dense_rank(-gene_stat_PDX),
    best_rank = pmin(tumor_rank,PDX_rank,CL_rank, na.rm=T)) %>%
  dplyr::left_join(gene_stats, by = 'Gene')

# take genes that are ranked in the top 1000 from either dataset, used for finding mutual nearest neighbors
DE_gene_set <- DE_genes %>%
  dplyr::filter(best_rank < global$top_DE_genes_per) %>%
  .[['Gene']]

cov_diff_eig_TCGA <- run_cPCA(TCGA_obj, PDXe_obj, global$fast_cPCA)
cov_diff_eig_CCLE <- run_cPCA(CCLE_obj, PDXe_obj, global$fast_cPCA)

if(is.null(global$fast_cPCA)) {
  TCGA_cur_vecs <- cov_diff_eig_TCGA$vectors[, global$remove_cPCA_dims, drop = FALSE]
  CCLE_cur_vecs <- cov_diff_eig_CCLE$vectors[, global$remove_cPCA_dims, drop = FALSE]
} else {
  TCGA_cur_vecs <- cov_diff_eig_TCGA$rotation[, global$remove_cPCA_dims, drop = FALSE]
  CCLE_cur_vecs <- cov_diff_eig_CCLE$rotation[, global$remove_cPCA_dims, drop = FALSE]
}

rownames(TCGA_cur_vecs) <- colnames(TCGA_mat)
rownames(CCLE_cur_vecs) <- colnames(TCGA_mat)

TCGA_cor <- resid(lm(t(TCGA_mat) ~ 0 + TCGA_cur_vecs)) %>% t()
PDXe_cor <- resid(lm(t(PDXe_mat) ~ 0 + TCGA_cur_vecs)) %>% t()
CCLE_cor <- resid(lm(t(CCLE_mat) ~ 0 + CCLE_cur_vecs)) %>% t()

mnn_res_TCGA <- run_MNN(PDXe_cor, TCGA_cor,  k1 = global$mnn_k_tumor, k2 = global$mnn_k_CL, ndist = global$mnn_ndist, 
                   subset_genes = DE_gene_set)
mnn_res_CCLE <- run_MNN(PDXe_cor, CCLE_cor,  k1 = global$mnn_k_tumor, k2 = global$mnn_k_CL, ndist = global$mnn_ndist, 
                        subset_genes = DE_gene_set)
combined_mat <- rbind(mnn_res_TCGA$corrected, PDXe_cor,mnn_res_CCLE$corrected)
comb_ann <- rbind(
  TCGA_ann %>% dplyr::select(sampleID, lineage, subtype, `Primary/Metastasis`) %>%
    dplyr::mutate(type = 'tumor'),
  PDXe_ann %>% dplyr::select(sampleID, lineage, subtype, `Primary/Metastasis`) %>%
    dplyr::mutate(type = 'PDX'),
  CCLE_ann %>% dplyr::select(sampleID, lineage, subtype, `Primary/Metastasis`) %>%
    dplyr::mutate(type = 'CL')
)

comb_obj <- create_Seurat_object(combined_mat, comb_ann)
comb_obj <- cluster_data(comb_obj)
saveRDS(comb_obj,"../Results/Comb_seurat_obj.rds")

alignment<-comb_obj[["umap"]]@cell.embeddings
alignment.df<-data.frame(cbind(alignment[comb_ann$sampleID,],comb_ann))
#alignment.df$type<-as.factor(alignment.df$type)
pdf(file = "../Results/UMAP_embeddings.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 8)

ggplot2::ggplot(data=alignment.df, 
                ggplot2::aes(UMAP_1,UMAP_2,color=factor(type),shape=factor(type)))+
  ggplot2::geom_point(pch=21, alpha=0.7)  +
  #ggplot2::scale_color_manual(values=c(`CL`='black', `tumor`='white')) +
  #ggplot2::scale_size_manual(values=c(`CL`=1, `tumor`=0.75)) +
  ggplot2::theme_classic() + 
  ggplot2::theme(legend.position = 'bottom', 
                 text=ggplot2::element_text(size=8),
                legend.margin =ggplot2::margin(0,0,0,0)) +
  ggplot2::guides(fill=FALSE, color=FALSE) +
  ggplot2::scale_fill_manual(values=tissue_colors) +
  ggplot2::xlab("UMAP 1") +
  ggplot2::ylab("UMAP 2")
dev.off()

  
  