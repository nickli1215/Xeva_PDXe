library(here)
library(magrittr)
library(tidyverse)
library(PharmacoGx)
library(reticulate)
library(msigdbr)
library(umap)
library(Rtsne)
library(tidyr)
library(tibble)
library(dplyr)
###Commments on Expression Data
#TCGA data is log2(RSEM+1) normalized 
#Morag is log2(FPKM+0.001)
#PDXe is log2(FPKM+0.001)
#gCSI and CCLE are log2(TPM+0.001 )
fpkmToTpm <- function(fpkm)
{
  return(exp(log(fpkm) - log(sum(fpkm)) + log(1e6)))
}

fpkmToTpm_matrix <- function(fpkm_matrix) {
  fpkm_matrix_new <- apply(fpkm_matrix, 2, fpkmToTpm)  
  return(fpkm_matrix_new)
}

all_gene_sets<-msigdbr(species="human",category="C5")
immune_related_genes<-all_gene_sets[grepl("immune",all_gene_sets$gs_name,ignore.case=T),]

PDX_morag<-readRDS('../Data/pdxMorag_rnaseq.rda')
Morag.exprs<-exprs(PDX_morag)
Morag.metadata<-pData(PDX_morag)
Morag.gene.metadata<-fData(PDX_morag)
rownames(Morag.gene.metadata)<-Morag.gene.metadata$gene_id
Morag.gene.metadata<-Morag.gene.metadata[rownames(Morag.exprs),]
rownames(Morag.exprs)<-Morag.gene.metadata$gene_name
Morag.metadata<-Morag.metadata[!is.na(Morag.metadata$Source),]
Morag.metadata<-Morag.metadata[intersect(rownames(Morag.metadata),colnames(Morag.exprs)),]
Morag.exprs<-Morag.exprs[,rownames(Morag.metadata)]
Morag.exprs<-2**Morag.exprs - 0.001
Morag.exprs<-fpkmToTpm_matrix(Morag.exprs)


gene_vars<-rowVars(Morag.exprs)
all_genes<-rownames(Morag.exprs)
variance_data<-data.frame(cbind(all_genes,gene_vars,c(1:length(all_genes))))
variance_data<-variance_data %>% group_by(all_genes)%>% slice_max(n = 1, gene_vars, with_ties = FALSE)
indices<-variance_data$V3
Morag.exprs<-Morag.exprs[as.numeric(indices),]
Morag.exprs<-log2(Morag.exprs+1)
means<-rowMeans(t(Morag.exprs))
stdvs<-sqrt(rowVars(t(Morag.exprs)))
Morag.scaled.exprs<-t((t(Morag.exprs)-means)/stdvs)

EMT_genes<-gsub("C1ORD172","C1ORF172",EMT_genes)
EMT_genes<-gsub("ORF","orf",EMT_genes)
EMT_genes<-gsub("TACSTD1","EPCAM",EMT_genes)
EMT_genes<-gsub("GPR110","ADGRF1",EMT_genes)
EMT_genes<-gsub("MTAC2D1","TC2N",EMT_genes)
EMT_genes<-gsub("LRRC54","TSKU",EMT_genes)
EMT_subset_genes<-EMT_genes[EMT_genes%in%rownames(Morag.exprs)]
rownames(Morag.exprs)[grepl( "orf", rownames(Morag.exprs),ignore.case = TRUE, perl = FALSE,
                             fixed = FALSE, useBytes = FALSE)]

Morag.EMT.exprs<-Morag.exprs[EMT_subset_genes,]

Morag.noimmune.exprs<-Morag.exprs[!rownames(Morag.exprs) %in% immune_related_genes$gene_symbol,]
morag.umap <- umap(t(Morag.noimmune.exprs))
Morag.umap.df<-data.frame(UMAP1=morag.umap$layout[,1],
                          UMAP2=morag.umap$layout[,2])
Morag.umap.df<-cbind(Morag.umap.df,Morag.metadata)
colnames(Morag.umap.df)[17]<-"Disease.Pathology"

Morag.metadata[Morag.metadata$Source=="Patient",]$Source<-"tumor"

ggplot(Morag.umap.df,aes(x=UMAP1,y=UMAP2,shape=Source,color=Source))+
  geom_point()


pset_folder<-"G:/UHN/psets"
CCLE<-readRDS(paste0(pset_folder,"/CCLE.rds"))
CCLE.metadata<-cellInfo(CCLE)
breast.cells<-CCLE.metadata[CCLE.metadata$tissueid=="Breast",]$cellid
CCLE<-subsetTo(CCLE,breast.cells)
CCLE.rnaseq<-summarizeMolecularProfiles(CCLE,"Kallisto_0.46.1.rnaseq")
gene.metadata<-CCLE.rnaseq@elementMetadata
CCLE.exprs<-assay(CCLE.rnaseq)
CCLE.exprs<-CCLE.exprs-log2(0.001)+log2(1)
CCLE.exprs<-CCLE.exprs[,colSums(is.na(CCLE.exprs))==0]

rownames(gene.metadata)<-gene.metadata$gene_id
gene.metadata<-gene.metadata[rownames(CCLE.exprs),]
CCLE.exprs<-CCLE.exprs[rownames(gene.metadata),]
rownames(CCLE.exprs)<-gene.metadata$gene_name
gene_vars<-rowVars(CCLE.exprs)
all_genes<-rownames(CCLE.exprs)
variance_data<-data.frame(cbind(all_genes,gene_vars,c(1:length(all_genes))))
variance_data<-variance_data %>% group_by(all_genes)%>% slice_max(n = 1, gene_vars, with_ties = FALSE)
indices<-variance_data$V3
CCLE.exprs<-CCLE.exprs[as.numeric(indices),]
#means<-rowMeans(t(CCLE.exprs))
#stdvs<-sqrt(rowVars(t(CCLE.exprs)))
#CCLE.exprs<-t((t(CCLE.exprs)-means)/stdvs)


gCSI<-readRDS(paste0(pset_folder,"/gCSI2.rds"))
gCSI.metadata<-cellInfo(gCSI)
breast.cells<-gCSI.metadata[gCSI.metadata$tissueid=="Breast",]$cellid
gCSI<-subsetTo(gCSI,breast.cells)
gCSI.rnaseq<-summarizeMolecularProfiles(gCSI,"Kallisto_0.46.1.rnaseq")
gene.metadata<-gCSI.rnaseq@elementMetadata
gCSI.exprs<-assay(gCSI.rnaseq)
gCSI.exprs<-gCSI.exprs-log2(0.001)+log2(1)
gCSI.exprs<-gCSI.exprs[,colSums(is.na(gCSI.exprs))==0]

rownames(gene.metadata)<-gene.metadata$gene_id
gene.metadata<-gene.metadata[rownames(gCSI.exprs),]
gCSI.exprs<-gCSI.exprs[rownames(gene.metadata),]
rownames(gCSI.exprs)<-gene.metadata$gene_name
gene_vars<-rowVars(gCSI.exprs)
all_genes<-rownames(gCSI.exprs)
variance_data<-data.frame(cbind(all_genes,gene_vars,c(1:length(all_genes))))
variance_data<-variance_data %>% group_by(all_genes)%>% slice_max(n = 1, gene_vars, with_ties = FALSE)
indices<-variance_data$V3
gCSI.exprs<-gCSI.exprs[as.numeric(indices),]
#means<-rowMeans(t(gCSI.exprs))
#stdvs<-sqrt(rowVars(t(gCSI.exprs)))
#gCSI.exprs<-t((t(gCSI.exprs)-means)/stdvs)


TCGA<-readRDS(paste0(pset_folder,"/TCGA/BRCA.rds"))
rownames(TCGA)<-TCGA[,1]
TCGA<-TCGA[,-1]
TCGA<-TCGA[as.numeric(indices),]
#means<-rowMeans(t(TCGA))
#stdvs<-sqrt(rowVars(t(TCGA)))
#TCGA<-t((t(TCGA)-means)/stdvs)

PDXe<-readRDS("../Data/PDXe_RNASeq_FPKM.rds")
PDXe.metadata<-pData(PDXe)
PDXe.breast.pdxs<-PDXe.metadata[PDXe.metadata$tissue.name=="Breast cancer",]$biobase.id
PDXe.exprs<-log2(fpkmToTpm_matrix(2**exprs(PDXe)-1)[,PDXe.breast.pdxs]+1)
gene_vars<-rowVars(PDXe.exprs)
all_genes<-rownames(PDXe.exprs)
variance_data<-data.frame(cbind(all_genes,gene_vars,c(1:length(all_genes))))
variance_data<-variance_data %>% group_by(all_genes)%>% slice_max(n = 1, gene_vars, with_ties = FALSE)
indices<-variance_data$V3
PDXe.exprs<-PDXe.exprs[as.numeric(indices),]
#means<-rowMeans(t(PDXe.exprs))
#stdvs<-sqrt(rowVars(t(PDXe.exprs)))
#PDXe.exprs<-t((t(PDXe.exprs)-means)/stdvs)


gene_list<-intersect(intersect(intersect(intersect(rownames(CCLE.exprs),rownames(TCGA)),rownames(PDXe.exprs)),rownames(Morag.exprs)),rownames(gCSI.exprs))
#Morag.exprs<-Morag.exprs[Morag.gene.metadata$gene_name,]

Morag.gene.metadata<-Morag.gene.metadata[Morag.gene.metadata$gene_type=="protein_coding",]

TCGA.exprs<-TCGA[gene_list,]
PDXe.exprs<-PDXe.exprs[gene_list,]
CCLE.exprs<-CCLE.exprs[gene_list,]
Morag.exprs<-Morag.exprs[gene_list,]
gCSI.exprs<-gCSI.exprs[gene_list,]

comb_ann <- cbind.data.frame(`sampleID` = c(colnames(PDXe.exprs),colnames(Morag.exprs),colnames(CCLE.exprs),colnames(gCSI.exprs)),
                             `type` = c(#rep('tumor',ncol(TCGA.exprs)),
                                        rep('PDX',ncol(PDXe.exprs)),
                                        as.vector(Morag.metadata[colnames(Morag.exprs),]$Source),
                                        rep("CL",ncol(CCLE.exprs)+ncol(gCSI.exprs))),
                             `source`=c(#rep('TCGA',ncol(TCGA.exprs)),
                                        rep('PDXe',ncol(PDXe.exprs)),
                                        rep('Morag',ncol(Morag.exprs)),
                                        rep("CCLE",ncol(CCLE.exprs)),
                                        rep("gCSI",ncol(gCSI.exprs))))
rownames(comb_ann)<-comb_ann$sampleID
original_combined_obj <-  CreateSeuratObject(cbind(TCGA.exprs,PDX.exprs,CCLE.exprs),
                                               min.cells = 0,
                                               min.features = 0,
                                               meta.data = comb_ann %>%
                                                 magrittr::set_rownames(comb_ann$sampleID))

#all.exprs<-cbind(TCGA.exprs,PDXe.exprs,Morag.exprs,CCLE.exprs,gCSI.exprs)
all.exprs<-cbind(PDXe.exprs,Morag.exprs,CCLE.exprs,gCSI.exprs)

all.umap = umap(t(all.exprs))
all.umap.df<-data.frame(UMAP1=all.umap$layout[,1],
                          UMAP2=all.umap$layout[,2],
                          Source=comb_ann$type,
                        dataset=comb_ann$source)
ggplot(all.umap.df,aes(x=UMAP1,y=UMAP2,shape=Source,color=Source))+
  geom_point()
ggplot(all.umap.df,aes(x=UMAP1,y=UMAP2,shape=Source,color=dataset))+
  geom_point()




