library(Biobase)
library(GGally)
library(immunedeconv)
library(tidyr)
library(tibble)
library(dplyr)
library(reshape2)
library(stringr)
library(ggpubr)
library(gridExtra)
library(PharmacoGx)
library(grDevices)
library(ComplexHeatmap)

fpkmToTpm <- function(fpkm)
{
  return(exp(log(fpkm) - log(sum(fpkm)) + log(1e6)))
}

fpkmToTpm_matrix <- function(fpkm_matrix) {
  fpkm_matrix_new <- apply(fpkm_matrix, 2, fpkmToTpm)  
  return(fpkm_matrix_new)
}

set_cibersort_binary("G:/UHN/CIBERSORT/CIBERSORT.R")
set_cibersort_mat("G:/UHN/CIBERSORT/LM22.txt")
pset_folder<-"G:/UHN/psets"
CCLE<-readRDS(paste0(pset_folder,"/CCLE.rds"))
CCLE.breast.cells<-rownames(CCLE@cell[CCLE@cell$tissueid=="Breast",])
CCLE<-subsetTo(CCLE,cells=CCLE.breast.cells)
CCLE.experiment<-summarizeMolecularProfiles(CCLE,mDataType="Kallisto_0.46.1.rnaseq")
#CCLE.experiment<-summarizeMolecularProfiles(CCLE,mDataType="rna")

gene_metadata<-CCLE.experiment@elementMetadata
#####CCLE##########
CCLE.exprs<-assay(CCLE.experiment)
CCLE.exprs<-CCLE.exprs[gene_metadata$gene_id,]
rownames(CCLE.exprs)<-gene_metadata$gene_name
CCLE.exprs<-CCLE.exprs[gene_metadata$gene_type=="protein_coding",]
CCLE.exprs<-CCLE.exprs[,colSums(is.na(CCLE.exprs))==0]
CCLE.exprs<-2**CCLE.exprs - 0.001 
gene_vars<-rowVars(CCLE.exprs)
all_genes<-rownames(CCLE.exprs)
variance_data<-data.frame(cbind(all_genes,gene_vars,c(1:length(all_genes))))
variance_data<-variance_data %>% group_by(all_genes)%>% slice_max(n = 1, gene_vars, with_ties = FALSE)
indices<-variance_data$V3
CCLE.exprs<-CCLE.exprs[as.numeric(indices),]

#####GDSC##########
GDSC<-readRDS(paste0(pset_folder,"/GDSC2.rds"))
gdsc.breast.cells<-rownames(GDSC@cell[GDSC@cell$tissueid=="Breast",])
GDSC<-subsetTo(GDSC,cells=gdsc.breast.cells)
GDSC.experiment<-summarizeMolecularProfiles(GDSC,mDataType="Kallisto_0.46.1.rnaseq")
GDSC.exprs<-assay(GDSC.experiment)
GDSC.exprs<-GDSC.exprs[gene_metadata$gene_id,]
rownames(GDSC.exprs)<-gene_metadata$gene_name
GDSC.exprs<-GDSC.exprs[gene_metadata$gene_type=="protein_coding",]
GDSC.exprs<-GDSC.exprs[,colSums(is.na(GDSC.exprs))==0]
GDSC.exprs<-2**GDSC.exprs - 0.001 
gene_vars<-rowVars(GDSC.exprs)
all_genes<-rownames(GDSC.exprs)
variance_data<-data.frame(cbind(all_genes,gene_vars,c(1:length(all_genes))))
variance_data<-variance_data %>% group_by(all_genes)%>% slice_max(n = 1, gene_vars, with_ties = FALSE)
indices<-variance_data$V3
GDSC.exprs<-GDSC.exprs[as.numeric(indices),]


#####gCSI##########
gCSI<-readRDS(paste0(pset_folder,"/gCSI2.rds"))
gCSI.breast.cells<-rownames(gCSI@cell[gCSI@cell$tissueid=="Breast",])
gCSI<-subsetTo(gCSI,cells=gCSI.breast.cells)
gCSI.experiment<-summarizeMolecularProfiles(gCSI,mDataType="Kallisto_0.46.1.rnaseq")
gCSI.exprs<-assay(gCSI.experiment)

gCSI.exprs<-gCSI.exprs[gene_metadata$gene_id,]
rownames(gCSI.exprs)<-gene_metadata$gene_name
gCSI.exprs<-gCSI.exprs[gene_metadata$gene_type=="protein_coding",]
gCSI.exprs<-gCSI.exprs[,colSums(is.na(gCSI.exprs))==0]
gCSI.exprs<-2**gCSI.exprs - 0.001 
gene_vars<-rowVars(gCSI.exprs)
all_genes<-rownames(gCSI.exprs)
variance_data<-data.frame(cbind(all_genes,gene_vars,c(1:length(all_genes))))
variance_data<-variance_data %>% group_by(all_genes)%>% slice_max(n = 1, gene_vars, with_ties = FALSE)
indices<-variance_data$V3
gCSI.exprs<-gCSI.exprs[as.numeric(indices),]


#####GRAY##########
GRAY<-readRDS(paste0(pset_folder,"/GRAY2017.rds"))
GRAY.breast.cells<-rownames(GRAY@cell[GRAY@cell$tissueid=="Breast",])
GRAY<-subsetTo(GRAY,cells=GRAY.breast.cells)
GRAY.experiment<-summarizeMolecularProfiles(GRAY,mDataType="Kallisto_0.46.1.rnaseq")
GRAY.exprs<-assay(GRAY.experiment)
gene_metadata<-GRAY.experiment@elementMetadata

GRAY.exprs<-GRAY.exprs[gene_metadata$gene_id,]
rownames(GRAY.exprs)<-gene_metadata$gene_name
GRAY.exprs<-GRAY.exprs[gene_metadata$gene_type=="protein_coding",]
GRAY.exprs<-GRAY.exprs[,colSums(is.na(GRAY.exprs))==0]
GRAY.exprs<-2**GRAY.exprs - 0.001 
gene_vars<-rowVars(GRAY.exprs)
all_genes<-rownames(GRAY.exprs)
variance_data<-data.frame(cbind(all_genes,gene_vars,c(1:length(all_genes))))
variance_data<-variance_data %>% group_by(all_genes)%>% slice_max(n = 1, gene_vars, with_ties = FALSE)
indices<-variance_data$V3
GRAY.exprs<-GRAY.exprs[as.numeric(indices),]


#####PDXe##########
PDXe<-readRDS("../Data/PDXe_RNASeq_FPKM.rds")
metadata<-pData(PDXe)
tissue.types<-rownames(table(metadata$tissue)[table(metadata$tissue)>15])
full_pdxs<-metadata[metadata$tissue %in% tissue.types,]$biobase.id
exprs_data<-exprs(PDXe)
exprs_data.not_log<-2**exprs_data -1
exprs_data.not_log.tpm<-fpkmToTpm_matrix(exprs_data.not_log)
breast_pdxs<-metadata[metadata$tissue=="BRCA",]$biobase.id
PDXe.exprs<-exprs_data.not_log.tpm[,breast_pdxs]
gene_vars<-rowVars(PDXe.exprs)
all_genes<-rownames(PDXe.exprs)
variance_data<-data.frame(cbind(all_genes,gene_vars,c(1:length(all_genes))))
variance_data<-variance_data %>% group_by(all_genes)%>% slice_max(n = 1, gene_vars, with_ties = FALSE)
indices<-variance_data$V3
PDXe.exprs<-PDXe.exprs[as.numeric(indices),]

#####TCGA##########
rawdata_folder<-"G:/UHN/RawData"
TCGA<-readRDS(paste0(pset_folder,"/TCGA/GDCData/BRCA.rds"))
rownames(TCGA)<-TCGA[,1]
TCGA.exprs<-2**TCGA[,-1]-1
gene_vars<-rowVars(TCGA.exprs)
all_genes<-rownames(TCGA.exprs)
variance_data<-data.frame(cbind(all_genes,gene_vars,c(1:length(all_genes))))
variance_data<-variance_data %>% group_by(all_genes)%>% slice_max(n = 1, gene_vars, with_ties = FALSE)
indices<-variance_data$V3
TCGA.exprs<-TCGA.exprs[as.numeric(indices),]

#####Morag##########
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
Morag.PDX.exprs<-Morag.exprs[,rownames(Morag.metadata[Morag.metadata$Source=="PDX",])]
Morag.Tumor.exprs<-Morag.exprs[,rownames(Morag.metadata[Morag.metadata$Source=="Patient",])]