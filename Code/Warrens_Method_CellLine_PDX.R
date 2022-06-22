library(PharmacoGx)
library(reticulate)
conda_list()
use_condaenv(condaenv = 'base', required = TRUE)

pset_folder<-"G:/UHN/psets"
CCLE<-readRDS(paste0(pset_folder,"/CCLE.rds"))
CCLE.exprs<-assay(summarizeMolecularProfiles(CCLE,"Kallisto_0.46.1.rnaseq"))
CCLE.exprs<-CCLE.exprs-log2(0.001)+log2(1)
CCLE.metadata<-cellInfo(CCLE)

PDXe<-readRDS('../Data/PDXE_microArray.rds')
PDX.exprs<-exprs(PDXe)
metadata<-pData(PDXe)

TCGA<-readRDS(paste0(pset_folder,"/TCGA_BRCA.rds"))
TCGA.exprs<-assay(TCGA[["BRCA_RNASeq2GeneNorm-20160128"]])
TCGA.exprs<-log2(TCGA.exprs+1)

