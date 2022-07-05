library(here)
library(magrittr)
library(tidyverse)
library(PharmacoGx)
library(reticulate)
conda_list()
use_condaenv(condaenv = 'base', required = TRUE)
library(Seurat)
source('analysis_helpers.R')
source('global_params.R')
source('Celligner_helpers.R')

library(umap)

###Commments on Expression Data
#TCGA data is RSEM normalized
#

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
Morag.exprs<-Morag.exprs-log2(0.001)+log2(1)
morag.umap <- umap(t(Morag.exprs))
Morag.umap.df<-data.frame(UMAP1=morag.umap$layout[,1],
                          UMAP2=morag.umap$layout[,2],
                          Source=Morag.metadata$Source)
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


TCGA<-readRDS(paste0(pset_folder,"/TCGA/BRCA.rds"))
rownames(TCGA)<-TCGA[,1]
TCGA<-TCGA[,-1]


PDXe<-readRDS("../Data/PDXe_RNASeq_FPKM.rds")
PDXe.metadata<-pData(PDXe)
PDXe.breast.pdxs<-PDXe.metadata[PDXe.metadata$tissue.name=="Breast cancer",]$biobase.id
PDXe.exprs<-exprs(PDXe)[,PDXe.breast.pdxs]


gene_list<-intersect(intersect(intersect(intersect(rownames(CCLE.exprs),rownames(TCGA)),rownames(PDXe.exprs)),rownames(Morag.exprs)),rownames(gCSI.exprs))
TCGA.exprs<-TCGA[gene_list,]
PDXe.exprs<-PDXe.exprs[gene_list,]
CCLE.exprs<-CCLE.exprs[gene_list,]
Morag.exprs<-Morag.exprs[gene_list,]
gCSI.exprs<-gCSI.exprs[gene_list,]

comb_ann <- cbind.data.frame(`sampleID` = c(colnames(TCGA.exprs),colnames(PDXe.exprs),colnames(Morag.exprs),colnames(CCLE.exprs),colnames(gCSI.exprs)),
                             `type` = c(rep('tumor',ncol(TCGA.exprs)),
                                        rep('PDX',ncol(PDXe.exprs)),
                                        as.vector(Morag.metadata[colnames(Morag.exprs),]$Source),
                                        rep("CL",ncol(CCLE.exprs)+ncol(gCSI.exprs))),
                             `source`=c(rep('TCGA',ncol(TCGA.exprs)),
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

all.exprs<-cbind(TCGA.exprs,PDXe.exprs,Morag.exprs,CCLE.exprs,gCSI.exprs)

all.umap = umap(t(all.exprs))
all.umap.df<-data.frame(UMAP1=all.umap$layout[,1],
                          UMAP2=all.umap$layout[,2],
                          Source=comb_ann$type,
                        dataset=comb_ann$source)
ggplot(all.umap.df,aes(x=UMAP1,y=UMAP2,shape=Source,color=Source))+
  geom_point()
ggplot(all.umap.df,aes(x=UMAP1,y=UMAP2,shape=Source,color=dataset))+
  geom_point()

