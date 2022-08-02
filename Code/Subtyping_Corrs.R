library(PharmacoGx)
library(umap)
library(ggplot2)
library(ComplexHeatmap)
###Morag
Morag.assay<-readRDS("../Data/Morag/Morag_gene_exp.rds")
Morag.rnaseq<-readRDS("../Data/Morag/pdxMorag_rnaseq.rda")
Morag.source<-pData(Morag.rnaseq)[,"Source"]
Morag.tumours<-rownames(pData(Morag.rnaseq))[Morag.source=="Patient"&!is.na(Morag.source)]
Morag.PDXs<-rownames(pData(Morag.rnaseq))[Morag.source=="PDX"&!is.na(Morag.source)]

####Cescon
Cescon.sbt<-readRDS("../Data/Cescon/subtype.rds")
Cescon.gene.exp<-readRDS("../Data/Cescon/gene_exp.rds")
####CCLE
CCLE.breast<-readRDS("../Data/CCLE/CCLE_breast.rds")
CCLE.experiment<-summarizeMolecularProfiles(CCLE.breast,"RNASeq.TPM")
CCLE.cell<-colData(CCLE.experiment)
CCLE.geneinfo<-rowData(CCLE.experiment)
cells<-colSums(is.na(assay(CCLE.experiment)))==0
CCLE.experiment<-CCLE.experiment[,cells]
TCGA.FPKM<-read.table("../Data/TCGA/TCGA-BRCA.htseq_fpkm.tsv",,sep='\t',header=TRUE,row.names=1)
TCGA.subtypes<-readRDS("../Data/TCGA/subtypeinfo.rds")
TCGA.subtypes<-TCGA.subtypes$subtype
gene_subset<-intersect(rownames(CCLE.experiment),rownames(exprs(Morag.assay)))
gene_subset<-intersect(gene_subset,rownames(TCGA.FPKM))
CCLE.experiment<-CCLE.experiment[gene_subset,]
TCGA.FPKM<-TCGA.FPKM[gene_subset,]
Morag.exprs<-exprs(Morag.assay)[gene_subset,]
Cescon.exprs<-exprs(Cescon.gene.exp)[gene_subset,]


CCLE.corrs<-cor(assay(CCLE.experiment),TCGA.FPKM,method="spearman")
Morag.corrs<-cor(Morag.exprs,TCGA.FPKM,method="spearman")
Morag.tumour.corrs<-Morag.corrs[Morag.tumours,]
Morag.PDX.corrs<-Morag.corrs[Morag.PDXs,]
Cescon.corrs<-cor(Cescon.exprs,TCGA.FPKM,method="spearman")

CCLE.subtype.vec<-colData(CCLE.experiment)$subtype[rownames(CCLE.corrs)]
TCGA.subtype.vec<-TCGA.subtypes[colnames(Morag.corrs)]
Cescon.subtype.vec<-Cescon.sbt$subtype

subtypes<-c("Basal","Her2","LumB","LumA","Normal")

Morag.PDX.subtype.vec<-pData(Morag.assay)[Morag.PDXs,"subtype"]
Morag.Tumour.subtype.vec<-pData(Morag.assay)[Morag.PDXs,"subtype"]



Morag.Tumour.subtype.corr.mat<-matrix(,nrow=5,ncol=5)
rownames(Morag.PDX.subtype.corr.mat)<-subtypes
colnames(Morag.Tumour.subtype.corr.mat)<-subtypes


CCLE.subtype.corr.mat<-matrix(,nrow=5,ncol=5)
rownames(CCLE.subtype.corr.mat)<-subtypes
colnames(CCLE.subtype.corr.mat)<-subtypes



getSubtype_Corrs<-function(corrs,col.subtype.vec,row.subtype.vec,subtypes){
  corr.mat<-matrix(,nrow=5,ncol=5)
  rownames(corr.mat)<-subtypes
  colnames(corr.mat)<-subtypes
  for (y.subtype in subtypes){
    sampleCols<-col.subtype.vec[col.subtype.vec==y.subtype]
    for (x.subtype in subtypes){
      sampleRows<-row.subtype.vec[row.subtype.vec==x.subtype]
      subset.corrs<-corrs[sampleRows,sampleCols]
      corr.mat[x.subtype,y.subtype]<-mean(subset.corrs,na.rm=TRUE)
    }
  }
  return(corr.mat)
}


getSubtype_Corrs(Morag.PDX.corrs,TCGA.subtype.vec,Morag.PDX.subtype.vec,subtypes)
getSubtype_Corrs(Morag.tumour.corrs,TCGA.subtype.vec,Morag.Tumour.subtype.vec,subtypes)
getSubtype_Corrs(CCLE.corrs,TCGA.subtype.vec,CCLE.subtype.vec,subtypes)



subtype.mapping<-c(1,2,3,4,5)
##Check HER2, PR and ER and see matching 
##Is the subtype same or not 
##PDXe microarray data for subtyping and see if passaging is maintained
names(subtype.mapping)<-c("Basal","Her2","LumB","LumA","Normal")

TCGA.subtype.numeric.vec<-subtype.mapping[TCGA.subtype.vec]
Morag.subtype.numeric.vec<-subtype.mapping[Morag.subtype.vec]
Morag.TCGA.grid<-meshgrid(Morag.subtype.numeric.vec,TCGA.subtype.numeric.vec)



gene_list<-rownames(assay(CCLE.experiment)[rank(-rowVars(assay(CCLE.experiment)))<=1000,])
CCLE.experiment<-CCLE.experiment[gene_list,]
Morag.exprs<-exprs(Morag.assay)[gene_list,]
TCGA.FPKM<-TCGA.FPKM[gene_list,]
correlations<-cor(assay(CCLE.experiment),Morag.exprs,method="spearman")
morag.subtypes<-pData(Morag.assay)[colnames(correlations),"subtype"]
CCLE.subtypes<-colData(CCLE.experiment)$subtype


Morag.umap<-umap(t(Morag.exprs))

Morag.pData<-pData(Morag.assay)
Morag.umap.df<-cbind(Morag.pData,Morag.umap$layout)
Morag.umap.df<-data.frame(Morag.umap.df)
colnames(Morag.umap.df)[c(66,67)]<-c("umap1","umap2")
ggplot(data=Morag.umap.df,ggplot2::aes(umap1,umap2,color=factor(subtype),shape=factor(subtype)))+
         ggplot2::geom_point(pch=21, alpha=0.7)


