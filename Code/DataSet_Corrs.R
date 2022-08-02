##Check HER2, PR and ER and see matching 
##Is the subtype same or not 
##PDXe microarray data for subtyping and see if passaging is maintained
##Look at other subtyping methods as well (other than pam50)

##grid of violinplots figure with transcriptomic correlations but for datasets
##Remove 0 variance genes and compute common genes
##Separate morag into 2 datasets
##Do this before batch correction 
##Is there a good correlation 
library(PharmacoGx)
library(umap)
library(ggplot2)
library(ComplexHeatmap)
library(tidyr)
library(reshape)
library(genefu)
library(sva)


color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  
  dev.new(width=1.75, height=5)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}
cPalette<-colorRampPalette(c("light green", "yellow", "orange", "red"))(200)

Morag.assay<-readRDS("../Data/Morag/Morag_gene_exp.rds")
Morag.assay<-Morag.assay[rowVars(assay(Morag.assay),na.rm=TRUE)!=0,]

CCLE.breast<-readRDS("../Data/CCLE/CCLE_breast.rds")

CCLE.experiment<-summarizeMolecularProfiles(CCLE.breast,"RNASeq.TPM")
CCLE.cell<-colData(CCLE.experiment)
CCLE.geneinfo<-rowData(CCLE.experiment)
cells<-colSums(is.na(assay(CCLE.experiment)))==0
CCLE.experiment<-CCLE.experiment[rowVars(assay(CCLE.experiment),na.rm=TRUE)!=0,cells]

Cescon.assay<-readRDS("../Data/Cescon/gene_count.rds")
Cescon.exprs<-assay(Cescon.assay)
Cescon.exprs<-Cescon.exprs[rowVars(Cescon.exprs,na.rm=TRUE)!=0,]

TCGA.FPKM<-read.table("../Data/TCGA/TCGA-BRCA.htseq_fpkm.tsv",,sep='\t',header=TRUE,row.names=1)
TCGA.FPKM<-TCGA.FPKM[rowVars(TCGA.FPKM)>0,]
TCGA.subtypes<-readRDS("../Data/TCGA/subtypeinfo.rds")
TCGA.subtypes<-TCGA.subtypes$subtype
gene_subset<-intersect(rownames(CCLE.experiment),rownames(exprs(Morag.assay)))
gene_subset<-intersect(gene_subset,rownames(TCGA.FPKM))
gene_subset<-intersect(gene_subset,rownames(Cescon.exprs))

Morag.rnaseq<-readRDS("../Data/Morag/pdxMorag_rnaseq.rda")
Morag.pdata<-pData(Morag.rnaseq)
Morag.tumours<-rownames(Morag.pdata[Morag.pdata$Source=="Patient",])
Morag.tumours<-Morag.tumours[-22]
Morag.PDXs<-rownames(Morag.pdata[Morag.pdata$Source=="PDX",])
Morag.PDXs<-Morag.PDXs[-28]

Morag.PDX.assay<-Morag.assay[,Morag.PDXs]
Morag.tumours.assay<-Morag.assay[,Morag.tumours]

CCLE.experiment<-CCLE.experiment[gene_subset,]
TCGA.FPKM<-TCGA.FPKM[gene_subset,]

Morag.tumour.exprs<-assay(Morag.tumours.assay)
Morag.tumour.exprs<-Morag.tumour.exprs[gene_subset,]

Morag.PDX.exprs<-assay(Morag.PDX.assay)
Morag.PDX.exprs<-Morag.PDX.exprs[gene_subset,]

Morag.exprs<-exprs(Morag.assay)[gene_subset,]
Cescon.exprs<-Cescon.exprs[gene_subset,]

datasets<-list("Morag Tumours"=Morag.tumour.exprs,
               "Morag PDXs"=Morag.PDX.exprs,
               "TCGA"=TCGA.FPKM,
               "Cescon"=Cescon.exprs,
               "CCLE"=assay(CCLE.experiment))

plots <- matrix(list(), 5, 5)
for (i in seq(1,5)){
  dataset1<-names(datasets)[i]
  data1<-datasets[[dataset1]]
  for (j in seq(1,5)){
    dataset2<-names(datasets)[j]
    if(j==i){
      text=dataset1
      p<-ggplot() + 
        annotate("text",x = 4, y = 25,size=6, label = text) + 
        theme_void()
      plots[[i,j]]<-p
      next
    }
    if(j>i){
      next
      }
    data2<-datasets[[dataset2]]
    corrs<-cor(data1,data2,method="spearman")
    p<-ggplot(melt(corrs),aes(x="",y=value)) + 
      geom_violin(fill = "blue", alpha = 0.3)+ 
      xlab("") + ylab("Correlation")
    plots[[i,j]]<-p
    
    text=round(median(corrs),digits=3)
    colour<-cPalette[round(median(corrs),digits=2)*100]
    p<-ggplot() +  
      annotate("text",x = 4, y = 25,size=8, label = text) + 
      theme_void()
    plots[[j,i]]<-p
    
  }
}
pdf("../Results/Dataset_Comparison.pdf",width=10,height=10)
grid.arrange(arrangeGrob(grobs = plots))
dev.off()
layout <-'
#ABCD
1#EFG
23#HI
456#J
789Z#
'

rowVars(Morag.exprs)+rowVars(Cescon.exprs)+rowVars(assay(CCLE.experiment))+rowVars(as.matrix(TCGA.FPKM))

##########BATCH CORRECTED#########
Morag.data<-readRDS("../Data/Morag/pdx_68_samples_rnaseq_gene_count.rda")
Morag.assay<-Morag.data$data
Morag.assay<-Morag.assay[rowVars(assay(Morag.assay),na.rm=TRUE)!=0,]
Morag.counts<-assay(Morag.assay)*1e6


CCLE.breast<-readRDS("../Data/CCLE/CCLE_breast.rds")

CCLE.experiment<-summarizeMolecularProfiles(CCLE.breast,"Kallisto_0.46.1.rnaseq.counts")
CCLE.cell<-colData(CCLE.experiment)
CCLE.geneinfo<-rowData(CCLE.experiment)
cells<-colSums(is.na(assay(CCLE.experiment)))==0
CCLE.experiment<-CCLE.experiment[rowVars(assay(CCLE.experiment),na.rm=TRUE)!=0,cells]
CCLE.counts<-2**assay(CCLE.experiment) - 1

Cescon.assay<-readRDS("../Data/Cescon/gene_count.rds")
Cescon.counts<-assay(Cescon.assay)*1e6
Cescon.counts<-Cescon.counts[rowVars(Cescon.counts,na.rm=TRUE)!=0,]

TCGA.counts<-read.table("../Data/TCGA/TCGA-BRCA.htseq_counts.tsv",sep='\t',header=TRUE,row.names=1)
TCGA.counts<-2**TCGA.counts[rowVars(as.matrix(TCGA.counts))>0,]-1
TCGA.subtypes<-readRDS("../Data/TCGA/subtypeinfo.rds")
TCGA.subtypes<-TCGA.subtypes$subtype
gene_subset<-intersect(rownames(CCLE.counts),rownames(Morag.counts))
gene_subset<-intersect(gene_subset,rownames(TCGA.counts))
gene_subset<-intersect(gene_subset,rownames(Cescon.counts))

Morag.rnaseq<-readRDS("../Data/Morag/pdxMorag_rnaseq.rda")
Morag.pdata<-pData(Morag.rnaseq)
Morag.tumours<-rownames(Morag.pdata[Morag.pdata$Source=="Patient",])
Morag.tumours<-Morag.tumours[-22]
Morag.PDXs<-rownames(Morag.pdata[Morag.pdata$Source=="PDX",])
Morag.PDXs<-Morag.PDXs[-28]

Morag.PDX.assay<-Morag.assay[,Morag.PDXs]
Morag.Tumour.assay<-Morag.assay[,Morag.tumours]

Morag.Tumour.assay<-Morag.Tumour.assay[gene_subset,]
Morag.PDX.assay<-Morag.PDX.assay[gene_subset,]
TCGA.counts<-TCGA.counts[gene_subset,]
Cescon.counts<-Cescon.counts[gene_subset,]
CCLE.counts<-CCLE.counts[gene_subset,]

all.exprs<-cbind(assay(Morag.Tumour.assay),
      assay(Morag.PDX.assay),
      TCGA.counts,
      Cescon.counts,
      CCLE.counts)
batch <- c(rep("Morag", ncol(assay(Morag.Tumour.assay))),
           rep("Morag", ncol(assay(Morag.PDX.assay))),
           rep("TCGA",ncol(TCGA.counts)),
           rep("Cescon",ncol(Cescon.counts)),
           rep("CCLE",ncol(CCLE.counts)))

adjusted_counts <- ComBat_seq(as.matrix(all.exprs), batch=batch, group=NULL)
gene_lengths<-fData(Morag.rnaseq)[gene_subset,"end"] - fData(Morag.rnaseq)[gene_subset,"start"]
adjusted_counts.FPKM<-convertCounts(
  adjusted_counts,
  "FPKM",
  gene_lengths,
  log = FALSE,
  normalize = "none",
  prior.count = NULL
)

sample.dataset<-c(rep("Morag Tumour", ncol(assay(Morag.Tumour.assay))),
                  rep("Morag PDX", ncol(assay(Morag.PDX.assay))),
                  rep("TCGA",ncol(TCGA.counts)),
                  rep("Cescon",ncol(Cescon.counts)),
                  rep("CCLE",ncol(CCLE.counts)))


plots <- matrix(list(), 5, 5)
for (i in seq(1,5)){
  dataset1<-unique(sample.dataset)[i]
  data1<-all.exprs[,sample.dataset==dataset1]
  for (j in seq(1,5)){
    dataset2<-unique(sample.dataset)[j]
    if(j==i){
      text=dataset1
      p<-ggplot() + 
        annotate("text",x = 4, y = 25,size=6, label = text) + 
        theme_void()
      plots[[i,j]]<-p
      next
    }
    if(j>i){
      next
    }
    data2<-all.exprs[,sample.dataset==dataset2]
    corrs<-cor(data1,data2,method="spearman")
    p<-ggplot(melt(corrs),aes(x="",y=value)) + 
      geom_violin(fill = "blue", alpha = 0.3)+ 
      coord_cartesian(ylim=c(0, 1))+ 
      xlab("") + ylab("Correlation")
    plots[[i,j]]<-p
    
    text=round(median(corrs),digits=3)
    #colour<-cPalette[round(median(corrs),digits=2)*100]
    p<-ggplot() +  
      annotate("text",x = 4, y = 25,size=8, label = text) + 
      theme_void()
    plots[[j,i]]<-p
    
  }
}
pdf("../Results/Dataset_Comparison_Corrected.pdf",width=10,height=10)
grid.arrange(arrangeGrob(grobs = plots))
dev.off()

saveRDS(adjusted_counts,"../Data/Combined_batch_Corrected_gene_exp.rds")
saveRDS(adjusted_counts.FPKM,"../Data/Combined_batchCorrected_FPKM.rds")

dataset.umap<-umap(t(adjusted_counts.FPKM))
dataset.umap.df<-data.frame(dataset.umap$layout)
dataset.umap.df$batch <-batch
dataset.umap.df$dataset<-sample.dataset
corrs<-cor(adjusted_counts.FPKM,method="spearman")
full.cluster<-hclust(as.dist((1-corrs)/2))
Heatmap(corrs,show_row_names = FALSE,show_column_names = FALSE,show_column_dend = FALSE)
dataset.umap.df$cluster<-as.factor(cutree(full.cluster,k=5))

pdf("../Results/Adjusted_Scatterplot.PDF",width=8,height=8)
ggplot(dataset.umap.df, aes(x=X1, y=X2, color=dataset,shape=cluster)) +
  geom_point()
dev.off()

feature.data<-fData(Morag.assay)
colnames(feature.data)[17]<-"Entrez.GeneID"
rownames(adjusted_counts.FPKM) <- feature.data[gene_subset, 'Entrez.GeneID']
all.sbt<- molecular.subtyping(sbt.model="pam50", data=t(log2(adjusted_counts.FPKM+1)),
                                annot=feature.data, do.mapping=F)


