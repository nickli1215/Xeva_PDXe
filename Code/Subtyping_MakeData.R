library(genefu)
library(Biobase)
library(PharmacoGx)
data(pam50.robust)

CCLE<-readRDS("../../psets/CCLE.rds")
CCLE.breast.cells<-rownames(CCLE@cell[CCLE@cell$tissueid=="Breast",])
CCLE.breast<-PharmacoGx::subsetTo(CCLE,CCLE.breast.cells)
saveRDS(CCLE.breast,"../Data/CCLE/CCLE_breast.rds")

CCLE.breast<-readRDS("../Data/CCLE/CCLE_breast.rds")
CCLE.experiment<-summarizeMolecularProfiles(CCLE.breast,mDataType="Kallisto_0.46.1.rnaseq")
CCLE.exprs<-t(assay(CCLE.experiment))
CCLE.exprs<-(CCLE.exprs-rowMeans(CCLE.exprs))/sqrt(rowVars(CCLE.exprs))

CCLE.pdata<-colData(CCLE.experiment)
CCLE.feature.data<-rowData(CCLE.experiment)



colnames(CCLE.exprs) <- CCLE.feature.data[colnames(CCLE.exprs), 'gene_name']
rownames(CCLE.feature.data)<-CCLE.feature.data$gene_name
CCLE.feature.data$EntrezGene.ID<-CCLE.feature.data$gene_id
CCLE.sbt <- molecular.subtyping(sbt.model="pam50", data=CCLE.exprs,
                                 annot=CCLE.feature.data, do.mapping=F)
CCLE.pdata[["subtype"]]<-CCLE.sbt$subtype
CCLE.pdata<-cbind(CCLE.pdata,CCLE.sbt$subtype.proba)
colData(CCLE.experiment)<-CCLE.pdata
CCLE.breast@molecularProfiles[["RNASeq.TPM"]]<-CCLE.experiment
saveRDS(CCLE.breast,"../Data/CCLE/CCLE_breast.rds")

CCLE.TPM<-read.table(file='../Data/CCLE/CCLE_TPM.tsv', sep = '\t', header = TRUE)
TCGA.TPM<-read.table(file='../Data/TCGA/TCGA_BRCA_tpm.tsv',sep='\t',header=TRUE,row.names = 1)
TCGA.IDs<-read.csv(file="../Data/TCGA/TCGA_ID_MAP.csv",header=TRUE,row.names=1)
colnames(TCGA.TPM)<-gsub("X","",colnames(TCGA.TPM))
colnames(TCGA.TPM)<-gsub("[[:punct:]]","-",colnames(TCGA.TPM))

TCGA.BRCA.IDs<-TCGA.IDs[colnames(TCGA.TPM),]
colnames(TCGA.TPM)<-TCGA.BRCA.IDs$AliquotBarcode

TCGA.FPKM<-read.table("../Data/TCGA/TCGA-BRCA.htseq_fpkm.tsv",,sep='\t',header=TRUE,row.names=1)
gene.data<-fData(Morag.data$data)
TCGA.FPKM<-t(TCGA.FPKM)
subset.genes<-colnames(TCGA.FPKM)[colnames(TCGA.FPKM)%in%rownames(gene.data)]
TCGA.FPKM<-TCGA.FPKM[,subset.genes]
colnames(TCGA.FPKM)<-gene.data[subset.genes, 'gene_name']
TCGA.FPKM<-(TCGA.FPKM-rowMeans(TCGA.FPKM))/rowVars(TCGA.FPKM)

TCGA.sbt <- molecular.subtyping(sbt.model="pam50", data=TCGA.FPKM,
                                 annot=gene.data, do.mapping=F)
saveRDS(TCGA.sbt,"../Data/TCGA_TPM_subtyping.rds")


###Morag#####
Morag.data<-readRDS(file='../Data/Morag/pdx_68_samples_rnaseq_gene_exp.rda')
Morag.metadata<-Morag.data$meta
Morag.assay<-Morag.data$data
Morag.exprs<-exprs(Morag.assay)
Morag.pdata<-pData(Morag.assay)
feature.data<-fData(Morag.assay)


Morag.exprs <- t(exprs(Morag.assay)) -log2(0.001)
Morag.exprs<-(Morag.exprs - rowMeans(Morag.exprs))/sqrt(rowVars(Morag.exprs))

colnames(Morag.exprs) <- feature.data[colnames(Morag.exprs), 'gene_name']
Morag.sbt <- molecular.subtyping(sbt.model="pam50", data=Morag.exprs,
                           annot=feature.data, do.mapping=F)
Morag.pdata[["subtype"]]<-Morag.sbt$subtype
Morag.pdata<-cbind(Morag.pdata,Morag.sbt$subtype.proba)
pData(Morag.assay)<-Morag.pdata

saveRDS(Morag.pdata,"../Data/Morag/Morag_phenoData.rds")
saveRDS(Morag.assay,"../Data/Morag/Morag_gene_exp.rds")

#####Cescon#####
load("../Data/Cescon_PDX_gene.exp.kallisto.hg38.RData")
saveRDS(gene.count,"../Data/Cescon/gene_count.rds")
saveRDS(gene.exp,"../Data/Cescon/gene_exp.rds")
saveRDS(transcript.count,"../Data/Cescon/transcript_count.rds")
saveRDS(transcript.exp,"../Data/Cescon/transcript_exp.rds")
Cescon.exprs<-assay(gene.exp)
gene.data<-fData(Morag.data$data)
Cescon.exprs<-t(Cescon.exprs) - log2(0.001)
subset.genes<-colnames(Cescon.exprs)[colnames(Cescon.exprs)%in%rownames(gene.data)]
Cescon.exprs<-Cescon.exprs[,subset.genes]
colnames(Cescon.exprs)<-gene.data[subset.genes, 'gene_name']
Cescon.exprs<-(Cescon.exprs-rowMeans(Cescon.exprs))/sqrt(rowVars(Cescon.exprs))

Cescon.sbt <- molecular.subtyping(sbt.model="pam50", data=Cescon.exprs,
                                annot=gene.data, do.mapping=F)
saveRDS(Cescon.sbt,"../Data/Cescon/subtype.rds")




###PDXe#####
PDXe.data<-readRDS(file='../Data/PDXE_RNASeq_FPKM.rds')
PDXe.exprs<-exprs(PDXe.data)
PDXe.pdata<-pData(PDXe.data)
feature.data<-fData(PDXe.data)


PDXe.exprs <- t(exprs(PDXe.data)) 
PDXe.exprs<-(PDXe.exprs-rowMeans(PDXe.exprs))/sqrt(rowVars(PDXe.exprs))


colnames(PDXe.exprs) <- feature.data[colnames(PDXe.exprs), 'geneName']
PDXe.sbt <- molecular.subtyping(sbt.model="pam50", data=PDXe.exprs,
                                 annot=feature.data, do.mapping=F)
PDXe.pdata[["subtype"]]<-PDXe.sbt$subtype
PDXe.pdata<-cbind(PDXe.pdata,PDXe.sbt$subtype.proba)
pData(PDXe.data)<-PDXe.pdata
PDXe.pdata<-PDXe.pdata[PDXe.pdata$tissue=="BRCA",]
saveRDS(PDXe.pdata,"../Data/PDXe/PDXe_phenoData.rds")
saveRDS(PDXe.data,"../Data/PDXe/PDXe_gene_exp.rds")

PDXe.subset.pdata<-PDXe.pdata[,c("biobase.id","subtype","Basal","Her2","LumA","LumB","Normal")]
colnames(PDXe.subset.pdata)[1]<-"Sample.ID"
PDXe.subset.pdata$dataset<-"PDXe"

Morag.subset.pdata<-Morag.pdata[,c("Sample.ID","subtype","Basal","Her2","LumA","LumB","Normal")]
Morag.subset.pdata$dataset<-"Morag"


Cescon.subset.pdata<-data.frame(cbind(Cescon.sbt$subtype.proba,as.character(Cescon.sbt$subtype)))
colnames(Cescon.subset.pdata)[6]<-"subtype"
Cescon.subset.pdata$Sample.ID<-rownames(Cescon.subset.pdata)
Cescon.subset.pdata<-Cescon.subset.pdata[,c("Sample.ID","subtype","Basal","Her2","LumA","LumB","Normal")]
Cescon.subset.pdata$dataset<-"Cescon"

full.pdata<-data.frame(rbind(Cescon.subset.pdata,Morag.subset.pdata,PDXe.subset.pdata))
write.csv(full.pdata,"../Results/PDX_subtyping.csv")





