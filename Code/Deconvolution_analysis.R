library(Biobase)
library(GGally)
library(immunedeconv)
library(tidyr)
library(dplyr)
library(reshape2)
set_cibersort_binary("G:/UHN/CIBERSORT/CIBERSORT.R")
set_cibersort_mat("G:/UHN/CIBERSORT/LM22.txt")

PDX_exprs<-readRDS('../Data/PDXE_microArray.rds')
metadata<-pData(PDX_exprs)
breast_pdxs<-metadata[metadata$tumor.type=="breast",]$biobase.id
exprs_data<-exprs(PDX_exprs)[,breast_pdxs]
exprs_data.not_log<-exp(exprs_data)
deconvolution<-deconvolute(exprs_data.not_log,"cibersort")
breast_metadata<-metadata[breast_pdxs,]
breast_pdx_deconvolution<-rbind(t(breast_metadata)[,breast_pdxs],deconvolution[,-1][,breast_pdxs])
rownames(breast_pdx_deconvolution)<-c(unlist(colnames(breast_metadata)),deconvolution$cell_type)
breast_pdx_deconvolution.df<-as.data.frame(t(breast_pdx_deconvolution[-c(1,4),]))
multiple.samples<-names(table(breast_pdx_deconvolution.df$patient.id)[table(breast_pdx_deconvolution.df$patient.id)>1])
breast_pdx_deconvolution.df<-breast_pdx_deconvolution.df[breast_pdx_deconvolution.df$patient.id %in% multiple.samples,]
breast_pdx_deconvolution.df<-melt(breast_pdx_deconvolution.df,id.vars=c("patient.id","passage"),variable.name="cell.type",value.name="Proportion")

plotting.df<-dcast(breast_pdx_deconvolution.df,patient.id+cell.type ~ passage,function(x) mean(na.omit(as.numeric(x))))
plotting.df[,"cell.type"]<-as.character(plotting.df$cell.type)
rownames(plotting.df)<-plotting.df$patient.id
ggparcoord(plotting.df[,-1],
           columns = 2:6,groupColumn=1
) 

dcast(breast_pdx_deconvolution.df, patient.id+cell.type ~ passage,value.var="Proportion")


ggparcoord(breast_pdx_deconvolution.df,
           columns = 1:4
) 
#biobase.id, patient.id, passage, tumor.type
#Each passage is basically a new PDX, 0th is direct sample
#Subset and remove patients that have only 1 passage/sample
#Subset only to breast cancer
#Look at cibersort results for different passages over time
#natural log