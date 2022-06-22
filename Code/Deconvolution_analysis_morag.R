library(Biobase)
library(GGally)
library(immunedeconv)
library(tidyr)
library(tibble)
library(dplyr)
library(reshape2)
library(stringr) 
library(ggpubr)
library(PharmacoGx)
library(gridExtra)
library(grDevices)
library(ComplexHeatmap)

set_cibersort_binary("G:/UHN/CIBERSORT/CIBERSORT.R")
set_cibersort_mat("G:/UHN/CIBERSORT/LM22.txt")
pset_folder<-"G:/UHN/psets"


plot_cellwise_results<-function(deconvolution_df,deconvolution_method,results_dir){
  all_cell_boxplots<-list()
  
  for (cell in unique(deconvolution_df$cell.type)){
    cell_data<-deconvolution_df[deconvolution_df$cell.type==cell,]
    cell_data$Source<-as.factor(cell_data$Source)
    cell_data$Proportion<-as.numeric(as.vector(cell_data$Proportion))
    cell_data<-cell_data[!is.na(cell_data$Source),]
    
    mean_duplicate<-mean(cell_data[cell_data$patient.id=="V2013-1738"&cell_data$Source=="PDX","Proportion"])
    cell_data[cell_data$patient.id=="V2013-1738"&cell_data$Source=="PDX","Proportion"]=mean_duplicate
    cell_data<-cell_data %>% distinct(patient.id,Source,.keep_all=TRUE)
    
    patient_pdxs<-cell_data[cell_data$Source=="Patient",]
    rownames(patient_pdxs)<-patient_pdxs$patient.id
    derived_pdxs<-cell_data[cell_data$Source=="PDX",]
    rownames(derived_pdxs)<-derived_pdxs$patient.id
    
    matched_patients<-intersect(rownames(derived_pdxs),rownames(patient_pdxs))
    patient_pdxs<-patient_pdxs[matched_patients,]
    derived_pdxs<-derived_pdxs[matched_patients,]
    
    my_comparisons<-list("Patient"=patient_pdxs$Proportion,"Source"=derived_pdxs$Proportion)
    bplot<-ggboxplot(cell_data, x = "Source", y = "Proportion",
                 color="Source",id="patient.id",
                 palette = "jco",add = "jitter")+
      font("title",size=12)+ggtitle(cell)+
      #stat_compare_means(label = "p.format",comparisons=my_comparisons,
      #method="wilcox.test",paired=TRUE,label.y=max(cell_data$Proportion))+ theme_bw()
      stat_compare_means(method="wilcox.test")
    
    
    all_cell_boxplots[[cell]]<-bplot
    
    heatmap_data<-as.matrix(cbind(patient_pdxs$Proportion,derived_pdxs$Proportion))
    colnames(heatmap_data)<-c("Patient","PDX")
    rownames(heatmap_data)<-patient_pdxs$patient.id
    
    hmap<-Heatmap(heatmap_data,cluster_rows=FALSE,
                     column_order=c("Patient","PDX"),
                     na_col = "black",
                     show_column_dend = FALSE,
                     column_title = cell)
    if(grepl( "/",cell, fixed = TRUE)){
      cell<-str_replace_all(cell,"/","+")
    }
    mainDir<-paste0(results_dir,"/",deconvolution_method)
    dir.create(mainDir, showWarnings = FALSE)
    dir.create(paste0(mainDir,"/heatmaps"),showWarnings = FALSE)
    pdf(paste0(mainDir,"/","heatmaps/",cell,"_heatmap.pdf"),height=10,width=5)
    print(hmap)
    dev.off()
    
  }
  ml<-ggarrange(plotlist=all_cell_boxplots,nrow=3,ncol=2,common.legend=TRUE)
  
  mainDir<-paste0(results_dir,"/",deconvolution_method)
  dir.create(mainDir, showWarnings = FALSE)
  pdf(paste0(mainDir,"/",deconvolution_method,"_Full_boxplot.pdf"),height=9,width=6)
  print(ml)
  dev.off() 
  
}

plot_patientwise_results<-function(deconvolution_df,deconvolution_method,results_dir){
  patient_pdxs<-deconvolution_df[deconvolution_df$Source=="Patient",]
  derived_pdxs<-deconvolution_df[deconvolution_df$Source=="PDX",]
  matched_patients<-intersect(derived_pdxs$patient.id,patient_pdxs$patient.id)
  for (patient in matched_patients){
    patient_data<-deconvolution_df[deconvolution_df$patient.id==patient,]
    patient_data$Proportion<-as.numeric(patient_data$Proportion)
    p<-ggplot(patient_data,aes(fill=cell.type,x=Source,y=Proportion))+ 
      geom_bar(position="fill", stat="identity")+
      geom_col(width=1)
    
    mainDir<-paste0(results_dir,"/",deconvolution_method)
    dir.create(mainDir, showWarnings = FALSE)
    dir.create(paste0(mainDir,"/barplots"),showWarnings = FALSE)
    pdf(paste0(mainDir,"/","barplots/",patient,"_barplot.pdf"),height=10,width=5)
    print(p)
    dev.off()
    
  }
}


get_PDX_results<-function(dat_exprs,pdxs,metadata,method,results_dir){
  dir.create(results_dir, showWarnings = FALSE)
  if(method=="consensus_tme"){
    deconvolution<-deconvolute(dat_exprs[,pdxs],method,indications=rep("brca",length(pdxs)))
  }else{
    deconvolution<-deconvolute(dat_exprs[,pdxs],method)
  }
  metadata<-metadata[pdxs,]
  
  if(method=="mcp_counter"){
    #deconvolution[-1]<-deconvolution[-1]/colSums(deconvolution[-1])
  }
  breast_deconvolution<-rbind(t(metadata)[c("Anonymous.Patient.ID","Source"),pdxs],deconvolution[,-1][,pdxs])
  rownames(breast_deconvolution)<-c("patient.id","Source",deconvolution$cell_type)
  breast_deconvolution.df<-as.data.frame(t(breast_deconvolution))
  #multiple.samples<-names(table(breast_deconvolution.df$patient.id)[table(breast_deconvolution.df$patient.id)>1])
  #breast_deconvolution.df<-breast_deconvolution.df[breast_deconvolution.df$patient.id %in% multiple.samples,]
  breast_deconvolution.df<-melt(breast_deconvolution.df,id.vars=c("patient.id","Source"),variable.name="cell.type",value.name="Proportion")
  breast_deconvolution.df<-breast_deconvolution.df[breast_deconvolution.df$patient.id!=""&
                                                     !is.na(breast_deconvolution.df$patient.id),]
  print("Plotting results")
  plot_cellwise_results(breast_deconvolution.df,method,results_dir)
  plot_patientwise_results(breast_deconvolution.df,method,results_dir)
  #wilcox_test_results(breast_deconvolution.df,method)
  
  #return(breast_deconvolution.df)
}


PDX_exprs<-readRDS('../Data/pdxMorag_rnaseq.rda')
metadata<-pData(PDX_exprs)
gene_data<-fData(PDX_exprs)
exprs_data<-exprs(PDX_exprs)
pdxs<-intersect(metadata$Sample.ID,colnames(exprs_data))
exprs_data<-exprs_data[gene_data$gene_id,]
rownames(exprs_data)<-gene_data$Symbol
exprs_data<-exprs_data[gene_data$gene_type=="protein_coding",]


gene_vars<-rowVars(exprs_data)
all_genes<-rownames(exprs_data)
variance_data<-data.frame(cbind(all_genes,gene_vars,c(1:length(all_genes))))
variance_data<-variance_data %>% group_by(all_genes)%>% slice_max(n = 1, gene_vars, with_ties = FALSE)
indices<-variance_data$V3
exprs_data<-exprs_data[as.numeric(indices),]



exprs_data.not_log<-2**exprs_data


get_PDX_results(exprs_data.not_log,pdxs,metadata,"cibersort","../Results/PDX_morag")
get_PDX_results(exprs_data.not_log,pdxs,metadata,"cibersort_abs","../Results/PDX_morag")
get_PDX_results(exprs_data.not_log,pdxs,metadata,"xcell","../Results/PDX_morag")
get_PDX_results(exprs_data.not_log,pdxs,metadata,"quantiseq","../Results/PDX_morag")
get_PDX_results(exprs_data.not_log,pdxs,metadata,"mcp_counter","../Results/PDX_morag") #Convert to percentages
get_PDX_results(exprs_data.not_log,pdxs,metadata,"epic","../Results/PDX_morag")
get_PDX_results(exprs_data.not_log,pdxs,metadata,"consensus_tme","../Results/PDX_morag")