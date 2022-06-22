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
library(gridExtra)
library(grDevices)
library(ComplexHeatmap)


set_cibersort_binary("G:/UHN/CIBERSORT/CIBERSORT.R")
set_cibersort_mat("G:/UHN/CIBERSORT/LM22.txt")
pset_folder<-"G:/UHN/psets"

PDX_exprs<-readRDS('../Data/PDXE_microArray.rds')
metadata<-pData(PDX_exprs)
tissue.types<-rownames(table(metadata$tumor.type)[table(metadata$tumor.type)>15])
full_pdxs<-metadata[metadata$tumor.type %in% tissue.types,]$biobase.id

breast_pdxs<-metadata[metadata$tumor.type=="breast",]$biobase.id
exprs_data<-exprs(PDX_exprs)

breast_exprs_data<-exprs_data[,breast_pdxs]
gene_vars<-rowVars(breast_exprs_data)
all_genes<-rownames(breast_exprs_data)
variance_data<-data.frame(cbind(all_genes,gene_vars,c(1:length(all_genes))))
variance_data<-variance_data %>% group_by(all_genes)%>% slice_max(n = 1, gene_vars, with_ties = FALSE)
indices<-variance_data$V3
breast_exprs_data<-breast_exprs_data[as.numeric(indices),]
breast_exprs_data.not_log<-2**breast_exprs_data


exprs_data<-exprs(PDX_exprs)
gene_vars<-rowVars(exprs_data)
all_genes<-rownames(exprs_data)
variance_data<-data.frame(cbind(all_genes,gene_vars,c(1:length(all_genes))))
variance_data<-variance_data %>% group_by(all_genes)%>% slice_max(n = 1, gene_vars, with_ties = FALSE)
indices<-variance_data$V3
exprs_data<-exprs_data[as.numeric(indices),]
exprs_data.not_log<-2**exprs_data

breast_metadata<-metadata[breast_pdxs,]
breast_metadata$passage<-ifelse(breast_metadata$passage>=1,"PDX","Patient")

metadata$passage<-as.vector(symnum(metadata$passage,cutpoints=c(0,1.1,2.1,6),symbols=c("Patient","1-2","3+")))



get_pvalue<-function(x,data){
  group1<-x[1]
  group2<-x[2]
  
  group_1_data<-data[data$passage==group1,]
  group_1_data<-aggregate(x=group_1_data$Proportion, by = list(group_1_data$patient.id), FUN = mean)
  colnames(group_1_data)<-c("patient.id",'Proportion')
  group_2_data<-data[data$passage==group2,]
  group_2_data<-aggregate(x=group_2_data$Proportion, by = list(group_2_data$patient.id), FUN = mean)
  colnames(group_2_data)<-c("patient.id",'Proportion')
  
  rownames(group_1_data)<-group_1_data$patient.id
  rownames(group_2_data)<-group_2_data$patient.id
  
  patients<-intersect(group_1_data$patient.id,group_2_data$patient.id)
  
  group_1_data<-group_1_data[patients,]
  group_2_data<-group_2_data[patients,]
  
  wilc.test<-wilcox.test(group_1_data$Proportion,group_2_data$Proportion,paired=TRUE)
  return(wilc.test$p.value)
}

plot_results<-function(deconvolution_df,deconvolution_method,results_dir,xlab="Passage no."){
  all_cell_plots<-list()
  
  for (cell in unique(deconvolution_df$cell.type)){
    cell_data<-deconvolution_df[deconvolution_df$cell.type==cell,]
    cell_data$passage<-as.factor(cell_data$passage)
    cell_data$Proportion<-as.numeric(as.vector(cell_data$Proportion))
    cell_data_long<-melt(cell_data[,-3],id=c("patient.id","passage"),value.name = "Proportion")
    cell_data_long<-cell_data_long[!is.na(cell_data_long$Proportion),]
    ylim1 = boxplot.stats(cell_data_long$Proportion)$stats[c(1, 5)]
    
    comparisons<-list()
    passages<-str_sort(unique(cell_data_long$passage))
    print(passages)
    if("Patient" %in% passages){
      passages<-passages[passages!="Patient"]
      passages<-append("Patient",passages)
    }
    
    base<-as.character(passages[1])
    compared_passages<-passages[-1]
    for(i in 1:length(compared_passages)){
      comparisons[[i]]<-c(base,as.character(compared_passages[i]))
    }

    groups<-do.call(rbind, comparisons)
    
    pvalues<-lapply(comparisons,function(x){return(get_pvalue(x,cell_data_long))})
    pvalues<-round(unlist(pvalues),5)
    .y.<-rep("Proportion",length(pvalues))
    if(deconvolution_method=="mcp_counter"){
      .y.<-rep("Score",length(pvalues))
    }
    p.signif<-symnum(pvalues,cutpoints=c(0,0.0001,0.001,0.01,0.05,1),symbols=c("****","***","**","*","ns"))
    p.signif<-as.vector(p.signif)
    
    
    stats<-data.frame(.y.=.y.,group1=as.vector(groups[,1]),
                      group2=as.vector(groups[,2]),
                      p=unlist(pvalues),p.signif=p.signif)
    stats<-as_tibble(stats)
    #ggplot(data=cell_data, aes(x=passage, y=Proportion, group=patient.id)) +
      #geom_line() + geom_point()+
    #ggplot(data=cell_data, aes(x=passage, y=Proportion,)) + 
      #geom_violin() +
      #scale_color_brewer(palette="Paired")+
     # theme_minimal() + ggtitle(paste0(deconvolution_method,": ",cell))
    #p<-ggpaired(cell_data, x = "passage", y = "Proportion",color="passage",id="patient.id")+
    print(passages)
    print(as.vector(passages))
    cell_data_long$passage<-factor(cell_data_long$passage,levels=as.vector(passages))
    p0<-ggboxplot(cell_data_long, x = "passage", y = "Proportion",
                  color="passage",id="patient.id",
                  palette = "jco")+
      font("title",size=12)+ggtitle(cell)+
      #stat_compare_means(aes(group=Dataset),method="wilcox.test",comparisons=my_comparisons,paired=TRUE)+
      stat_pvalue_manual(data=stats,label="p.signif",x="group2", y.position = ylim1[2])+
      theme_bw()
    
    if(deconvolution_method %in% c("mcp_counter","consensus_tme")){
      cell_data_long$variable<-"Score"
      colnames(cell_data_long)[4]<-"Score"
      p0<-ggboxplot(cell_data_long, x = "passage", y = "Score",
                    color="passage",id="patient.id",
                    palette = "jco")+
        font("title",size=12)+ggtitle(cell)+
        #stat_compare_means(aes(group=Dataset),method="wilcox.test",comparisons=my_comparisons,paired=TRUE)+
        stat_pvalue_manual(data=stats,label="p.signif",x="group2", y.position = ylim1[2])+
        theme_bw()
    }    
    all_cell_plots[[cell]]<-p0
    
    #cell_data<-aggregate(cell_data$Proportion, by=list(cell_data$passage,cell_data$patient.id),
                       #     FUN=mean, na.rm=TRUE)
    cell_data_aggregate<-aggregate(x=cell_data$Proportion,by=list(cell_data$patient.id,cell_data$passage), FUN = mean)
    colnames(cell_data_aggregate)<-c("patient.id","passage","Proportion")
    heatmap_data<-reshape(cell_data_aggregate,idvar="patient.id",timevar = "passage",direction="wide")
    rownames(heatmap_data)<-heatmap_data[,1]
    heatmap_data<-heatmap_data[,-1]
    colnames(heatmap_data)<-str_remove_all(colnames(heatmap_data),"Proportion.")
    
    ordered_columns<-passages
    heatmap_data<-heatmap_data[rowSums(!is.na(heatmap_data))>1,]
    
    hmap<-Heatmap(as.matrix(heatmap_data),cluster_rows=FALSE,
                  column_order=ordered_columns,
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
  #ml <- marrangeGrob(all_cell_plots, nrow=2, ncol=2, top=deconvolution_method)
  ml<-ggarrange(plotlist=all_cell_plots,nrow=3,ncol=2,common.legend=TRUE)
  
  mainDir<-paste0(results_dir,"/",deconvolution_method)
  dir.create(mainDir, showWarnings = FALSE)
  pdf(paste0(mainDir,"/",deconvolution_method,"_Full_boxplot.pdf"),height=9,width=6)
  print(ml)
  dev.off()  
}
wilcox_test<-function(cell_data,combination){
  i<-combination[1]
  j<-combination[2]
  passage1_proportions<-cell_data[(cell_data$passage==i)&(!is.na(cell_data$Proportion)),]$Proportion 
  passage2_proportions<-cell_data[(cell_data$passage==j)&(!is.na(cell_data$Proportion)),]$Proportion 
  if(length(passage1_proportions)<1|length(passage2_proportions)<1){
    return(NA)
  }
  res<-wilcox.test(passage1_proportions, passage2_proportions)
  return(res$p.value)
}
wilcox_test_results<-function(deconvolution_df,deconvolution_method,results_dir){
  for (cell in unique(deconvolution_df$cell.type)){
    print(cell)
    cell_data<-deconvolution_df[deconvolution_df$cell.type==cell,]
    cell_data$passage<-as.factor(cell_data$passage)
    cell_data$Proportion<-as.numeric(as.vector(cell_data$Proportion))
    combinations<-t(as.matrix(combn(6,2)))
    res<-apply(combinations,1,function(x) wilcox_test(cell_data,x))
    res<-cbind(combinations,as.vector(unlist(res)))
    mainDir<-paste0(results_dir,"/",deconvolution_method)
    dir.create(mainDir, showWarnings = FALSE)
    saveRDS(res,paste0(mainDir,"/",cell,"Wilcox_results.rds"))
  }
}


plot_patientwise_results<-function(deconvolution_df,deconvolution_method,results_dir){
  patients<-names(table(deconvolution_df$patient.id)[table(deconvolution_df$patient.id)>8])

  
  for (patient in patients){
    patient_data<-deconvolution_df[deconvolution_df$patient.id==patient,]
    patient_data$Proportion<-as.numeric(patient_data$Proportion)
    patient_data<-aggregate(patient_data$Proportion, by=list(patient_data$passage,patient_data$cell.type),
              FUN=mean, na.rm=TRUE)
    colnames(patient_data)<-c("passage","cell.type","Proportion")
    
    p<-ggplot(patient_data,aes(fill=cell.type,x=passage,y=Proportion))+ 
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

get_PDXe_results<-function(dat_exprs,pdxs,metadata,method,results_dir){
  dir.create(results_dir, showWarnings = FALSE)
  if(method=="consensus_tme"){
    deconvolution<-deconvolute(dat_exprs[,pdxs],method,indications=rep("brca",length(pdxs)))
  }else{
    deconvolution<-deconvolute(dat_exprs[,pdxs],method)
  }
  if(method=="mcp_counter"){
    #deconvolution[-1]<-deconvolution[-1]/colSums(deconvolution[-1])
  }
  deconvolution_2<-rbind(t(metadata)[,pdxs],deconvolution[,-1])
  rownames(deconvolution_2)<-c(unlist(colnames(metadata)),deconvolution$cell_type)
  deconvolution.df<-as.data.frame(t(deconvolution_2[-c(1,4),]))
  #multiple.samples<-names(table(breast_deconvolution.df$patient.id)[table(breast_deconvolution.df$patient.id)>1])
  #breast_deconvolution.df<-breast_deconvolution.df[breast_deconvolution.df$patient.id %in% multiple.samples,]
  deconvolution.df<-melt(deconvolution.df,id.vars=c("patient.id","passage"),variable.name="cell.type",value.name="Proportion")
  
  print("Plotting results")
  plot_results(deconvolution.df,method,results_dir)
  #wilcox_test_results(breast_deconvolution.df,method)
  plot_patientwise_results(deconvolution.df,method,results_dir)
  saveRDS(deconvolution.df,paste0(results_dir,"/Deconvolution_results.rds"))
}




#deconvolution<-deconvolute(exprs_data.not_log,"mcp_counter") #Proportion of total immune cells (not including tumour/endothelial tissue), sums to 1
#deconvolution<-colSums(deconvolute(exprs_data.not_log,"mcp_counter")[-1]) #Raw count
#colSums(deconvolute(exprs_data.not_log,"xcell")[-1]) #Proportion not including uncharacterized (sums to less than 1)
#colSums(deconvolute(exprs_data.not_log,"quantiseq")[-1]) #Proportion including uncharacterized
#colSums(deconvolute(exprs_data.not_log,"epic")[-1]) #Proportion including uncharacterized
#deconvolution<-deconvolute(exprs_data.not_log,"cibersort_abs") #Proportion not including uncharacterized (sums to less than 1)

get_PDXe_results(breast_exprs_data.not_log,breast_pdxs,breast_metadata,"cibersort","../Results/PDXe/Breast")
get_PDXe_results(breast_exprs_data.not_log,breast_pdxs,breast_metadata,"cibersort_abs","../Results/PDXe/Breast")
get_PDXe_results(breast_exprs_data.not_log,breast_pdxs,breast_metadata,"xcell","../Results/PDXe/Breast")
get_PDXe_results(breast_exprs_data.not_log,breast_pdxs,breast_metadata,"quantiseq","../Results/PDXe/Breast")
get_PDXe_results(breast_exprs_data.not_log,breast_pdxs,breast_metadata,"mcp_counter","../Results/PDXe/Breast") #Convert to percentages
get_PDXe_results(breast_exprs_data.not_log,breast_pdxs,breast_metadata,"epic","../Results/PDXe/Breast")
get_PDXe_results(breast_exprs_data.not_log,breast_pdxs,breast_metadata,"consensus_tme","../Results/PDXe/Breast")

get_PDXe_results(exprs_data.not_log,full_pdxs,metadata,"cibersort","../Results/PDXe/PanCancer")
get_PDXe_results(exprs_data.not_log,full_pdxs,metadata,"cibersort_abs","../Results/PDXe/PanCancer")
get_PDXe_results(exprs_data.not_log,full_pdxs,metadata,"xcell","../Results/PDXe/PanCancer")
get_PDXe_results(exprs_data.not_log,full_pdxs,metadata,"quantiseq","../Results/PDXe/PanCancer")
get_PDXe_results(exprs_data.not_log,full_pdxs,metadata,"mcp_counter","../Results/PDXe/PanCancer") #Convert to percentages
get_PDXe_results(exprs_data.not_log,full_pdxs,metadata,"epic","../Results/PDXe/PanCancer")

get_PDXe_results(exprs_data.not_log,full_pdxs,metadata,"cibersort","../Results/PDXe/PanCancer_Stratified")
get_PDXe_results(exprs_data.not_log,full_pdxs,metadata,"cibersort_abs","../Results/PDXe/PanCancer_Stratified")
get_PDXe_results(exprs_data.not_log,full_pdxs,metadata,"xcell","../Results/PDXe/PanCancer_Stratified")
get_PDXe_results(exprs_data.not_log,full_pdxs,metadata,"quantiseq","../Results/PDXe/PanCancer_Stratified")
get_PDXe_results(exprs_data.not_log,full_pdxs,metadata,"mcp_counter","../Results/PDXe/PanCancer_Stratified") #Convert to percentages
get_PDXe_results(exprs_data.not_log,full_pdxs,metadata,"epic","../Results/PDXe/PanCancer_Stratified")

gene_vars<-rowVars(naturally_logged)
all_genes<-rownames(naturally_logged)
variance_data<-data.frame(cbind(all_genes,gene_vars,c(1:length(all_genes))))
variance_data<-variance_data %>% group_by(all_genes)%>% slice_max(n = 1, gene_vars, with_ties = FALSE)
indices<-variance_data$V3
naturally_logged<-naturally_logged[as.numeric(indices),]


exprData_notlog<-exp(naturally_logged)
deconvolution<-deconvolute(exprData_notlog,"mcp_counter")

get_results(exprData_notlog,ccle_cells,breast_metadata,"cibersort","../Results/PDXe/Breast")
get_results(exprData_notlog,ccle_cells,breast_metadata,"cibersort_abs","../Results/PDXe/Breast")
get_results(exprData_notlog,ccle_cells,breast_metadata,"xcell","../Results/PDXe/Breast")
get_results(exprData_notlog,ccle_cells,breast_metadata,"quantiseq","../Results/PDXe/Breast")
get_results(exprData_notlog,ccle_cells,breast_metadata,"mcp_counter","../Results/PDXe/Breast") #Convert to percentages
get_results(exprData_notlog,ccle_cells,breast_metadata,"epic","../Results/PDXe/Breast")

#biobase.id, patient.id, passage, tumor.type
#Each passage is basically a new PDX, 0th is direct sample
#Subset and remove patients that have only 1 passage/sample
#Subset only to breast cancer
#Look at cibersort results for different passages over time
#natural log