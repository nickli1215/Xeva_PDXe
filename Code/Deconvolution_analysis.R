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

set_cibersort_binary("G:/UHN/CIBERSORT/CIBERSORT.R")
set_cibersort_mat("G:/UHN/CIBERSORT/LM22.txt")
pset_folder<-"G:/UHN/psets"

PDX_exprs<-readRDS('../Data/PDXE_microArray.rds')
metadata<-pData(PDX_exprs)

breast_pdxs<-metadata[metadata$tumor.type=="breast",]$biobase.id
exprs_data<-exprs(PDX_exprs)[,breast_pdxs]
exprs_data.not_log<-exp(exprs_data)
breast_metadata<-metadata[breast_pdxs,]

PDX_exprs<-readRDS('../Data/pdxMorag_rnaseq.rda')
metadata<-pData(PDX_exprs)

plot_results<-function(deconvolution_df,deconvolution_method,results_dir,xlab="Passage no."){
  all_cell_plots<-list()
  
  for (cell in unique(deconvolution_df$cell.type)){
    cell_data<-deconvolution_df[deconvolution_df$cell.type==cell,]
    cell_data$passage<-as.factor(cell_data$passage)
    cell_data$Proportion<-as.numeric(as.vector(cell_data$Proportion))

    #ggplot(data=cell_data, aes(x=passage, y=Proportion, group=patient.id)) +
      #geom_line() + geom_point()+
    #ggplot(data=cell_data, aes(x=passage, y=Proportion,)) + 
      #geom_violin() +
      #scale_color_brewer(palette="Paired")+
     # theme_minimal() + ggtitle(paste0(deconvolution_method,": ",cell))
    p<-ggboxplot(cell_data, x = "passage", y = "Proportion", 
              ylab = "Cell Proportion", xlab =xlab,title=cell)+ font("title",size=8) 
    
    all_cell_plots[[cell]]<-p
  }
  ml <- marrangeGrob(all_cell_plots, nrow=2, ncol=2, top=deconvolution_method)
  
  
  mainDir<-paste0(results_dir,"/",deconvolution_method)
  dir.create(mainDir, showWarnings = FALSE)
  ggsave(paste0(mainDir,"/",deconvolution_method,"_Full_boxplot.pdf"),plot=ml, width = 4, height = 4)
    
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

get_PDXe_results<-function(dat_exprs,pdxs,metadata,method,results_dir){
  dir.create(results_dir, showWarnings = FALSE)
  deconvolution<-deconvolute(dat_exprs,method)
  if(method=="mcp_counter"){
    #deconvolution[-1]<-deconvolution[-1]/colSums(deconvolution[-1])
  }
  breast_deconvolution<-rbind(t(metadata)[,pdxs],deconvolution[,-1][,pdxs])
  rownames(breast_deconvolution)<-c(unlist(colnames(metadata)),deconvolution$cell_type)
  breast_deconvolution.df<-as.data.frame(t(breast_deconvolution[-c(1,4),]))
  #multiple.samples<-names(table(breast_deconvolution.df$patient.id)[table(breast_deconvolution.df$patient.id)>1])
  #breast_deconvolution.df<-breast_deconvolution.df[breast_deconvolution.df$patient.id %in% multiple.samples,]
  breast_deconvolution.df<-melt(breast_deconvolution.df,id.vars=c("patient.id","passage"),variable.name="cell.type",value.name="Proportion")
  
  print("Plotting results")
  plot_results(breast_deconvolution.df,method,results_dir)
  #wilcox_test_results(breast_deconvolution.df,method)
  
  #return(breast_deconvolution.df)
}

plot_cellline_results<-function(deconvolution_df,deconvolution_method){
  for (cell in unique(deconvolution_df$cell.type)){
    cell_data<-deconvolution_df[deconvolution_df$cell.type==cell,]
    #ggplot(data=cell_data, aes(x=passage, y=Proportion, group=patient.id)) +
    #geom_line() + geom_point()+
    #ggplot(data=cell_data, aes(x=passage, y=Proportion,)) + 
    #geom_violin() +
    #scale_color_brewer(palette="Paired")+
    # theme_minimal() + ggtitle(paste0(deconvolution_method,": ",cell))
    ggboxplot(cell_data, x = "passage", y = "Proportion", 
              ylab = "Cell Proportion", xlab = "Passage no.")
    
    mainDir<-paste0("../Results/",deconvolution_method)
    dir.create(mainDir, showWarnings = FALSE)
    print(cell)
    if(grepl( "/",cell, fixed = TRUE)){
      cell<-str_replace_all(cell,"/","+")
    }
    ggsave(paste0(mainDir,"/",cell,"_boxplot.pdf"), width = 4, height = 4)
    
  }
}


#deconvolution<-deconvolute(exprs_data.not_log,"mcp_counter") #Proportion of total immune cells (not including tumour/endothelial tissue), sums to 1
#deconvolution<-colSums(deconvolute(exprs_data.not_log,"mcp_counter")[-1]) #Raw count
#colSums(deconvolute(exprs_data.not_log,"xcell")[-1]) #Proportion not including uncharacterized (sums to less than 1)
#colSums(deconvolute(exprs_data.not_log,"quantiseq")[-1]) #Proportion including uncharacterized
#colSums(deconvolute(exprs_data.not_log,"epic")[-1]) #Proportion including uncharacterized
#deconvolution<-deconvolute(exprs_data.not_log,"cibersort_abs") #Proportion not including uncharacterized (sums to less than 1)

get_PDXe_results(exprs_data.not_log,breast_pdxs,breast_metadata,"cibersort","../Results/PDXe")
get_PDXe_results(exprs_data.not_log,breast_pdxs,breast_metadata,"cibersort_abs","../Results/PDXe")
get_PDXe_results(exprs_data.not_log,breast_pdxs,breast_metadata,"xcell","../Results/PDXe")
get_PDXe_results(exprs_data.not_log,breast_pdxs,breast_metadata,"quantiseq","../Results/PDXe")
get_PDXe_results(exprs_data.not_log,breast_pdxs,breast_metadata,"mcp_counter","../Results/PDXe") #Convert to percentages
get_PDXe_results(exprs_data.not_log,breast_pdxs,breast_metadata,"epic","../Results/PDXe")


ccle<-readRDS(paste0(pset_folder,"/CCLE.rds"))
breast_cell_info <- ccle@cell[ccle@cell$unique.tissueid.fromstudies=="breast",]
ccle_cells <- breast_cell_info[breast_cell_info$tissueid=="Breast","cellid"]
ccle_cells<-ccle_cells[!is.na(ccle_cells)]
ccle_subset<-subsetTo(ccle,cells=ccle_cells,molecular.data.cells=ccle_cells)
exprData<-molecularProfiles(ccle_subset,mDataType="Kallisto_0.46.1.rnaseq")
naturally_logged<-exprData*log(2)
exprData_notlog<-exp(naturally_logged)
deconvolution<-deconvolute(exprData_notlog,"mcp_counter")

get_results(exprData_notlog,ccle_cells,breast_metadata,"cibersort","../Results/PDXe")
get_results(exprData_notlog,ccle_cells,breast_metadata,"cibersort_abs","../Results/PDXe")
get_results(exprData_notlog,ccle_cells,breast_metadata,"xcell","../Results/PDXe")
get_results(exprData_notlog,ccle_cells,breast_metadata,"quantiseq","../Results/PDXe")
get_results(exprData_notlog,ccle_cells,breast_metadata,"mcp_counter","../Results/PDXe") #Convert to percentages
get_results(exprData_notlog,ccle_cells,breast_metadata,"epic","../Results/PDXe")

#biobase.id, patient.id, passage, tumor.type
#Each passage is basically a new PDX, 0th is direct sample
#Subset and remove patients that have only 1 passage/sample
#Subset only to breast cancer
#Look at cibersort results for different passages over time
#natural log