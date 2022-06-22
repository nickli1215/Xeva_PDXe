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

set_cibersort_binary("G:/UHN/CIBERSORT/CIBERSORT.R")
set_cibersort_mat("G:/UHN/CIBERSORT/LM22.txt")
pset_folder<-"G:/UHN/psets"
CCLE<-readRDS(paste0(pset_folder,"/CCLE.rds"))
CCLE.breast.cells<-rownames(CCLE@cell[CCLE@cell$tissueid=="Breast",])
CCLE<-subsetTo(CCLE,cells=CCLE.breast.cells)
CCLE.experiment<-summarizeMolecularProfiles(CCLE,mDataType="Kallisto_0.46.1.rnaseq")
gene_metadata<-CCLE.experiment@elementMetadata

CCLE.exprs<-assay(CCLE.experiment)
CCLE.exprs<-CCLE.exprs[gene_metadata$gene_id,]
rownames(CCLE.exprs)<-gene_metadata$gene_name
CCLE.exprs<-CCLE.exprs[gene_metadata$gene_type=="protein_coding",]
CCLE.exprs<-CCLE.exprs[,colSums(is.na(CCLE.exprs))==0]
CCLE.exprs<-2**CCLE.exprs - 0.001 


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




get_pvalue<-function(x,data){
  group1<-x[1]
  group2<-x[2]
  group_1_data<-data[data$Dataset==group1,]
  group_2_data<-data[data$Dataset==group2,]
  rownames(group_1_data)<-group_1_data$cellline
  rownames(group_2_data)<-group_2_data$cellline
  
  cells<-intersect(group_1_data$cellline,group_2_data$cellline)
  
  group_1_data<-group_1_data[cells,]
  group_2_data<-group_2_data[cells,]
  
  wilc.test<-wilcox.test(group_1_data$Proportion,group_2_data$Proportion,paired=TRUE)
  return(wilc.test$p.value)
}



deconvolution_methods<-c("cibersort_abs","xcell","quantiseq","mcp_counter","epic","consensus_tme")
for(deconvolution_method in deconvolution_methods){
  if(deconvolution_method=="consensus_tme"){
    
    CCLE.deconvolution<-deconvolute(CCLE.exprs,deconvolution_method,indications=rep("brca",ncol(CCLE.exprs)))
    rownames(CCLE.deconvolution)<-CCLE.deconvolution$cell_type
    GDSC.deconvolution<-deconvolute(GDSC.exprs,deconvolution_method,indications=rep("brca",ncol(GDSC.exprs)))
    rownames(GDSC.deconvolution)<-GDSC.deconvolution$cell_type
    gCSI.deconvolution<-deconvolute(gCSI.exprs,deconvolution_method,indications=rep("brca",ncol(gCSI.exprs)))
    rownames(gCSI.deconvolution)<-gCSI.deconvolution$cell_type
    GRAY.deconvolution<-deconvolute(GRAY.exprs,deconvolution_method,indications=rep("brca",ncol(GRAY.exprs)))
    rownames(GRAY.deconvolution)<-GRAY.deconvolution$cell_type
  }else{
    CCLE.deconvolution<-deconvolute(CCLE.exprs,deconvolution_method)
    rownames(CCLE.deconvolution)<-CCLE.deconvolution$cell_type
    GDSC.deconvolution<-deconvolute(GDSC.exprs,deconvolution_method)
    rownames(GDSC.deconvolution)<-GDSC.deconvolution$cell_type
    gCSI.deconvolution<-deconvolute(gCSI.exprs,deconvolution_method)
    rownames(gCSI.deconvolution)<-gCSI.deconvolution$cell_type
    GRAY.deconvolution<-deconvolute(GRAY.exprs,deconvolution_method)
    rownames(GRAY.deconvolution)<-GRAY.deconvolution$cell_type
  }
  
  all.cells<-unique(c(gCSI.breast.cells,gdsc.breast.cells,CCLE.breast.cells,GRAY.breast.cells))
  all.cells<-all.cells[all.cells %in% colnames(gCSI.deconvolution) + 
                         all.cells %in% colnames(GDSC.deconvolution) + 
                         all.cells %in% colnames(CCLE.deconvolution) + 
                         all.cells %in% colnames(GRAY.deconvolution) >1]
  
  all_df<-matrix(ncol = 6, nrow = 0)
  x <- c("cellline","cell_type","CCLE","GDSC","gCSI","GRAY")
  colnames(all_df) <- x
  
  for (cellline in all.cells){
    if(cellline %in% colnames(CCLE.deconvolution)){
      CCLE_data<-CCLE.deconvolution[,cellline]
    }else{
      CCLE_data<-rep(NA,length(rownames(CCLE.deconvolution)))
    }
    if(cellline %in% colnames(GDSC.deconvolution)){
      GDSC_data<-GDSC.deconvolution[,cellline]
    }else{
      GDSC_data<-rep(NA,length(rownames(CCLE.deconvolution)))
    }
    if(cellline %in% colnames(gCSI.deconvolution)){
      gCSI_data<-gCSI.deconvolution[,cellline]
    }else{
      gCSI_data<-rep(NA,length(rownames(CCLE.deconvolution)))
    }
    if(cellline %in% colnames(GRAY.deconvolution)){
      GRAY_data<-GRAY.deconvolution[,cellline]
    }else{
      GRAY_data<-rep(NA,length(rownames(CCLE.deconvolution)))
    }
    cell_types<-CCLE.deconvolution$cell_type
    all_data<-cbind(rep(cellline,length(rownames(CCLE.deconvolution))),cell_types,CCLE_data,GDSC_data,gCSI_data,GRAY_data)
    colnames(all_data)<-x
    all_df<-rbind(all_df,all_data)
  }
  
  
  all_cell_plots<-list()
  
  for (cell in unique(all_df$cell_type)){
    cell_data<-all_df[all_df$cell_type==cell,]
    cell_data_long<-melt(cell_data,id=c("cellline","cell_type"),value.name = "Proportion")
    colnames(cell_data_long)[3]<-"Dataset"
    cell_data_long<-cell_data_long[!is.na(cell_data_long$Proportion),]
    ylim1 = boxplot.stats(cell_data_long$Proportion)$stats[c(1, 5)]
    
    my_comparisons<-list (c("CCLE","GDSC"), c("GDSC", "gCSI"),  c("gCSI", "GRAY"), 
                          c("CCLE","gCSI"),c("GDSC","GRAY"),
                          c("CCLE","GRAY"))
    groups<-do.call(rbind, my_comparisons)
    pvalues<-lapply(my_comparisons,function(x){return(get_pvalue(x,cell_data_long))})
    pvalues<-round(unlist(pvalues),5)
    .y.<-rep("Proportion",length(pvalues))
    y.position<-c(rep(ylim1[2]+0.05,3),ylim1[2]+0.1,ylim1[2]+0.15,ylim1[2]+0.2)
    y.position<-as.numeric(y.position)
    
    stats<-data.frame(.y.=.y.,group1=as.vector(groups[,1]),
                      group2=as.vector(groups[,2]),
                      p=unlist(pvalues),
                      y.position=as.vector(y.position))
    stats<-as_tibble(stats)
    p0<-ggboxplot(cell_data_long, x = "Dataset", y = "Proportion",
                 color="Dataset",id="cellline",
                 palette = "jco")+
      font("title",size=12)+ggtitle(cell)+
      #stat_compare_means(aes(group=Dataset),method="wilcox.test",comparisons=my_comparisons,paired=TRUE)+
      stat_pvalue_manual(data=stats,label="p")+
     theme_bw()
    ylim1 = boxplot.stats(cell_data_long$Proportion)$stats[c(1, 5)]
    #p0 = p0 + coord_cartesian(ylim = ylim1*1.05)
    
    all_cell_plots[[cell]]<-p0
  }
  ml<-ggarrange(plotlist=all_cell_plots,nrow=3,ncol=2,common.legend=TRUE)
  results_dir<-"../Results/CellLines"
  mainDir<-paste0(results_dir,"/",deconvolution_method)
  dir.create(mainDir, showWarnings = FALSE)
  pdf(paste0(mainDir,"/",deconvolution_method,"_Full_boxplot.pdf"),height=9,width=6)
  print(ml)
  dev.off() 
  
  
  
  all_celllines_plots<-list()
  for (cell in unique(all_df$cell_type)){
    cell_data<-all_df[all_df$cell_type==cell,]
    plotting_data<-as.matrix(cell_data[,c("CCLE","GDSC","gCSI","GRAY")])
    rownames(plotting_data)<-cell_data$cellline
    p<-Heatmap(plotting_data,cluster_rows=FALSE,
               column_order=c("CCLE","gCSI","GDSC","GRAY"),
               na_col = "black",
               show_column_dend = FALSE,
               column_title = cell)
    if(grepl( "/",cell, fixed = TRUE)){
      cell<-str_replace_all(cell,"/","+")
    }
    dir.create(paste0(mainDir,"/heatmaps"),showWarnings = FALSE)
    pdf(paste0(mainDir,"/","heatmaps/",cell,"_heatmap.pdf"),height=10,width=5)
    print(p)
    dev.off()
  }
  
  for (cellline in all.cells){
    cellline_data<-all_df[all_df$cellline==cellline,]
    cellline_data<-melt(cellline_data)
    colnames(cellline_data)<-c("Cell-line","cell.type","Dataset","Proportion")
    cellline_data<-cellline_data[!is.na(cellline_data$Proportion),]
    cellline_data$Proportion<-as.numeric(cellline_data$Proportion)
    p<-ggplot(cellline_data,aes(fill=cell.type,x=Dataset,y=Proportion))+ 
      geom_bar(position="fill", stat="identity")+
      geom_col(width=1)
    
    mainDir<-paste0(results_dir,"/",deconvolution_method)
    dir.create(mainDir, showWarnings = FALSE)
    dir.create(paste0(mainDir,"/barplots"),showWarnings = FALSE)
    pdf(paste0(mainDir,"/","barplots/",cellline,"_barplot.pdf"),height=10,width=5)
    print(p)
    dev.off()
  }
  
}
