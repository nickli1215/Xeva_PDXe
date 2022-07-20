source('LoadData.R')

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



#deconvolution_methods<-c("cibersort_abs","xcell","quantiseq","mcp_counter","epic","consensus_tme")
deconvolution_methods<-c("epic","xcell","quantiseq","mcp_counter")

for(deconvolution_method in deconvolution_methods){
  if(deconvolution_method=="consensus_tme"){
    
    CCLE.deconvolution<-deconvolute(CCLE.exprs,deconvolution_method,indications=rep("brca",ncol(CCLE.exprs)))
    CCLE.deconvolution<-as.data.frame(CCLE.deconvolution)
    rownames(CCLE.deconvolution)<-CCLE.deconvolution$cell_type
    CCLE.deconvolution<-melt(CCLE.deconvolution,
                             id.vars="cell_type",
                             value.name="Score",
                             variable.name="Sample.ID")
    CCLE.deconvolution$dataset<-"CCLE"
    
    GDSC.deconvolution<-deconvolute(GDSC.exprs,deconvolution_method,indications=rep("brca",ncol(GDSC.exprs)))
    GDSC.deconvolution<-as.data.frame(GDSC.deconvolution)
    rownames(GDSC.deconvolution)<-GDSC.deconvolution$cell_type
    GDSC.deconvolution<-melt(GDSC.deconvolution,
                             id.vars="cell_type",
                             value.name="Score",
                             variable.name="Sample.ID")
    GDSC.deconvolution$dataset<-"GDSC"
    
    gCSI.deconvolution<-deconvolute(gCSI.exprs,deconvolution_method,indications=rep("brca",ncol(gCSI.exprs)))
    gCSI.deconvolution<-as.data.frame(gCSI.deconvolution)
    rownames(gCSI.deconvolution)<-gCSI.deconvolution$cell_type
    gCSI.deconvolution<-melt(gCSI.deconvolution,
                             id.vars="cell_type",
                             value.name="Score",
                             variable.name="Sample.ID")
    gCSI.deconvolution$dataset<-"gCSI"
    
    GRAY.deconvolution<-deconvolute(GRAY.exprs,deconvolution_method,indications=rep("brca",ncol(GRAY.exprs)))
    GRAY.deconvolution<-as.data.frame(GRAY.deconvolution)
    rownames(GRAY.deconvolution)<-GRAY.deconvolution$cell_type
    GRAY.deconvolution<-melt(GRAY.deconvolution,
                             id.vars="cell_type",
                             value.name="Score",
                             variable.name="Sample.ID")
    GRAY.deconvolution$dataset<-"GRAY"
    
    PDXe.deconvolution<-deconvolute(PDXe.exprs,deconvolution_method,indications=rep("brca",ncol(PDXe.exprs)))
    PDXe.deconvolution<-as.data.frame(PDXe.deconvolution)
    rownames(PDXe.deconvolution)<-PDXe.deconvolution$cell_type
    PDXe.deconvolution<-melt(PDXe.deconvolution,
                             id.vars="cell_type",
                             value.name="Score",
                             variable.name="Sample.ID")
    PDXe.deconvolution$dataset<-"PDXe"

    #Morag.deconvolution<-deconvolute(Morag.exprs,deconvolution_method,indications=rep("brca",ncol(Morag.exprs)))
    Morag.PDX.deconvolution<-deconvolute(Morag.PDX.exprs,deconvolution_method,indications=rep("brca",ncol(Morag.PDX.exprs)))
    Morag.PDX.deconvolution<-as.data.frame(Morag.PDX.deconvolution)
    rownames(Morag.PDX.deconvolution)<-Morag.PDX.deconvolution$cell_type
    
    Morag.Tumor.deconvolution<-deconvolute(Morag.Tumor.exprs,deconvolution_method,indications=rep("brca",ncol(Morag.Tumor.exprs)))
    Morag.Tumor.deconvolution<-as.data.frame(Morag.Tumor.deconvolution)
    rownames(Morag.Tumor.deconvolution)<-Morag.Tumor.deconvolution$cell_type
    
    #Morag.PDX.deconvolution<-Morag.deconvolution[,c("cell_type",rownames(Morag.metadata[Morag.metadata$Source=="PDX",]))]
    #Morag.Tumor.deconvolution<-Morag.deconvolution[,c("cell_type",rownames(Morag.metadata[Morag.metadata$Source=="Patient",]))]
    
    #Morag.deconvolution<-Morag.deconvolution[,-1]
    Morag.PDX.deconvolution<-melt(Morag.PDX.deconvolution,
                             id.vars="cell_type",
                             value.name="Score",
                             variable.name="Sample.ID")
    Morag.PDX.deconvolution$dataset<-"Morag.PDX"
    
    Morag.Tumor.deconvolution<-melt(Morag.Tumor.deconvolution,
                                  id.vars="cell_type",
                                  value.name="Score",
                                  variable.name="Sample.ID")
    Morag.Tumor.deconvolution$dataset<-"Morag.Tumor"
    
    #TCGA.deconvolution<-deconvolute(TCGA.exprs,deconvolution_method,indications=rep("brca",ncol(TCGA.exprs)))
    #TCGA.deconvolution<-as.data.frame(TCGA.deconvolution)
    #rownames(TCGA.deconvolution)<-TCGA.deconvolution$cell_type
    #TCGA.deconvolution<-TCGA.deconvolution[,-1]
    #TCGA.deconvolution<-melt(TCGA.deconvolution,
                        #      id.vars="cell_type",
                        #      value.name="Score",
                        #      variable.name="Sample.ID")
    #TCGA.deconvolution$dataset<-"TCGA"
    
  }else{
    CCLE.deconvolution<-deconvolute(CCLE.exprs,deconvolution_method)
    CCLE.deconvolution<-as.data.frame(CCLE.deconvolution)
    rownames(CCLE.deconvolution)<-CCLE.deconvolution$cell_type
    CCLE.deconvolution<-melt(CCLE.deconvolution,
                             id.vars="cell_type",
                             value.name="Score",
                             variable.name="Sample.ID")
    CCLE.deconvolution$dataset<-"CCLE"
    
    GDSC.deconvolution<-deconvolute(GDSC.exprs,deconvolution_method)
    GDSC.deconvolution<-as.data.frame(GDSC.deconvolution)
    rownames(GDSC.deconvolution)<-GDSC.deconvolution$cell_type
    GDSC.deconvolution<-melt(GDSC.deconvolution,
                             id.vars="cell_type",
                             value.name="Score",
                             variable.name="Sample.ID")
    GDSC.deconvolution$dataset<-"GDSC"
    
    gCSI.deconvolution<-deconvolute(gCSI.exprs,deconvolution_method)
    gCSI.deconvolution<-as.data.frame(gCSI.deconvolution)
    rownames(gCSI.deconvolution)<-gCSI.deconvolution$cell_type
    gCSI.deconvolution<-melt(gCSI.deconvolution,
                             id.vars="cell_type",
                             value.name="Score",
                             variable.name="Sample.ID")
    gCSI.deconvolution$dataset<-"gCSI"
    
    GRAY.deconvolution<-deconvolute(GRAY.exprs,deconvolution_method)
    GRAY.deconvolution<-as.data.frame(GRAY.deconvolution)
    rownames(GRAY.deconvolution)<-GRAY.deconvolution$cell_type
    GRAY.deconvolution<-melt(GRAY.deconvolution,
                             id.vars="cell_type",
                             value.name="Score",
                             variable.name="Sample.ID")
    GRAY.deconvolution$dataset<-"GRAY"
    
    PDXe.deconvolution<-deconvolute(PDXe.exprs,deconvolution_method)
    PDXe.deconvolution<-as.data.frame(PDXe.deconvolution)
    rownames(PDXe.deconvolution)<-PDXe.deconvolution$cell_type
    PDXe.deconvolution<-melt(PDXe.deconvolution,
                             id.vars="cell_type",
                             value.name="Score",
                             variable.name="Sample.ID")
    PDXe.deconvolution$dataset<-"PDXe"
    
    Morag.PDX.deconvolution<-deconvolute(Morag.PDX.exprs,deconvolution_method)
    Morag.PDX.deconvolution<-as.data.frame(Morag.PDX.deconvolution)
    rownames(Morag.PDX.deconvolution)<-Morag.PDX.deconvolution$cell_type
    Morag.PDX.deconvolution<-melt(Morag.PDX.deconvolution,
                             id.vars="cell_type",
                             value.name="Score",
                             variable.name="Sample.ID")
    Morag.PDX.deconvolution$dataset<-"Morag.PDX"
    
    Morag.Tumor.deconvolution<-deconvolute(Morag.Tumor.exprs,deconvolution_method)
    Morag.Tumor.deconvolution<-as.data.frame(Morag.Tumor.deconvolution)
    rownames(Morag.Tumor.deconvolution)<-Morag.Tumor.deconvolution$cell_type
    Morag.Tumor.deconvolution<-melt(Morag.Tumor.deconvolution,
                             id.vars="cell_type",
                             value.name="Score",
                             variable.name="Sample.ID")
    Morag.Tumor.deconvolution$dataset<-"Morag.Tumor"
    
    #Morag.PDX.deconvolution<-Morag.deconvolution[,c("cell_type",rownames(Morag.metadata[Morag.metadata$Source=="PDX",]))]
    #Morag.Tumor.deconvolution<-Morag.deconvolution[,c("cell_type",rownames(Morag.metadata[Morag.metadata$Source=="Patient",]))]
    
    #Morag.deconvolution<-Morag.deconvolution[,-1]

    
    #TCGA.deconvolution<-deconvolute(TCGA.exprs,deconvolution_method)
    #TCGA.deconvolution<-as.data.frame(TCGA.deconvolution)
    #rownames(TCGA.deconvolution)<-TCGA.deconvolution$cell_type
    #TCGA.deconvolution<-melt(TCGA.deconvolution,
                                    #id.vars="cell_type",
                                    #value.name="Score",
                                    #variable.name="Sample.ID")
    #TCGA.deconvolution$dataset<-"TCGA"
    
  }
  df_list <- list(#TCGA.deconvolution,
                  PDXe.deconvolution,
                  Morag.PDX.deconvolution,
                  Morag.Tumor.deconvolution,
                  GRAY.deconvolution,
                  CCLE.deconvolution,
                  GDSC.deconvolution,
                  gCSI.deconvolution)
  all.deconvolution<-do.call("rbind", df_list)
  
  #all.cells<-unique(c(gCSI.breast.cells,gdsc.breast.cells,CCLE.breast.cells,GRAY.breast.cells))
 # all.cells<-all.cells[all.cells %in% colnames(gCSI.deconvolution) + 
                         #all.cells %in% colnames(GDSC.deconvolution) + 
                        # all.cells %in% colnames(CCLE.deconvolution) + 
                         #all.cells %in% colnames(GRAY.deconvolution) >1]
  
  #all_df<-matrix(ncol = 6, nrow = 0)
  #x <- c("cellline","cell_type","CCLE","GDSC","gCSI","GRAY")
  #colnames(all_df) <- x
  
  #for (cellline in all.cells){
    #if(cellline %in% colnames(CCLE.deconvolution)){
      #CCLE_data<-CCLE.deconvolution[,cellline]
    #}else{
     # CCLE_data<-rep(NA,length(rownames(CCLE.deconvolution)))
   # }
    #if(cellline %in% colnames(GDSC.deconvolution)){
      #GDSC_data<-GDSC.deconvolution[,cellline]
    #}else{
      #GDSC_data<-rep(NA,length(rownames(CCLE.deconvolution)))
   # }
    #if(cellline %in% colnames(gCSI.deconvolution)){
     # gCSI_data<-gCSI.deconvolution[,cellline]
    #}else{
      #gCSI_data<-rep(NA,length(rownames(CCLE.deconvolution)))
    #}
    #if(cellline %in% colnames(GRAY.deconvolution)){
      #GRAY_data<-GRAY.deconvolution[,cellline]
    #}else{
      #GRAY_data<-rep(NA,length(rownames(CCLE.deconvolution)))
    #}
    #cell_types<-CCLE.deconvolution$cell_type
    #all_data<-cbind(rep(cellline,length(rownames(CCLE.deconvolution))),cell_types,CCLE_data,GDSC_data,gCSI_data,GRAY_data)
    #colnames(all_data)<-x
    #all_df<-rbind(all_df,all_data)
  #}
  
  print("Plotting Results")

  all_cell_plots<-list()
  for (cell in unique(all.deconvolution$cell_type)){
    cell.data<-all.deconvolution[all.deconvolution$cell_type==cell,]
    p0<-ggboxplot(cell.data, x = "dataset", y = "Score",
      color="dataset",
      palette = "jco")+
      font("title",size=12)+ggtitle(cell)+
      stat_compare_means(label = "p.signif",
                         aes(group=dataset),
                         ref.group = "Morag.Tumor",
                         method="wilcox.test")+
      theme_bw()
    all_cell_plots[[cell]]<-p0
  }
  ml<-ggarrange(plotlist=all_cell_plots,nrow=3,ncol=2,common.legend=TRUE)
  results_dir<-"../Results/CellLines"
  mainDir<-paste0(results_dir,"/",deconvolution_method)
  dir.create(mainDir, showWarnings = FALSE)
  pdf(paste0(mainDir,"/",deconvolution_method,"_Full_boxplot_separated_noTCGA.pdf"),height=9,width=6)
  print(ml)
  dev.off() 
  
  #for (cell in unique(all_df$cell_type)){
    #cell_data<-all_df[all_df$cell_type==cell,]
    #cell_data_long<-melt(cell_data,id=c("cellline","cell_type"),value.name = "Proportion")
    #colnames(cell_data_long)[3]<-"Dataset"
    #cell_data_long<-cell_data_long[!is.na(cell_data_long$Proportion),]
    #ylim1 = boxplot.stats(cell_data_long$Proportion)$stats[c(1, 5)]
    
    #my_comparisons<-list (c("CCLE","GDSC"), c("GDSC", "gCSI"),  c("gCSI", "GRAY"), 
                          #c("CCLE","gCSI"),c("GDSC","GRAY"),
                          #c("CCLE","GRAY"))
    #my_comparisons
    #groups<-do.call(rbind, my_comparisons)
    #pvalues<-lapply(my_comparisons,function(x){return(get_pvalue(x,cell_data_long))})
    #pvalues<-round(unlist(pvalues),5)
    #.y.<-rep("Proportion",length(pvalues))
    #y.position<-c(rep(ylim1[2]+0.05,3),ylim1[2]+0.1,ylim1[2]+0.15,ylim1[2]+0.2)
    #y.position<-as.numeric(y.position)
    
   # stats<-data.frame(.y.=.y.,group1=as.vector(groups[,1]),
                      #group2=as.vector(groups[,2]),
                     # p=unlist(pvalues),
                      #y.position=as.vector(y.position))
    #stats<-as_tibble(stats)
    #p0<-ggboxplot(cell_data_long, x = "Dataset", y = "Proportion",
                 #color="Dataset",id="cellline",
                 #palette = "jco")+
      #font("title",size=12)+ggtitle(cell)+
      #stat_compare_means(aes(group=Dataset),method="wilcox.test",comparisons=my_comparisons,paired=TRUE)+
      #stat_pvalue_manual(data=stats,label="p")+
     #theme_bw()
    #ylim1 = boxplot.stats(cell_data_long$Proportion)$stats[c(1, 5)]
    #p0 = p0 + coord_cartesian(ylim = ylim1*1.05)
    
    #all_cell_plots[[cell]]<-p0
  #}

  
  
  #all_celllines_plots<-list()
  #for (cell in unique(all_df$cell_type)){
    #cell_data<-all_df[all_df$cell_type==cell,]
    #plotting_data<-as.matrix(cell_data[,c("CCLE","GDSC","gCSI","GRAY")])
    #rownames(plotting_data)<-cell_data$cellline
    #p<-Heatmap(plotting_data,cluster_rows=FALSE,
               #column_order=c("CCLE","gCSI","GDSC","GRAY"),
               #na_col = "black",
               #show_column_dend = FALSE,
               #column_title = cell)
    #if(grepl( "/",cell, fixed = TRUE)){
      #cell<-str_replace_all(cell,"/","+")
    #}
    #dir.create(paste0(mainDir,"/heatmaps"),showWarnings = FALSE)
    #pdf(paste0(mainDir,"/","heatmaps/",cell,"_heatmap.pdf"),height=10,width=5)
    #print(p)
    #dev.off()
  #}
  
  #for (cellline in all.cells){
    #cellline_data<-all_df[all_df$cellline==cellline,]
    #cellline_data<-melt(cellline_data)
    #colnames(cellline_data)<-c("Cell-line","cell.type","Dataset","Proportion")
    #cellline_data<-cellline_data[!is.na(cellline_data$Proportion),]
    #cellline_data$Proportion<-as.numeric(cellline_data$Proportion)
   # p<-ggplot(cellline_data,aes(fill=cell.type,x=Dataset,y=Proportion))+ 
      #geom_bar(position="fill", stat="identity")+
     # geom_col(width=1)
    
    #mainDir<-paste0(results_dir,"/",deconvolution_method)
    #dir.create(mainDir, showWarnings = FALSE)
    #dir.create(paste0(mainDir,"/barplots"),showWarnings = FALSE)
    #pdf(paste0(mainDir,"/","barplots/",cellline,"_barplot.pdf"),height=10,width=5)
    #print(p)
    #dev.off()
  #}
  
}
