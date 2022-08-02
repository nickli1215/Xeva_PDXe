library(genefu)
library(Biobase)
library(biomaRt)
library(org.Hs.eg.db)

data(pam50.robust)

PDXe.data<-readRDS("../Data/PDXE_microArray.rds")
PDXe.pdata<-pData(PDXe.data)
PDXe.data<-PDXe.data[,PDXe.pdata$tumor.type=="breast"]
PDXe.exprs<-exprs(PDXe.data)
feature.data<-fData(PDXe.data)
PDXe.pdata<-pData(PDXe.data)


PDXe.exprs <- t(exprs(PDXe.data)) 
#PDXe.exprs<-(PDXe.exprs-rowMeans(PDXe.exprs))/sqrt(rowVars(PDXe.exprs))
entrez.ids<-mapIds(org.Hs.eg.db, colnames(PDXe.exprs), 'ENTREZID', 'SYMBOL')
feature.data<-feature.data[colnames(PDXe.exprs),]
feature.data$EntrezGene.ID<-unname(entrez.ids)

PDXe.sbt <- molecular.subtyping(sbt.model="pam50", data=PDXe.exprs,
                                annot=feature.data, do.mapping=F)
PDXe.pdata[["subtype"]]<-PDXe.sbt$subtype


#PDXe.pdata<- PDXe.pdata %>% group_by (patient.id)

PDXe.subtype.info<-reshape(PDXe.pdata[,c("patient.id","passage","subtype")],
        idvar="patient.id",
        timevar="passage",
        direction="wide")
PDXe.subtype.info<-PDXe.subtype.info[,c("patient.id","subtype.0","subtype.1","subtype.4","subtype.5")]
rownames(PDXe.subtype.info)<-PDXe.subtype.info$patient.id
PDXe.subtype.4.5<-PDXe.subtype.info %>% mutate(mycol=coalesce(subtype.4,subtype.5)) %>%
  select(patient.id, mycol)
four.five.vec <- PDXe.subtype.4.5$mycol
PDXe.subtype.info<-PDXe.subtype.info[,c("patient.id","subtype.0","subtype.1")]
colnames(PDXe.subtype.info)<-c("Patient ID","0","1-3")

PDXe.subtype.info$`4-5`<-four.five.vec
PDXe.subtype.info<-PDXe.subtype.info[rowSums(!is.na(PDXe.subtype.info))>2,c(2,3,4)]
colors = structure(c(7,5,2,3), names = c("Basal","LumB","LumA","Her2")) # black, red, green, blue
pdf("../Results/PDXe_passage_subtyping.PDF",width=5,height=10)
Heatmap(PDXe.subtype.info,
        col=colors,
        name="subtype",
        row_title = "Patient ID",
        row_names_side = "left",
        column_title = "PDX Passage no.", 
        column_title_side = "top",
        column_names_rot = 0,
        column_names_side="top")
dev.off()


