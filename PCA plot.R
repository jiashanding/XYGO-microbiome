a1<-read.csv('./Xiangya RNA-seq-tpm.csv',header = T)
a1<-a1[duplicated(a1$symbol)==F,]
a1<-a1[a1$type=='protein_coding',]
row.names(a1)<-a1$symbol
a2<-read_xlsx('./Xiangya cohort-RNAseq&microbiome clinicaldata.xlsx')%>%as.data.frame()
row.names(a2)<-a2$id
library(ggplot2)
library(factoextra)
library(FactoMineR)

a2<-cbind(a2,t(a1[,row.names(a2)]))
a2<-a2[,c(-1,-3)]
a2$Group[a2$Group=="LVSI(+)&N(+)"]<-"LVSI(+)&LM(+)"
a2$Group[a2$Group=="LVSI(+)&N(-)"]<-"LVSI(+)&LM(-)"
a2$Group[a2$Group=="LVSI(-)&N(-)"]<-"LVSI(-)&LM(-)"
table(a2$Group)
# a2$Group<-factor(a2$Group,levels = c("LVSI(-)&LM(-)","LVSI(+)&LM(-)","LVSI(+)&LM(+)"))
a2.pca<- PCA(a2[,-1], graph = FALSE)
pdf('./PCA.plot.pdf',height = 5,width = 6.6,paper = 'a4')
fviz_pca_ind(a2.pca,
             geom.ind = "point",
             pointsize =2.5,
             pointshape = 21,
             fill.ind = a2[,1],
             palette = c("#F5A974","#73C876","#B09CD0"),
             addEllipses = TRUE, 
             legend.title = "Type",
             title="")+
  # scale_y_continuous(limits=c(-15,15))+
  # scale_x_continuous(limits=c(-25,30))+
  theme_bw() +
  theme(text=element_text(size=12,face="plain",color="black"),
        axis.title=element_text(size=13,face="plain",color="black"),
        axis.text = element_text(size=12,face="plain",color="black"),
        legend.title = element_text(size=13,face="plain",color="black"),
        legend.text = element_text(size=12,face="plain",color="black"),
        legend.background = element_blank(),
        panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())
dev.off()



