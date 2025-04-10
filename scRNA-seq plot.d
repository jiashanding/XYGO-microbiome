setwd('E:/宫颈癌瘤内菌-最终/CC_scRNA_microbiome')
library(Seurat)
library(dplyr)
library(Matrix)
library(ggpubr)
library(magrittr)
library(cowplot)
library(ggsci)
library(ggplot2)
library(clustree)
library(harmony)
library(hdf5r)
library(reshape2)
library(readxl)
library(magrittr)
pbmc<-readRDS('./scRNA-seq metadata.rds')
i=c('KLF6','ARRDC3','DUSP6','PDZD8','VWA5A','ZNF215')
i=c('OVCH2','OXT','KCNT1','IFITM10','CHCHD10')
FeaturePlot(pbmc,reduction = "rna_harmony_umap",features=i,min.cutoff=0,max.cutoff=5,#pt.size = 0.8,
            cols = c("lightgrey","#053061"))+xlab('UMAP 1')+ylab('UMAP 2')

View(pbmc@meta.data)
colnames(pbmc@meta.data)
table(pbmc@meta.data$celltype)
pbmc.t<-subset(pbmc,subset = celltype=='T cell')
FeaturePlot(pbmc.t,reduction = "rna_harmony_umap",features=i,min.cutoff=0,max.cutoff=5,#pt.size = 0.8,
            cols = c("lightgrey","#053061"))+xlab('UMAP 1')+ylab('UMAP 2')

pdf('./annotation-distribution_harmony.pdf',width = 6.6,height = 5)
DimPlot(pbmc.t,reduction = "rna_harmony_umap",label=F,#split.by = "sample.type",
        group.by = "celltype2",
        cols = c("#BEBADA","#BC80BD","#4292C6","#8DD3C7","#FDB462","#B3DE69","#FB8072","#FCCDE5","#FBB4AE","#bcbc81"
                 ,'grey'))+
  # NoLegend()+
  # scale_y_continuous(limits=c(-4,11))+
  # scale_x_continuous(limits=c(-9,9))+
  xlab('UMAP 1')+
  ylab('UMAP 2')+
  ggtitle("")
dev.off()
table(pbmc@meta.data$tissue.type)
pbmc.t1<-subset(pbmc.t,subset = tissue.type=='Early satge')
a1<-AverageExpression(pbmc.t1,group.by="celltype2")
a1<-a1$RNA
a1<-as.data.frame(a1)
a1.2<-a1[c('KLF6','ARRDC3','DUSP6','PDZD8','VWA5A','ZNF215','OVCH2','OXT','KCNT1','IFITM10','CHCHD10'),]
a1.2<-as.matrix(a1.2)

pbmc.t2<-subset(pbmc.t,subset = tissue.type=='Advanced satge')
a2<-AverageExpression(pbmc.t2,group.by="celltype2")
a2<-a2$RNA
a2<-as.data.frame(a2)
a2.2<-a2[c('KLF6','ARRDC3','DUSP6','PDZD8','VWA5A','ZNF215','OVCH2','OXT','KCNT1','IFITM10','CHCHD10'),]
a2.2<-as.matrix(a2.2)

library(reshape2)
library(pheatmap)
a3<-cbind(a1.2[,c('CD4+ CCR7+ T','CD4+ Treg','CD8+ Tn','CD8+ GZMK+ T','CD8+ Tex')],
          a2.2[,c('CD4+ CCR7+ T','CD4+ Treg','CD8+ Tn','CD8+ GZMK+ T','CD8+ Tex')])
pdf('./heatmap-tcell-11genes.pdf',height = 4,width =4,paper = 'a4')
pheatmap(as.matrix(a3),
         color = colorRampPalette(c('#1A488E','white','#EE3239'))(20),
         # annotation_col = y4.2,
         # annotation_colors = col1,
         # treeheight_row = 50,
         gaps_col = 5,
         show_rownames = T,
         # name = "Exp",
         scale = 'row',
         show_colnames = T,
         fontsize_row = 8,
         angle_col = 45,
         cluster_rows = F,
         cluster_cols = F)
dev.off()






