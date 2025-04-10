a1.1<-read.csv('../RNA-seq/整理/microbiome/xiangya RNAseq-taxonomy.normalized-LVSI(-)&N(-)-Genus.csv',sep = ',',header = T,row.names = 1)
a2.1<-read.csv('../RNA-seq/整理/microbiome/xiangya RNAseq-taxonomy.normalized-LVSI(+)&N(-)-Genus.csv',sep = ',',header = T,row.names = 1)
a3.1<-read.csv('../RNA-seq/整理/microbiome/xiangya RNAseq-taxonomy.normalized-LVSI(+)&N(+)-Genus.csv',sep = ',',header = T,row.names = 1)

a1<-rbind(a1.1,a1.2,a1.3)
a2<-rbind(a2.1,a2.2,a2.3)
a3<-rbind(a3.1,a3.2,a3.3)

rm1 <- apply(a1, 1, function(x){sum(x != 0 ) > ncol(a1)*0.5})
a1.rm<-a1[rm1==T,]
rm2 <- apply(a2, 1, function(x){sum(x != 0 ) > ncol(a2)*0.5})
a2.rm<-a2[rm2==T,]
rm3 <- apply(a3, 1, function(x){sum(x != 0 ) > ncol(a3)*0.5})
a3.rm<-a3[rm3==T,]

x<-list(`LVSI(-)&N(-)`=row.names(a1.rm),
        `LVSI(+)&N(-)`=row.names(a2.rm),
        `LVSI(+)&N(+)`=row.names(a3.rm))
venn <- Venn(x)
data <- process_data(venn)
saveRDS(data,'./Vnn-Xiangya RNAseq.rds')

pdf('./Vnn-Xiangya RNAseq.pdf',height = 5.5,width = 5.5,paper = 'a4')
ggVennDiagram(x, label = c("count"),edge_size = 0.6,
              size=0.5,lty="longdash",
              color="red") + 
  scale_fill_gradient(name="Count",low="#D1E5F0",high = "#2166AC") +
  hrbrthemes::theme_ipsum(base_family = "sans") +
  labs(title = "",
       subtitle = "",
       caption = "") +
  theme(plot.title = element_text(hjust = 0.5,vjust = .5,color = "black",face = 'bold',
                                  size = 20, margin = margin(t = 1, b = 12)),
        plot.subtitle = element_text(hjust = 0,vjust = .5,size=15),
        plot.caption = element_text(face = 'bold',size = 12),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),legend.position = 'none'
  )
dev.off()

v6<-c(data[["regionData"]][["item"]][[1]],data[["regionData"]][["item"]][[2]],data[["regionData"]][["item"]][[3]])#unique

i=row.names(a1)[1]
f0<-matrix(NA,ncol = 3,nrow = nrow(a1))
row.names(f0)<-row.names(a1)
colnames(f0)<-c('LVSI(-)&N(-)','LVSI(+)&N(-)','LVSI(+)&N(+)')
for (i in row.names(a1)) {
  f0[i,'LVSI(-)&N(-)']<-length(which(a1[i,]!=0))/ncol(a1)
  f0[i,'LVSI(+)&N(-)']<-length(which(a2[i,]!=0))/ncol(a2)
  f0[i,'LVSI(+)&N(+)']<-length(which(a3[i,]!=0))/ncol(a3)
}
f0<-as.data.frame(f0)
f0.v6<-f0[v6,]

pdf(file="./pheatmap-Xiangya RNAseq-frequency.pdf",width = 3,height = 5,paper = 'a4')
pheatmap(as.matrix(f0.v6),
         color = colorRampPalette(c('white','#1A488E'))(10),
         # annotation_col = y4.2,
         # annotation_colors = col1,
         # treeheight_row = 50,
         # gaps_col = 120,
         show_rownames = F,
         # name = "Exp",
         scale = 'row',
         show_colnames = T,
         fontsize_row = 7,
         angle_col = 45,
         cluster_rows = F,
         cluster_cols = F)
dev.off()

f1<-cbind(a1,a2,a3)
for (i in colnames(f1)) {
  f1[,i][f1[,i]>0.5]<-0.5
}
f1.v6<-f1[v6,]
f1.v5<-f1[v5,]
f1.v4<-f1[v4,]
ann.v6<-matrix(c(rep('LVSI(-)&LM(-)',ncol(a1)),rep('LVSI(+)&LM(-)',ncol(a2)),rep('LVSI(+)&LM(+)',ncol(a3))),ncol = 1,nrow = ncol(f1)
)
row.names(ann.v6)<-c(colnames(a1),colnames(a2),colnames(a3))
ann.v6<-as.data.frame(ann.v6)
names(ann.v6)[1]<-'Type'
col.v6<-list(Type=c(`LVSI(-)&LM(-)`="#F5A974",`LVSI(+)&LM(-)`="#73C876",`LVSI(+)&LM(+)`="#B09CD0"))

pdf('./fig1-TCGA分析/颜色adj-pheatmap-score-宫颈癌新测RNAseq-unique.pdf',width = 5,height = 5,paper = 'a4')
pheatmap(as.matrix(f1.v6),
         color = colorRampPalette(c('white','#1A488E'))(10),
         annotation_col = ann.v6,
         annotation_colors = col.v6,
         # treeheight_row = 50,
         # gaps_col = 120,
         show_rownames = F,
         # name = "Exp",
         # scale = 'column',
         show_colnames = F,
         fontsize_row = 7,
         angle_col = 0,
         cluster_rows = F,
         cluster_cols = F)
dev.off()
write.csv(f1.v6,'./fig1-TCGA分析/颜色adj-pheatmap-score-宫颈癌新测RNAseq-unique.csv')

pdf('./fig1-TCGA分析/pheatmap-score-宫颈癌新测RNAseq-共有.pdf',width = 6.5,height = 5,paper = 'a4')
pheatmap(as.matrix(f1.v5),
         color = colorRampPalette(c('white','#1A488E'))(10),
         annotation_col = ann.v6,
         annotation_colors = col.v6,
         # treeheight_row = 50,
         # gaps_col = 120,
         show_rownames = F,
         # name = "Exp",
         # scale = 'column',
         show_colnames = F,
         fontsize_row = 7,
         angle_col = 0,
         cluster_rows = F,
         cluster_cols = F)
dev.off()


f2<-matrix(NA,ncol = 3,nrow = nrow(a1))
f2[,1]<-rowMeans(a1)
f2[,2]<-rowMeans(a2)
f2[,3]<-rowMeans(a3)
row.names(f2)<-row.names(a1)
colnames(f2)<-c('LVSI(-)&LM(-)','LVSI(+)&LM(-)','LVSI(+)&LM(+)')
for (i in colnames(f2)) {
  f2[,i][f2[,i]>0.5]<-0.5
}
f2.v6<-f2[v6,]

pdf(file="./pheatmap-Xiangya RNAseq-score.pdf",width = 3.7,height = 5,paper = 'a4')
pheatmap(as.matrix(f2.v6),
         color = colorRampPalette(c('white',"#2166AC"))(10),
         # annotation_col = y4.2,
         # annotation_colors = col1,
         # treeheight_row = 50,
         # gaps_col = 120,
         show_rownames = F,
         # name = "Exp",
         scale = 'row',
         show_colnames = T,
         fontsize_row = 7,
         angle_col = 0,
         cluster_rows = F,
         cluster_cols = F)
dev.off()

a6<-read.csv('./Abundance.filtered.anno.xls',sep = '\t',header = T,check.names = F) 
# a6<-a6[,1:7]

a6.1<-a6[,c("Phylum",a2[a2$Group=='LVSI(-)&N(-)',][,1])]
rm6.1 <- apply(a6.1[,2:ncol(a6.1)], 1, function(x){sum(x != 0 ) > 0})
a6.1<-a6.1[rm6.1==T,]
a7.1<-table(a6.1$Phylum) %>% as.data.frame()
a7.1$Type<-'LVSI(-)&LM(-)'
a7.1$Freq<-(a7.1$Freq/sum(a7.1$Freq))*100

a6.2<-a6[,c("Phylum",a2[a2$Group=='LVSI(+)&N(-)',][,1])]
rm6.2 <- apply(a6.2[,2:ncol(a6.2)], 1, function(x){sum(x != 0 ) > 0})
a6.2<-a6.2[rm6.2==T,]
a7.2<-table(a6.2$Phylum) %>% as.data.frame()
a7.2$Type<-'LVSI(+)&LM(-)'
a7.2$Freq<-(a7.2$Freq/sum(a7.2$Freq))*100

a6.3<-a6[,c("Phylum",a2[a2$Group=='LVSI(+)&N(+)',][,1])]
rm6.3 <- apply(a6.3[,2:ncol(a6.3)], 1, function(x){sum(x != 0 ) > 0})
a6.3<-a6.3[rm6.3==T,]
a7.3<-table(a6.3$Phylum) %>% as.data.frame()
a7.3$Type<-'LVSI(+)&LM(+)'
a7.3$Freq<-(a7.3$Freq/sum(a7.3$Freq))*100

a7<-rbind(a7.1,a7.2,a7.3)
colnames(a7)
table(a7$Type)
a7$Type<-factor(a7$Type,levels = c('LVSI(-)&LM(-)','LVSI(+)&LM(-)','LVSI(+)&LM(+)'))

pdf('./2bRAD-M-bar-Phylum.pdf',height = 4.5,width = 4.8,paper = 'a4')
ggplot(a7,aes(x=Type,y=Freq, fill=Var1))+
  #geom_histogram( position="fill",stat='identity',color="white",width=0.8)+
  geom_bar(stat='identity',position='stack',width = 0.8)+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        # axis.text.x = element_blank(),
        # axis.ticks.x  = element_blank(),
        # axis.text.y = element_blank(),
        # axis.ticks.y  = element_blank(),
        legend.position = 'right',
        axis.text.x = element_text(angle=30, hjust=1, vjust=1),
        axis.line = element_blank())+
  scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(limit=c('LVSI(-)&LM(-)','LVSI(+)&LM(-)','LVSI(+)&LM(+)'),
                   expand = c(0,0.4))+
  scale_fill_manual(values = c("#FBB4AE","#B3CDE3","#7cc0a7",'#85B8CB','#F9E54E',"#B2182B","#90a0c7",'#91D4C2','#5EAA5F','#199FB1',"#9daf91",
                               '#D1DDDB',"#e5d4cd",'#b59283'))+
  #geom_text(aes(label=count2), nudge_y= +2000, color="black", size = 2)+
  ylab('')+
  xlab('')+
  guides(fill=guide_legend(title="Phylum"))
dev.off()

#####################
a3<-read_xlsx('./xiangya cohort-clinicaldata.xlsx') %>% as.data.frame()
a3$Type[row.names(a3) %in% a2[a2$Group=='LVSI(-)&N(-)',][,1]]<-'LVSI(-)&N(-)'
a3$Type[row.names(a3) %in% a2[a2$Group=='LVSI(+)&N(-)',][,1]]<-'LVSI(+)&N(-)'
a3$Type[row.names(a3) %in% a2[a2$Group=='LVSI(+)&N(+)',][,1]]<-'LVSI(+)&N(+)'
a3.melt<-reshape2::melt(a3,id='Type')
a3.melt$Type<-factor(a3.melt$Type,levels = c('LVSI(-)&N(-)','LVSI(+)&N(-)','LVSI(+)&N(+)'))

p1<-ggplot(a3.melt[a3.melt$variable=='shannon',],aes(x=Type,y=value,fill=Type,color=Type))+
  geom_boxplot(outlier.shape = NA,fatten = 1)+
  #geom_jitter(shape = 16)+
  scale_fill_manual(values = c("#F5A974","#73C976","#B19CD0"))+
  scale_color_manual(values = c("#F5A974","#73C976","#B19CD0"))+
  xlab('')+ylab('shannon')+
  stat_compare_means(label="p.format",method = 'kruskal.test')+
  #geom_signif(comparisons = list(c('LVSI(-)&N(-)','LVSI(+)&N(+)')),test = 'wilcox.test',map_signif_level=T,color="black")+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        # axis.text.x = element_blank(),
        # axis.ticks.x  = element_blank(),
        #axis.text.y = element_blank(),
        #axis.ticks.y  = element_blank(),
        axis.line = element_line(color="black", size=0.4, linetype="solid"),
        axis.text.x = element_text(angle=30, hjust=1, vjust=1),
        legend.position = 'none')
p2<-ggplot(a3.melt[a3.melt$variable=='chao1',],aes(x=Type,y=value,fill=Type,color=Type))+
  geom_boxplot(outlier.shape = NA,fatten = 1)+
  #geom_jitter(shape = 16)+
  scale_fill_manual(values = c("#F5A974","#73C976","#B19CD0"))+
  scale_color_manual(values = c("#F5A974","#73C976","#B19CD0"))+
  xlab('')+ylab('chao1')+
  stat_compare_means(label="p.format",method = 'kruskal.test')+
  # geom_signif(comparisons = list(c('LVSI(-)&N(-)','LVSI(+)&N(+)')),test = 'wilcox.test',map_signif_level=T,color="black")+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        # axis.text.x = element_blank(),
        # axis.ticks.x  = element_blank(),
        #axis.text.y = element_blank(),
        #axis.ticks.y  = element_blank(),
        axis.line = element_line(color="black", size=0.4, linetype="solid"),
        axis.text.x = element_text(angle=30, hjust=1, vjust=1),
        legend.position = 'none')
p3<-ggplot(a3.melt[a3.melt$variable=='simpson',],aes(x=Type,y=value,fill=Type,color=Type))+
  geom_boxplot(outlier.shape = NA,fatten = 1)+
  #geom_jitter(shape = 16)+
  scale_fill_manual(values = c("#F5A974","#73C976","#B19CD0"))+
  scale_color_manual(values = c("#F5A974","#73C976","#B19CD0"))+
  xlab('')+ylab('simpson')+
  stat_compare_means(label="p.format",method = 'kruskal.test')+
  # geom_signif(comparisons = list(c('LVSI(-)&N(-)','LVSI(+)&N(+)')),test = 'wilcox.test',map_signif_level=T,color="black")+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        # axis.text.x = element_blank(),
        # axis.ticks.x  = element_blank(),
        #axis.text.y = element_blank(),
        #axis.ticks.y  = element_blank(),
        axis.line = element_line(color="black", size=0.4, linetype="solid"),
        axis.text.x = element_text(angle=30, hjust=1, vjust=1),
        legend.position = 'none')
pdf('./2bRAD-M-shannon.pdf',height =3.5,width=2.1,paper = 'a4')
print(p1)
dev.off()
pdf('./2bRAD-M-chao1.pdf',height =3.5,width=2.1,paper = 'a4')
print(p2)
dev.off()
pdf('./2bRAD-M-simpson.pdf',height =3.5,width=2.1,paper = 'a4')
print(p3)
dev.off()






