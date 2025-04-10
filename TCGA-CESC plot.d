library(stringr)
library(tidyr)
library(readxl)
library(magrittr)
library(stringr)
library("maxstat")
library("survival")
library('survminer')
library(reshape2)
library(pheatmap)
library(data.table)
library(DESeq2)
library(ggVennDiagram)
library(ggplot2)
library(hrbrthemes)

#######################TCGA-survival
t1<-read.csv('../TCGA/TCGA-CESC/TCGA-CESC-taxonomy.normalized-01A-Genus.csv',check.names = F,header = T,row.names = 1)
colnames(t1)<-substr(colnames(t1),1,12)

t1$sum <- apply(t1, 1, sum)
t1<-t1[t1$sum!=0,]
t1<-t1[,-ncol(t1)]

t2<-read_xlsx('./TCGA-CESC.ClinicalData.xlsx') %>% as.data.frame()
table(t2$lymphovascular.invasion.indicator)
t2$`LSVI&N`[t2$N.stage=='N0'&t2$lymphovascular.invasion.indicator=='ABSENT']<-'LVSI(-)&LM(-)'
t2$`LSVI&N`[t2$N.stage=='N0'&t2$lymphovascular.invasion.indicator=='PRESENT']<-'LVSI(+)&LM(-)'
t2$`LSVI&N`[t2$N.stage=='N1'&t2$lymphovascular.invasion.indicator=='PRESENT']<-'LVSI(+)&LM(+)'
table(t2$`LSVI&N`)

t2$OS.time<-t2$OS.time/30
t2$PFI.time<-t2$PFI.time/30
t22<-t2[is.na(t2$`LSVI&N`)==F,]
t22$`LSVI&N`<-factor(t22$`LSVI&N`,levels = c('LVSI(-)&LM(-)','LVSI(+)&LM(-)','LVSI(+)&LM(+)'))
t22$group<-t22$`LSVI&N`
fit1 <- survfit(Surv(OS.time,OS)~group, data =t22)
p1<-ggsurvplot(fit1,
               fun = "pct", size = 1,
               xlab='Months',
               ylab='OS Probability (%)',
               break.time.by=24,
               risk.table = T,risk.table.y.text.col = T,risk.table.y.text = F,
               #surv.median.line = "hv",
               pval = T,
               risk.table.fontsize=4,
               font.legend=c(12),
               font.xlab = c(12),font.y = c(12),
               font.tickslab = 12,
               palette = c("#F5A974","#73C876","#B09CD0"),
               legend.labs=c('LVSI(-)&LM(-)','LVSI(+)&LM(-)','LVSI(+)&LM(+)'),
               legend.title="")
fit2 <- survfit(Surv(PFI.time,PFI)~group, data =t22)
p2<-ggsurvplot(fit2,
               fun = "pct", size = 1,
               xlab='Months',
               ylab='PFS Probability (%)',
               break.time.by=24,
               risk.table = T,risk.table.y.text.col = T,risk.table.y.text = F,
               #surv.median.line = "hv",
               pval = T,
               risk.table.fontsize=4,
               font.legend=c(12),
               font.xlab = c(12),font.y = c(12),
               font.tickslab = 12,
               palette = c("#F5A974","#73C876","#B09CD0"),
               legend.labs=c('LVSI(-)&LM(-)','LVSI(+)&LM(-)','LVSI(+)&LM(+)'),
               legend.title="")
p3<-ggarrange(p2$plot, p2$table,heights = c(3, 1),
              labels = c("",""),align = "v",
              ncol = 1, nrow = 2)
p4<-ggarrange(p1$plot, p1$table,heights = c(3, 1),
              labels = c("",""),align = "v",
              ncol = 1, nrow = 2)
p5<-ggarrange(p3, p4,labels = c("PFS","OS"),
              ncol = 1, nrow = 2)
pdf('./TCGA-LVSI&N-KM curve.pdf',width = 4.6,height = 12,onefile = F,paper = 'a4')
print(p5)
dev.off()

s1<-read.csv('./TCGA-survival.csv',row.names = 1,check.names = F,header = T)
s1<-s1[,c(6:10)]
s1<-s1[is.na(s1$cutoff_OS)==F,]
s1$Type<-'Non-sig'
s1$Type[s1$HR_cutoffOS_OS>1& s1$pvalue_cutoffOS_OS<0.05]<-'HR>1'
s1$Type[s1$HR_cutoffOS_OS<1& s1$pvalue_cutoffOS_OS<0.05]<-'HR<1'
for (i in s2) {
  s1$lable[row.names(s1)==i]<-substr(i,4,111)
}

pdf('./fig3/vol-TCGA-survival.pdf',height = 4.5,width = 4.2)
ggscatter(s1, x = 'HR_cutoffOS_OS', y = 'pvalue',
          color = "Type", 
          palette = c("#377EB8", "grey", "#E41A1C"),
          size = 1,
          label = s1$lable,
          font.label = 10, 
          repel = T,
          xlab = "HR value", 
          ylab = "-log10 (Pvalue)") + 
  # theme_base()+
  # ggtitle( this_tile )+
  # theme(plot.title = element_text(size=30,hjust = 0.5))+
  scale_x_continuous(limits = c(-2.5,5))+
  geom_hline(yintercept = c(-log10(0.05)), linetype="dashed") +
  geom_vline(xintercept = c(1), linetype="dashed")
dev.off()

library(ggraph)
library(tidytree)
t3<-read_xlsx('./phylogenetic tree.xlsx') %>% as.data.frame()

t4<-data.frame()
i=unique(t3$Genus)[1]
for (i in intersect(t3$Genus,substr(row.names(t1),4,111))) {
  t3.1<-t3[t3$Genus==i,]
  t4<-rbind(t4,t3.1)
}

t5<-data.frame()
for (i in 1:5) {
  t4.1<-t4[,i:(i+1)]
  names(t4.1)<-c('From','To')
  t5<-rbind(t5,t4.1)
}
mygraph <- graph_from_data_frame( t5 )
pdf('./phylogenetic tree.pdf',width = 8,height = 15,paper = 'a4')
ggraph(mygraph, layout = 'dendrogram', circular = F) + 
  geom_edge_elbow() +
  geom_node_point() +
  theme_graph()+coord_fixed()+coord_flip()+scale_y_reverse()+
  geom_node_text(aes(label=name),size=2.5)
dev.off()