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
library(ggthemes)
library(dplyr)

##logistic model
t1<-read.csv('../TCGA/TCGA-CESC/TCGA-CESC-taxonomy.normalized-01A-Genus.csv',check.names = F,header = T,row.names = 1)
colnames(t1)<-substr(colnames(t1),1,12)
row.names(t1)<-substr(row.names(t1),4,1111)

t2<-read_xlsx('TCGA-CESC clinicaldata.xlsx') %>% as.data.frame()
table(t2$lymphovascular.invasion.indicator)
row.names(t2)<-t2[,2]
t2$`LSVI&N`[t2$N.stage=='N0'&t2$lymphovascular.invasion.indicator=='ABSENT']<-'LSVI(-)&N(-)'
t2$`LSVI&N`[t2$N.stage=='N0'&t2$lymphovascular.invasion.indicator=='PRESENT']<-'LSVI(+)&N(-)'
t2$`LSVI&N`[t2$N.stage=='N1'&t2$lymphovascular.invasion.indicator=='PRESENT']<-'LSVI(+)&N(+)'
t2<-t2[,c(27:28,33:35)]
t2<-t2[t2$OS.time!=0,]
colnames(t2)

s1<-c("Malassezia","Actinomyces","Gordonia","Nocardioides","Sinorhizobium","Xanthomonas")#,"Sinorhizobium"
t3<-cbind(t2[intersect(row.names(t2),colnames(t1)),],t(t1[s1,intersect(row.names(t2),colnames(t1))]))
t3<-t3[is.na(t3$`LSVI&N`)==F,]
t3$Group<-t3$`LSVI&N`
t3$Group[t3$Group=="LSVI(+)&N(-)"]<-1
t3$Group[t3$Group=="LSVI(+)&N(+)"]<-1
t3$Group[t3$Group=="LSVI(-)&N(-)"]<-0
str(t3)
t3$Group<-as.numeric(t3$Group)

p3<-glm(Group~ Malassezia+Actinomyces+Gordonia+Nocardioides+Xanthomonas, family= binomial(), data=t3)
p4=summary(p3)
p5=cbind(coef=p4$coefficients[,'Estimate'],
         OR=exp(p4$coefficients[,'Estimate']),
         OR.95L=exp(p4$coefficients[,'Estimate']-p4$coefficients[,'Std. Error']),
         OR.95H=exp(p4$coefficients[,'Estimate']+p4$coefficients[,'Std. Error']),
         pvalue=p4$coefficients[,"Pr(>|z|)"])
p5<-as.data.frame(p5)
p5$coef<-as.numeric(p5$coef)
write.csv(p5,'./logistic-5genes.csv')

t3$riskscore<-t3[,"Malassezia"]*p5["Malassezia",1][[1]]+t3[,"Actinomyces"]*p5["Actinomyces",1][[1]]+t3[,"Gordonia"]*p5["Gordonia",1][[1]]+
  t3[,"Nocardioides"]*p5["Nocardioides",1][[1]]+t3[,"Xanthomonas"]*p5["Xanthomonas",1][[1]]#+t3[,"Sinorhizobium"]*p5["Sinorhizobium",1][[1]]

t3$OS.time<-t3$OS.time/30
t3$PFI.time<-t3$PFI.time/30
m1<-maxstat.test(Surv(OS.time,OS)~ riskscore, data = t3,smethod="LogRank",pmethod="HL")
plot(m1)
m2<-cbind(m1[["stats"]],m1[["cuts"]])
t3$group<-ifelse(t3$riskscore>0.18565989,'2high','1low')
fit1 <- survfit(Surv(OS.time,OS)~group, data =t3)
write.csv(t3,'./tcga-riskscore.csv.csv')

t55<-t3[order(t3$riskscore),]
t55$riskscore
p11 <- ggplot(data = t55) +
  geom_point(aes(x=seq(1:nrow(t55)), y=riskscore, color=group)) + 
  scale_color_manual(values = c("#A6CEE3","#1F78B4")) + 
  ggtitle("Risk Score") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  xlab("") + 
  ylab("Risk Score") +
  geom_vline(aes(xintercept=87),colour = "#BB0000",linetype = "dashed")+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        # axis.text.x = element_text(angle=30, hjust=1, vjust=1),
        axis.line = element_line(color="black", size=0.4, linetype="solid"),
        legend.position = 'none')
pdf(paste('./TCGA-riskscurve','.pdf',sep = ''),width = 3.3,height = 3.3,onefile = F,paper = 'a4')
print(p11)
dev.off()

outTab=data.frame()
cox <- coxph(Surv(OS.time,OS)~group, data =t3)
coxSummary = summary(cox)
coxP=coxSummary$coefficients[,"Pr(>|z|)"]
outTab=rbind(outTab,
             cbind(id='riskscore',coef=coxSummary$coefficients[,"coef"],
                   z=coxSummary$coefficients[,"z"],
                   HR=coxSummary$conf.int[,"exp(coef)"],
                   HR.95L=coxSummary$conf.int[,"lower .95"],
                   HR.95H=coxSummary$conf.int[,"upper .95"],
                   pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
)
write.csv(outTab,'./riskscore-OS-HR.csv',row.names = F)
pp1<-ggsurvplot(fit1,
                fun = "pct", size = 1,
                xlab='Months',
                ylab='OS Probability (%)',
                break.time.by=24,
                risk.table = T,risk.table.y.text.col = T,risk.table.y.text = T,
                #surv.median.line = "hv",
                pval = T,
                risk.table.fontsize=4,
                font.legend=c(12),
                font.xlab = c(12),font.y = c(12),
                font.tickslab = 12,
                palette = c("#A6CEE3","#1F78B4"),
                legend.labs=c('Low','High'),
                legend.title="Exp")
fit2 <- survfit(Surv(PFI.time,PFI)~group, data =t3)
pp2<-ggsurvplot(fit2,
                fun = "pct", size = 1,
                xlab='Months',
                ylab='PFS Probability (%)',
                break.time.by=24,
                risk.table = T,risk.table.y.text.col = T,risk.table.y.text = T,
                #surv.median.line = "hv",
                pval = T,
                risk.table.fontsize=4,
                font.legend=c(12),
                font.xlab = c(12),font.y = c(12),
                font.tickslab = 12,
                palette = c("#A6CEE3","#1F78B4"),
                legend.labs=c('Low','High'),
                legend.title="Exp")
pp3<-ggarrange(pp2$plot, pp2$table,heights = c(3, 1),
               labels = c("",""),align = "v",
               ncol = 1, nrow = 2)
pp4<-ggarrange(pp1$plot, pp1$table,heights = c(3, 1),
               labels = c("",""),align = "v",
               ncol = 1, nrow = 2)
pp5<-ggarrange(pp3, pp4,labels = c("PFS","OS"),
               ncol = 1, nrow = 2)
pdf(paste('./km-riskscore from LVSI&N -OS&PFI','.pdf',sep = ''),width = 4.8,height = 12,onefile = F,paper = 'a4')
print(pp5)
dev.off()

##Cox model
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
library(ggthemes)
t1<-read.csv('../TCGA/TCGA-CESC/TCGA-CESC-taxonomy.normalized-01A-Genus.csv',check.names = F,header = T,row.names = 1)
colnames(t1)<-substr(colnames(t1),1,12)
row.names(t1)<-substr(row.names(t1),4,1111)

t2<-read_xlsx('./TCGA-CESC clinicaldata.xlsx') %>% as.data.frame()
table(t2$lymphovascular.invasion.indicator)
row.names(t2)<-t2[,2]
t2$`LSVI&N`[t2$N.stage=='N0'&t2$lymphovascular.invasion.indicator=='ABSENT']<-'LSVI(-)&N(-)'
t2$`LSVI&N`[t2$N.stage=='N0'&t2$lymphovascular.invasion.indicator=='PRESENT']<-'LSVI(+)&N(-)'
t2$`LSVI&N`[t2$N.stage=='N1'&t2$lymphovascular.invasion.indicator=='PRESENT']<-'LSVI(+)&N(+)'
t2<-t2[,c(27:28,33:35)]
t2<-t2[t2$OS.time!=0,]
colnames(t2)

s1<-c("Malassezia","Actinomyces","Gordonia","Nocardioides","Sinorhizobium","Xanthomonas")#,"Sinorhizobium"
t3<-cbind(t2[intersect(row.names(t2),colnames(t1)),],t(t1[s1,intersect(row.names(t2),colnames(t1))]))

p3=coxph(Surv(OS.time, OS) ~ Malassezia+Actinomyces+Gordonia+Nocardioides+Xanthomonas, data = t3)
p4=summary(p3)
p5=cbind(coef=p4$coefficients[,"coef"],
         HR=p4$conf.int[,"exp(coef)"],
         HR.95L=p4$conf.int[,"lower .95"],
         HR.95H=p4$conf.int[,"upper .95"],
         pvalue=p4$coefficients[,"Pr(>|z|)"])
p5<-as.data.frame(p5)
p5$coef<-as.numeric(p5$coef)
i=row.names(t3)[1]
t3$riskscore<-t3[,"Malassezia"]*p5["Malassezia",1][[1]]+t3[,"Actinomyces"]*p5["Actinomyces",1][[1]]+t3[,"Gordonia"]*p5["Gordonia",1][[1]]+
  t3[,"Nocardioides"]*p5["Nocardioides",1][[1]]+t3[,"Xanthomonas"]*p5["Xanthomonas",1][[1]]#+t3[,"Sinorhizobium"]*p5["Sinorhizobium",1][[1]]

t3$OS.time<-t3$OS.time/30
t3$PFI.time<-t3$PFI.time/30
m1<-maxstat.test(Surv(OS.time,OS)~ riskscore, data = t3,smethod="LogRank",pmethod="HL")
plot(m1)
m2<-cbind(m1[["stats"]],m1[["cuts"]])
t3$group<-ifelse(t3$riskscore>0.170,'2high','1low')
fit1 <- survfit(Surv(OS.time,OS)~group, data =t3)
pp1<-ggsurvplot(fit1,
                fun = "pct", size = 1,
                xlab='Months',
                ylab='OS Probability (%)',
                break.time.by=24,
                risk.table = T,risk.table.y.text.col = T,risk.table.y.text = T,
                #surv.median.line = "hv",
                pval = T,
                risk.table.fontsize=4,
                font.legend=c(12),
                font.xlab = c(12),font.y = c(12),
                font.tickslab = 12,
                palette = c("#A6CEE3","#1F78B4"),
                legend.labs=c('Low','High'),
                legend.title="Exp")
fit2 <- survfit(Surv(PFI.time,PFI)~group, data =t3)
pp2<-ggsurvplot(fit2,
                fun = "pct", size = 1,
                xlab='Months',
                ylab='PFS Probability (%)',
                break.time.by=24,
                risk.table = T,risk.table.y.text.col = T,risk.table.y.text = T,
                #surv.median.line = "hv",
                pval = T,
                risk.table.fontsize=4,
                font.legend=c(12),
                font.xlab = c(12),font.y = c(12),
                font.tickslab = 12,
                palette = c("#A6CEE3","#1F78B4"),
                legend.labs=c('Low','High'),
                legend.title="Exp")
pp3<-ggarrange(pp2$plot, pp2$table,heights = c(3, 1),
               labels = c("",""),align = "v",
               ncol = 1, nrow = 2)
pp4<-ggarrange(pp1$plot, pp1$table,heights = c(3, 1),
               labels = c("",""),align = "v",
               ncol = 1, nrow = 2)
pp5<-ggarrange(pp3, pp4,labels = c("PFS","OS"),
               ncol = 1, nrow = 2)
pdf(paste('./Cox-km-riskscore-OS&PFI','.pdf',sep = ''),width = 4.8,height = 12,onefile = F,paper = 'a4')
print(pp5)
dev.off()










