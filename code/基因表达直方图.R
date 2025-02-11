rm(list=ls())
library(ggplot2)
library(tidyr)
library(cowplot)
library(gridExtra)

exp<-data
normalize_data <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

exp1<-exp[,4:25]
# 应用该函数到数据框的每一列
exp2<- as.data.frame(lapply(exp1, normalize_data))
exp3<-exp[,1:3]
data<-cbind(exp3,exp2)


p1<-ggplot(data,aes(x=group,y=CDKN2A,fill=tissue))+geom_boxplot()+theme(legend.position = "none")+theme(axis.title.x = element_blank())
pp2<-ggplot(data,aes(x=group,y=CILP2,fill=tissue))+geom_boxplot()+theme(legend.position = "none")+theme(axis.title.x = element_blank())
p3<-ggplot(data,aes(x=group,y=PLK1,fill=tissue))+geom_boxplot()+theme(legend.position = "none")+theme(axis.title.x = element_blank())
p4<-ggplot(data,aes(x=group,y=THBS2,fill=tissue))+geom_boxplot()+theme(legend.position = "none")+theme(axis.title.x = element_blank())
p5<-ggplot(data,aes(x=group,y=TMEM132A,fill=tissue))+geom_boxplot()+theme(legend.position = "none")+theme(axis.title.x = element_blank())
p6<-ggplot(data,aes(x=group,y=APOBR,fill=tissue))+geom_boxplot()+theme(legend.position = "none")+theme(axis.title.x = element_blank())
p7<-ggplot(data,aes(x=group,y=HMGA1,fill=tissue))+geom_boxplot()+theme(legend.position = "none")+theme(axis.title.x = element_blank())
p8<-ggplot(data,aes(x=group,y=TOMM6,fill=tissue))+geom_boxplot()+theme(legend.position = "none")+theme(axis.title.x = element_blank())
p9<-ggplot(data,aes(x=group,y=DERL3,fill=tissue))+geom_boxplot()+theme(legend.position = "none")+theme(axis.title.x = element_blank())
p10<-ggplot(data,aes(x=group,y=CBX2,fill=tissue))+geom_boxplot()+theme(legend.position = "none")+theme(axis.title.x = element_blank())
p11<-ggplot(data,aes(x=group,y=BCAT1,fill=tissue))+geom_boxplot()+theme(legend.position = "none")+theme(axis.title.x = element_blank())
p12<-ggplot(data,aes(x=group,y=MUC1,fill=tissue))+geom_boxplot()+theme(legend.position = "none")+theme(axis.title.x = element_blank())
p13<-ggplot(data,aes(x=group,y=RDH16,fill=tissue))+geom_boxplot()+theme(legend.position = "none")+theme(axis.title.x = element_blank())
p14<-ggplot(data,aes(x=group,y=ISG15,fill=tissue))+geom_boxplot()+theme(legend.position = "none")+theme(axis.title.x = element_blank())
p15<-ggplot(data,aes(x=group,y=LAD1,fill=tissue))+geom_boxplot()+theme(legend.position = "none")+theme(axis.title.x = element_blank())
p16<-ggplot(data,aes(x=group,y=FADS2,fill=tissue))+geom_boxplot()+theme(legend.position = "none")+theme(axis.title.x = element_blank())
p17<-ggplot(data,aes(x=group,y=SULT1E1,fill=tissue))+geom_boxplot()+theme(legend.position = "none")+theme(axis.title.x = element_blank())
p18<-ggplot(data,aes(x=group,y=S100A2,fill=tissue))+geom_boxplot()+theme(legend.position = "none")+theme(axis.title.x = element_blank())
p19<-ggplot(data,aes(x=group,y=SLC44A4,fill=tissue))+geom_boxplot()+theme(legend.position = "none")+theme(axis.title.x = element_blank())
p20<-ggplot(data,aes(x=group,y=APOBEC3A,fill=tissue))+geom_boxplot()+theme(legend.position = "none")+theme(axis.title.x = element_blank())
p21<-ggplot(data,aes(x=group,y=UBE2C,fill=tissue))+geom_boxplot()+theme(legend.position = "none")+theme(axis.title.x = element_blank())
p22<-ggplot(data,aes(x=group,y=NLRP2,fill=tissue))+geom_boxplot()+theme(legend.position = "none")+theme(axis.title.x = element_blank())




pdf("1.pdf",height = 20,width = 20)
plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,
          ncol = 5, # 设置列数
          nrow = 5, # 设置行数
          rel_widths = c(1, 1), # 设置每列的相对宽度
          rel_heights = rep(1, 1)) # 设置每行的相对高度
dev.off()

