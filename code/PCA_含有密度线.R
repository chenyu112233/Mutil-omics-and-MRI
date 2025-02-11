rm(list=ls())
library(FactoMineR)
library(ggplot2)
library(ggrepel)
library(readxl)
library(tidyverse)
library(aplot)
library(vegan)
library(patchwork)
library(cowplot)
library(grid)
library(ggplotify)
#需要读取差异分析的exp以及meta
list<-diff[diff$sig !="not",]
exp<-data_exp[data_exp$gene %in% list$gene,]
write.xlsx(exp,"Protein_DEP_exp.xlsx")


group<-meta[,1:2]
PCAmat<-exp[,2:ncol(exp)]
rownames(PCAmat)<-exp$gene
PCAmat<-t(PCAmat)#需要专置一下，行为样本，列为基因
gene.pca<-PCA(PCAmat,ncp=2,scale.unit = TRUE,graph=FALSE)

#提取样本在 PCA 前两轴中的坐标
pca_sample <- data.frame(gene.pca$ind$coord[ ,1:2])
pca_sample$patient_ID=row.names(pca_sample)
#提取 PCA 前两轴的贡献度
pca_eig1 <- round(gene.pca$eig[1,2], 2)
pca_eig2 <- round(gene.pca$eig[2,2],2 )
pca_sample <- merge(pca_sample, group,by="patient_ID")
head(pca_sample)


color=c("#5686C3","#75C500")
p <- ggplot(data = pca_sample, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(color = tissue), size = 2) +  #根据样本坐标绘制二维散点图
  scale_color_manual(values = color) +  #自定义颜色
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +  #去除背景和网格线
  labs(x =  paste('PCA1(', pca_eig1, '%)'), y = paste('PCA2(', pca_eig2, '%)'), color = '')  #将 PCA 轴贡献度添加到坐标轴标题中

###添加密度线
Fig1a.function.pc1.density<-
  ggplot(data = pca_sample)+
  geom_density(aes(x=Dim.1,group=tissue,fill=tissue),
               color="black",alpha=0.6,position='identity',
               show.legend=F)+
  scale_fill_manual(values=c("#5686C3","#75C500"))+
  scale_linetype_manual(values=c("solid","dashed"))+
  scale_y_discrete(expand=c(0,0.001))+
  labs(x=NULL,y=NULL)+
  theme_classic()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

Fig1a.function.pc2.density<-
  ggplot(data = pca_sample)+
  geom_density(aes(x=Dim.2,group=tissue,fill=tissue),
               color="black",alpha=0.6,position='identity',show.legend=F)+
  scale_fill_manual(values=c("#5686C3","#75C500"))+
  theme_classic()+
  scale_linetype_manual(values=c("solid","dashed"))+
  scale_y_discrete(expand=c(0,0.001))+
  labs(x=NULL,y=NULL)+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  coord_flip()

p2 <-p %>% 
  insert_top(Fig1a.function.pc1.density,height =0.3)%>%
  insert_right(Fig1a.function.pc2.density,width=0.3)%>%
  as.ggplot()
#添加名字
p <- p + 
  geom_text_repel(aes(label = patient_ID), size = 1, show.legend = FALSE, 
                  box.padding = unit(0.2, 'lines'))
#画95的置信区间
p<-p + stat_ellipse(aes(color = tissue), level = 0.95, show.legend = FALSE)

#保存为pdf
pdf("Protein_PCA.pdf",height=3,width = 4)
p
dev.off()

pdf("RNA_PCA_2.pdf",height=3,width = 4)
p2
dev.off()
