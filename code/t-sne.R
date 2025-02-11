rm(list=ls())

clin1<-read.table("./tumor/clin_meta.txt",sep="\t",header=TRUE)
clin2<-read.table("./tumor_normal/clin_clean.csv",sep=',',header=TRUE,row.names = 1)

raw<-read.table("./tumor_normal/tumor_normal_DIA_clean.csv",sep=',',header=T,row.names = 1)

normalize_raw<-normalize(raw)


#t-SNE,列名是基因，行名是样本
library(Rtsne)
set.seed(321)
newraw<-t(raw)
tsne_out<-Rtsne(newraw,dims=2,pca=T,max_iter=1000,theta=0.4,perplexity=20,verbose=F)

#可视化
library(ggplot2)
tsne_result<-as.data.frame(tsne_out$Y)
colnames(tsne_result)<-c("tSNE1","tSNE2")
#画图
ggplot()+
  geom_point(tsne_result,
             aes(x=tSNE1,y=tSNE2,color=data$new))+
  stat_ellipse(tsne_result,geom="polygon"，
                             aes(x=tSNE1,y=tSNE2,
                                 group=data$new,
                                 fill=data$new),
                             alpha=0.5,
                             lty="dashed",
                             color="black",
                             key_glyph="blank")+theme_bw()


ggplot()+
  geom_point(data=tsne_result,
             aes(x=tSNE1,y=tSNE2,color=data$new))+
  stat_ellipse(data=tsne_result,
               geom="polygon",
               aes(x=tSNE1,y=tSNE2,
                   group=data$new,
                   fill=data$new),
               alpha=0.5,
               lty="dashed",
               color="black",
               key_glyph="blank")+
  theme_bw()

