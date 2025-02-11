rm(list=ls())
library(readxl)
library(ggplot2)
library(openxlsx)
library(limma)
library(ggrepel)

meta<-read_excel("data/meta/all_meta.xlsx")
dia<-read.xlsx("./data/clean_data/exp_dia.xlsx",rowNames = TRUE)
phos<-read.xlsx("./data/clean_data/exp_phos.xlsx",rowNames = TRUE)

normal.meta<-meta[meta$tissue=="normal",]
tumor.meta<-meta[meta$tissue=="tumor",]

dia<-log2(dia+1)
phos<-log2(phos+1)


normal.dia<-dia[,normal.meta$ID]
normal.phos<-phos[,normal.meta$ID]
tumor.dia<-dia[,tumor.meta$ID]
tumor.phos<-phos[,tumor.meta$ID]
tumor.pCR<-tumor.meta$pCR
normal.pCR<-normal.meta$pCR



exp<-tumor.phos
y<-tumor.pCR
y[y==0]<-"RD"
y[y==1]<-"pCR"


group_list<-factor(y,ordered = F)

# 创建设计矩阵，指定组别信息
design<-model.matrix(~0+factor(group_list))
colnames(design) <- levels(factor(group_list))
rownames(design)<-colnames(exp)


#寻找差异性蛋白质
fit<-lmFit(exp,design)
#谁比较谁
con<-paste(unique(levels(group_list)), collapse = "-")
cont.matrix <- makeContrasts(contrasts = con,levels = design)
fit2 <- contrasts.fit(fit,cont.matrix) 
fit2 <- eBayes(fit2) 
tmpOut <- topTable(fit2,coef = 1,n = Inf,adjust.method = "BH",sort.by = "logFC") #coef指定哪个
limma.na <- na.omit(tmpOut)


# 设置pvalue
cut_off_pvalue=0.05
# 根据阈值分别为上调基因设置‘up’，下调基因设置‘Down’，无差异设置‘None’，保存到change列
tmpOut$sig <- "not"
tmpOut$sig[which((tmpOut$P.Value < 0.05) & (tmpOut$logFC > 1.5))] = "up"
tmpOut$sig[which((tmpOut$P.Value< 0.05) & (tmpOut$logFC < -1))] = "down"
#取出具有差异的
diff<-subset(tmpOut,sig %in% c("down","up"))
table(diff$sig)
write.csv(tmpOut,file="~/result/火山图/phos_all.csv",row.names = TRUE)
write.csv(diff,file="~/result/火山图/phos_DEPs.csv",row.names = TRUE)




################好看的火山图###########
p1<-ggplot(tmpOut, aes(x = logFC, y = -log10(P.Value)))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")+
  geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "black")+
  geom_point(aes(size = -log10(P.Value), 
                 color = -log10(P.Value)))+
  scale_color_gradientn(values = seq(0,1,0.2),
                        colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+
  scale_size_continuous(range = c(0,1))+
  theme_bw(base_size = 12)+#改变字体大小
  theme(panel.grid = element_blank(),
        legend.position = 'right',
        legend.justification = c(0,1))+
  # 设置图例
  guides(col = 
           guide_colorbar(title = "-Log10(P.value)",
                          ticks.colour = NA,
                          reverse = T,
                          title.vjust = 0.8,
                          barheight = 8,
                          barwidth = 1),
         size = "none") +
  # 添加标签：
  xlab("Log2FC")+
  ylab("-Log10(P.Value)")

##选择显著的加上标签
data<-tmpOut
data$row<-rownames(data)
data$label<-rep(NA,nrow(data))
data$label[order(abs(data$logFC), decreasing = T)[1:4]] <- data$row[order(abs(data$logFC), decreasing = T)[1:4]]
label<-c('AFTPH','NFAM1','ANKRD54','TFAP2A','NELFB','GALNT7')

geneList <- data[rownames(data) %in% label,]
geneList$label<-label

p2<-p1+ geom_text_repel(data = geneList,
               aes(x = logFC, y = -log10(P.Value), label = label),
                  #这里的filter很关键，筛选你想要标记的基因
              size = 3,color="black",
              box.padding = unit(0.4, "lines"), 
              segment.color = "black",   #连线的颜色
              segment.size = 0.2,  #连线的粗细
              max.overlaps = 50
)
ggsave("~/result/火山图/phos_火山图.pdf",plot = p2,width = 4, height = 3)





















