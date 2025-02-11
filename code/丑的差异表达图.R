rm(list=ls())
library(readxl)
library(ggplot2)
library(openxlsx)
library(limma)
library(ggrepel)

tpms<-data[,-1]

label<-meta$tissue

exp<-tpms
rownames(exp)<-data$gene
y<-label


group_list<-factor(y,levels = c("tumor", "normal"))

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


# 包装函数：
# 调整1: xlim和ylim得去掉
# 调整2: 修改图例的位置
plotVoc <- function(diff){
  color <- rep("#999999",nrow(diff))
  
  color[diff$P.Value <0.05 & diff$logFC > 1] <- "#FC4E07"
  color[diff$P.Value <0.05 & diff$logFC < -1] <- "#00AFBB"
  
  par(oma = c(0,2,0,0))
  
  pdf("phos_差异分析.pdf",width = 4,height=4)
  plot(diff$logFC,-log10(diff$P.Value),pch = 16,cex = 0.5,
       col = color, frame.plot = F,
       xlab = "log2FC", ylab = "-log10(Pvalue)", cex.axis = 1, cex.lab = 1.3)
  
  # 添加参考线：
  abline(h = -log10(0.05),lwd = 2, lty = 3)  # lwd设置线的宽度，lty设置线的类型；
  abline(v = c(-1,1),lwd = 2, lty = 3)  # lwd设置线的宽度，lty设置线的类型；
  
  # 添加图例
  #legend(x = 3, y = max(-log10(diff$P.Value)), legend = c("not","up","down"), 
         #bty = "n", # 去除边框
         #pch = 19,cex = 1, # 设置点的样式和大小
         #x.intersp = 1, # 设置字与点之间的距离；
         #y.intersp = 1, # 设置点与点的高度差，相当于行距；
         #col = c("#999999", "#FC4E07","#00AFBB"))
  
  # 添加标签：
  diff_label<-c("EIF4A3","ARCN1","SLC11A1","PLEKHF2","VARS","CYB5R3","MTMR1","S100A13")
  geneList <- diff[rownames(diff) %in% diff_label,]
  geneList$label<-diff_label
  color = c()
  color[which(geneList$sig == "up")] = "black"
  color[which(geneList$sig == "down")] = "black"
  text(geneList$logFC,-log10(geneList$P.Value),
       labels = rownames(geneList),
       adj = c(0,1.5),
       cex = 0.4,
       col = color)
  dev.off()
}

diff1<-tmpOut
plotVoc(diff1)



#保存信息
write.xlsx(diff1,"protein_DEPs_all.xlsx",rowNames=TRUE)

