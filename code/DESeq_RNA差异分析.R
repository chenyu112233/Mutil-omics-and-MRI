rm(list=ls())
library("DESeq2")
 #data是表达矩阵，meta是,二者的行名和列名需要一样。
data1<-as.data.frame(data[,2:258])
rownames(data1)<-data$gene

meta1<-as.data.frame(meta[,2])
rownames(meta1)<-meta$patient_ID

colnames(meta1)[1]<-"condition"
meta1$condition<-factor(meta1$condition,levels = c("normal", "tumor"))


#差异分析,建立dds矩阵
dds<-DESeqDataSetFromMatrix(countData = data1,
                            colData = meta1,
                            design = ~condition)
dds<-DESeq(dds)#dds标准化
dds$condition

#提取差异分析的结果
res<-results(dds)
res<-res[order(res$padj),]

cut_off_pvalue=0.05
# 根据阈值分别为上调基因设置‘up’，下调基因设置‘Down’，无差异设置‘None’，保存到change列
res$sig <- "not"
res$sig[which((res$pvalue < 0.05) & (res$log2FoldChange > 1))] = "up"
res$sig[which((res$pvalue< 0.05) & (res$log2FoldChange < -1))] = "down"
res1<-na.omit(res)
table(res1$sig)

#看一下谁是上调下调
a<-plotCounts(dds,gene = "COL11A1",intgroup = "condition",returnData = TRUE)
pdf("test1.pdf")
ggplot(a,aes(x=condition,y=count))+
  geom_point()+
  geom_line()+
  theme_bw()
dev.off()  
 
res2<-cbind(rownames(res1),res1) 
write.xlsx(res2,file="RNA_DESeq_all.xlsx",rowNames = TRUE)


plotVoc <- function(diff){
  color <- rep("#999999",nrow(diff))
  
  color[diff$pvalue <0.05 & diff$log2FoldChange > 1] <- "#FC4E07"
  color[diff$pvalue <0.05 & diff$log2FoldChange < -1] <- "#00AFBB"
  
  par(oma = c(0,2,0,0))
  
  pdf("RNA_差异分析.pdf",width = 4,height=4)
  plot(diff$log2FoldChange,-log10(diff$pvalue),pch = 16,cex = 0.5,
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
  diff_label<-c("COL11A1","GJB2","KIF15","COL10A1","GPAM","NDN","ACVR1C","ADRA1A")
  geneList <- diff[rownames(diff) %in% diff_label,]
  geneList$label<-diff_label
  color = c()
  color[which(geneList$sig == "up")] = "black"
  color[which(geneList$sig == "down")] = "black"
  text(geneList$log2FoldChange,-log10(geneList$pvalue),
       labels = rownames(geneList),
       adj = c(0,1.5),
       cex = 0.4,
       col = color)
  dev.off()
}
plotVoc(res1)


