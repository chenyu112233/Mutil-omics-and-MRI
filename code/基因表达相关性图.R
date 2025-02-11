rm(list=ls())
#前面两列是

cor.data <- read_excel("2_NAT_Project/result/bbj/RNA_Protein_SULT1E1.xlsx")
cor.data$RNA<-as.numeric(cor.data$RNA)
cor.data$Protein<-as.numeric(cor.data$Protein)
cor.data.Tumor<-cor.data[cor.data$tissue=="Tumor",]
cor.data.NAT<-cor.data[cor.data$tissue=="NAT",]

t.cor<-cor.test(cor.data.Tumor$RNA,cor.data.Tumor$Protein)
n.cor<-cor.test(cor.data.NAT$RNA,cor.data.NAT$Protein)
#我们观察到文献中点的大小还是变化的，所以计算散点图里点的大小，根据自己的数据来调整公式：
cor.data$diff <- abs(cor.data$RNA - cor.data$Protein) # 计算甲基化加速年龄
cor.data$size <- log10(cor.data$diff + 1) # 根据加速程度计算散点大小
# 计算图例里点的大小
cor.data$range<- cut(cor.data$size, breaks = quantile(cor.data$size), include.lowest = T) # 分配散点大小区间
cor.data$range2<-as.numeric(gsub("]","",sapply(strsplit(as.character(cor.data$range),","), "[",2), fixed = T)) # 取区间的后半部分，用于绘制图例


# 计算图例里需要绘制多少圆圈
num <- length(unique(cor.data$range2))
# 计算坐标轴的范围    
ylim <- range(cor.data$Protein) # y轴范围
xlim <- range(cor.data$RNA) # x轴范围
#必要的数据准备好了之后，接下来是生信花的个人show~
pdf("SULT1E1_RNA_Protein_表达.pdf", width = 4, height = 3.5)
par(bty="o", mgp = c(2,0.5,0), mar = c(4.1,4.1,2.1,4.1), tcl=-.25, font.main=3) # 画布基本设置
par(xpd=F) # 禁止显示超过画布的部分
plot(NULL, NULL, ylim = ylim, xlim = xlim, # 先绘制一个空的画布，仅有边框和坐标名
     xlab = "RNA ", ylab = "Protein",col="white",
     main = "")
rect(par("usr")[1], # 给画布设置背景色，掩盖边框
     par("usr")[3],
     par("usr")[2],
     par("usr")[4],
     col = "#EAE9E9",
     border = F)
grid(col = "white", lty = 1, lwd = 1.5) # 添加网格
## 画散点和回归线
# 在画布中添加一组（肿瘤组）的散点
tmp1 <- cor.data[which(cor.data$tissue == "Tumor"),]
reg1 <- lm(Protein~RNA, data=tmp1) # 计算回归线
points(tmp1$RNA, tmp1$Protein,
       pch = 19,    
       col = ggplot2::alpha("#E51718",0.8), # 重叠散点透明化
       cex = tmp1$size)
abline(reg1, lwd = 2, col = "#E51718") # 添加回归线
# 在画布中添加另一组（正常组）的散点
tmp2 <- cor.data[which(cor.data$tissue == "NAT"),]
reg2 <- lm(Protein~RNA, data=tmp2)
points(tmp2$RNA, tmp2$Protein,
       pch = 19,
       col = ggplot2::alpha("#1D2D60",0.8),
       cex = tmp2$size)
abline(reg2, lwd = 2, col = "#1D2D60")
# 如果有更多组，就按照以上格式依次添加。
## 画顶部和右侧地毯线
# 添加边际地毯线显示数据分布情况
rug(cor.data$RNA, col="black", lwd=1, side=3)
rug(cor.data$Protein, col="black", lwd=1, side=4)
## 添加相关性结果
# 手动把"~rho~" = 后面的数值修改为t.cor和n.cor的数值
text(1,17, adj = 0,expression("Tumour: N = 131; "~rho~" = 0.36; "~italic(P)~" < 0.001"), col = c("#E51718"), cex=0.5)
text(1,16,adj = 0,expression("Normal: N = 64; "~rho~" = 0.30; "~italic(P)~" = 0.013"), col = c("#1D2D60"), cex=0.5)
# 允许绘制超过画布的部分（用于添加图例）
par(xpd = T)
## 画图例
# 做散点的圆圈    
#points(rep(par("usr")[2]+0.5, num), y = seq(7,3,length.out = num),
#       pch = 19,
#       bty = "n",
#       cex = sort(unique(cor.data$range2)),
#      col = "black")
# 做散点图例的文字
#text(x = rep(par("usr")[2]+1, num + 1), y = c(10, seq(7,3,length.out = num)),
#     labels = c("Absolute\nVertical\nShift",
 #               round(10^(sort(unique(cor.data$range2))) - 1,0)), # 还原对数转化
#     adj = 0,cex = 0.8)
# 做分组的圆圈（肿瘤和正常）
points(
  x = rep(par("usr")[2]+0.5, 2), 
  y = c(16, 14.5),
       pch = 19,
       bty = "n",
       cex = 1,
       col = c("#E51718","#1D2D60"))
# 做分组图图例的文字
text(
  x = rep(par("usr")[2]+1, num + 1),
  y = c(16, 14.5),
     labels = c("Tumour","NAT"),
     adj = 0,cex = 0.8)
# 添加画布的边框
par(new = T, bty="o")
plot(-1, -1,
     col = "white",
     xlim = xlim, ylim = ylim,
     xlab = "", ylab = "",    
     xaxt = "n", yaxt = "n")
invisible(dev.off())

