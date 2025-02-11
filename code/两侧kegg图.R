rm(list = ls())
### R包载入：
library(org.Hs.eg.db)
library(clusterProfiler)
library(pathview)
library(topGO)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2) 
library(readxl)
library(openxlsx)
library(ggthemes)
library(ggplot2)
library(cols4all)
#将上下调的数据取出
exp<-data
colnames(exp)[1]<-"SYMBOL"
up<- as.character(exp[exp$sig=="up",]$SYMBOL)
down<- as.character(exp[exp$sig=="down",]$SYMBOL)
diff<-as.character(exp[exp$sig!="not",]$SYMBOL)

eg.up <- bitr(up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
eg.down <- bitr(down, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
eg.diff <- bitr(diff, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")


kk.up <- enrichKEGG(
  gene = eg.up$ENTREZID,
  keyType = "kegg",
  organism  = 'human',
  pvalueCutoff  = 0.01,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.01
)

kk.down <- enrichKEGG(
  gene = eg.down$ENTREZID,
  keyType = "kegg",
  organism  = 'human',
  pvalueCutoff  = 0.01,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.01
)



kk.diff <- enrichKEGG(
  gene = eg.diff$ENTREZID,
  keyType = "kegg",
  organism  = 'human',
  pvalueCutoff  = 0.01,
  pAdjustMethod  = "BH"
  #qvalueCutoff  = 0.9
)

#作为dataframe
kk.diff.dt<-as.data.frame(kk.diff)
kk.down.dt<-as.data.frame(kk.down)
kk.up.dt<-as.data.frame(kk.up)

down_kk<-kk.down.dt[kk.down.dt$pvalue<0.05,]
down_kk$group=-1
up_kk<-kk.up.dt[kk.up.dt$pvalue<0.05,]
up_kk$group=1

dat=rbind(up_kk,down_kk)
colnames(dat)
dat$p.adjust = -log10(dat$p.adjust)
dat$p.adjust=dat$p.adjust*dat$group 
dat=dat[order(dat$p.adjust,decreasing = F),]

#####绘图1######  
gk_plot <- ggplot(dat,aes(reorder(Description, p.adjust), y=p.adjust)) +
  geom_bar(aes(fill=factor((p.adjust>0)+1)),stat="identity", width=0.7, position=position_dodge(0.7)) +
  coord_flip() +
  scale_fill_manual(values=c("#0072B2", "#B20072"), guide=FALSE) +
  labs(x="", y="-log10(p.adjust)" ) +
  theme_pander()  +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.line.x = element_line(size = 0.3, colour = "black"),#x轴连线
        axis.ticks.length.x = unit(-0.20, "cm"),#修改x轴刻度的高度，负号表示向上
        axis.text.x = element_text(margin = margin(t = 0.3, unit = "cm")),##线与数字不要重叠
        axis.ticks.x = element_line(colour = "black",size = 0.3) ,#修改x轴刻度的线                         
        axis.ticks.y = element_blank(),
        axis.text.y  = element_text(hjust=0),
        panel.background = element_rect(fill=NULL, colour = 'white')
    )

pdf("phos_kegg.pdf",width = 10,height = 8)
gk_plot
dev.off()

######绘图2#####用dat绘图
down_kk$group="NAT"
up_kk$group="tumor"

up_kk1<-up_kk
down_kk1<-down_kk
down_kk1$LogP<- log(down_kk$pvalue)
up_kk1$LogP<- -log(up_kk$pvalue)

dat1=rbind(up_kk1,down_kk1)

#删掉重复的因子，如何没有重复可以不执行
table(duplicated(dat1$Description))
dup_row<-duplicated(dat1$Description)
dat1<-dat1[!dup_row,]

#指定因子；调整顺序：
dat1$Description <- factor(dat1$Description,levels = rev(dat1$Description))
#基础版上下调富集柱形图绘制：
p <- ggplot(dat1,
            aes(x =LogP, y = Description, fill = group)) + #数据映射
  geom_col() + #绘制添加条形图
  theme_bw()
p


#自定义主题调整：
mytheme <- theme(
  legend.position = 'none',
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  axis.line.x = element_line(color = 'grey60',size = 1.1),
  axis.text = element_text(size = 12)
)
p1 <- p + mytheme
p1
#先根据上下调标签拆分数据框：
up <- dat1[which(dat1$LogP > 0),]
down <- dat1[which(dat1$LogP < 0),]
#添加上调pathway标签：
p2 <- p1 +
  geom_text(data = up,
            aes(x = -0.2, y = Description, label = Description),
            size = 3.5,
            hjust = 1) #标签右对齐
p2
#添加下调pathway标签：
p3 <- p2 +
  geom_text(data = down,
            aes(x = 0.2, y = Description, label = Description),
            size = 3.5,
            hjust = 0) #标签左对齐
p3
#继续调整细节：
p4 <- p3 +
  scale_x_continuous(breaks=seq(-24, 24, 2)) + #x轴刻度修改
  labs(x = '-log(Pvalue)', y = ' ', title = 'Enriched Pathway') + #修改x/y轴标签、标题添加
  theme(plot.title = element_text(hjust = 0.5, size = 14)) #主标题居中、字号调整
p4
#颜色修改：

mycol1<-c("#4477AA","#EE6677" )
p5 <- p4 +
  scale_fill_manual(values = mycol1)
p5
#最后添加上下调提示标签：
p6 <- p5 +
  geom_text(x = 12, y = 15, label = "tumor", size = 6, color = '#EE6677') +
  geom_text(x = -12, y = 3, label = "NAT", size = 6, color = '#4477AA')
p6

pdf("phos_kegg_双侧图.pdf",width = 10,height = 8)
p6
dev.off()
