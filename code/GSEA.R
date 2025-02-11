rm(list=ls())
library(tidyr)
library(dplyr)
library(openxlsx)
library(ggplot2)
library(ggridges)
library(org.Hs.eg.db)
library(clusterProfiler)
library(gseaplot2)
library(enrichplot)
library(patchwork)
library(viridis)
library(viridisLite)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(GseaVis)

#Kegg需要ENTREZID
#GSEA需要基因的logFC
exp<-data[data$sig!="not",]

names(exp)[1]<-"SYMBOL"

eg = bitr(exp$SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

exp<-merge(exp,eg,by="SYMBOL")
geneList<-exp$logFC
names(geneList)<-exp$ENTREZID

geneList<-sort(geneList,decreasing = T)

kkgsea <- gseKEGG(geneList     = geneList,
                  organism     = 'hsa', 
                  minGSSize    = 10,
                  maxGSSize    = 500,
                  pvalueCutoff = 0.5,
                  pAdjustMethod = "none" )


gogsea <- gseGO(
  geneList,    # 根据logFC排序的基因集
  ont = "ALL",    # 可选"BP"、"MF"、"CC"三大类或"ALL"
  OrgDb = org.Hs.eg.db,    # 使用人的OrgDb
  keyType = "ENTREZID",    # 基因id类型
  pvalueCutoff = 1,      #如果是火山图这里就是1，如果不是那里就是0.05
  pAdjustMethod = "BH",    # p值校正方法
)

write.xlsx(gogsea,"差异分析_goGSEA.xlsx")
#########可视化

pdf("gsea_Phos.pdf",height = 8, width = 9)
ridgeplot(gogsea,10)
dev.off()

# 画图图
pdf("gsea_Phos_1.pdf",height = 8, width = 9)
dotplot(kkgsea, showCategory = 10, split = ".sign") + facet_grid(~.sign) +
  theme(plot.title = element_text(size = 10, color = "black", hjust = 0.5),#标题title
        axis.title = element_text(size =15,color = "black"), #GeneRatio
        axis.text = element_text(size = 10,color = "black"),
        axis.text.x = element_text(angle = 0, hjust = 1 ),
        legend.position = "right",
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        strip.text = element_text(size = 15))#分面字体大小
dev.off()

########GOGSEA#####的可视化分析

pdf("S100A2_表皮细胞分化GSEA.pdf",height = 4, width = 4)
gseaplot2(gogsea, title = gogsea$Description[2], geneSetID = 2)
dev.off()

pdf("S100A2_细胞代谢过程GSEA.pdf",height = 4, width = 4)
gseaplot2(gogsea, title = gogsea$Description[16], geneSetID = 16)
dev.off()

########GOGSEA多组可视化##########
pdf("Phos_多组.pdf",height = 8, width = 9)
ridgeplot(gogsea,showCategory = 15)
dev.off()


###############GOGSEA火山图##############
pdf("Phos_火山图.pdf",height = 6, width = 9)
volcanoGsea(data = gogsea,size = 4.5,nudge.y = c(-0.8,0.8))
dev.off()

#####关键基因热图######
genes <- as.data.frame(strsplit(gogsea$core_enrichment[2], "/"))
colnames(genes) <- 'gene'
ID1<-ID[ID$tissue=="Tumor",]
ID2<-rbind(ID1[ID1$sig=="high",],ID1[ID1$sig=="low",])
exp1 <- exp[exp$gene %in% genes$gene,]
exp2<-exp1[,colnames(exp) %in% ID2$patient_ID]
rownames(exp2)<-exp1$gene

expp <- t(scale(t(exp2)))

range(expp)
mycol <- colorRamp2(c(-1.1, 0, 2.5), c("#0571B0", "white", "#CA0020"))

#热图格子大小设置；
cellwidth = 0.5
cellheight = 0.5
cn = dim(expp)[2]
rn = dim(expp)[1]
w=cellwidth*cn
h=cellheight*rn

#热图绘制：
Heatmap(t(expp), #以其中30个基因展示为例
        name = 'Expression',
        col = mycol,
        cluster_rows = F)

#细节调整：
pdf("表皮细胞分化组的热图.pdf")
Heatmap(t(expp), #以其中30个基因展示为例
        name = 'Expression',
        col = mycol,
        rect_gp = gpar(col = "white", lwd = 1.5), #描边颜色和粗细
        column_names_gp = gpar(fontsize = 12, fontface = 'italic'), #列名/基因名斜体
        cluster_rows = F)
dev.off()
