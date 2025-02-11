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
# 读取文件，差异分析结果,只需要名字

exp<-data
DEG_list <- as.character(exp$gene)


#################################
# GO富集分析：
ego <- enrichGO(gene          = DEG_list,
                OrgDb         = org.Hs.eg.db,
                keyType = "SYMBOL",
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
)

# 转换为data.frame,方便查看：
ego_results<-as.data.frame(ego)

pdf("Protein.pdf",height = 7, width = 6)
## 最简单的方法绘制柱状图：直接barplot就可以搞定：
barplot(ego, showCategory=10, x = "GeneRatio")
dev.off()

# 保存富集分析结果：



pdf("GO_dotplot.pdf",height = 8, width = 9)
## 最简单的方法绘制柱状图：直接dotplot就可以搞定：
dotplot(ego,showCategory=10)
dev.off()

pdf("RNA_Protein_GO_dotplot_both_up.pdf",height = 8, width = 9)
dotplot(ego,
        x = "GeneRatio",
        color = "p.adjust",
        showCategory=10,
        split='ONTOLOGY',
        label_format = Inf)+#不换行
  #分面
  facet_grid(ONTOLOGY~.,
             space = 'free_y',#面板大小根据y轴自行调整
             scale='free_y'#子图坐标轴根据y轴自行调整
  )
dev.off()

#############把三种种类列出来###############
pdf("RNA_Protein_GO_网络_1.pdf",height = 10, width = 10)
enrichplot::cnetplot(ego, showCategory = 3,circular=TRUE, colorEdge = TRUE)#基因-通路关联网络图
dev.off()
#enrichplot::heatplot(GO,showCategory = 50)#基因-通路关联热图

pdf("RNA_Protein_GO_网络_2.pdf",height = 8, width = 10)
enrichplot::heatplot(ego,showCategory = 10)
dev.off()

#富集到的功能集/通路集之间的关联网络图：Jaccard系数,用来看两个样本中组成的相似性系数，showCategory是大家可以更改的部分，就是在图中显示多少个GO或者KEGG条目
pdf("RNA_Protein_GO_网络_3.pdf",height = 8, width = 9)
GO2 <- pairwise_termsim(ego)
enrichplot::emapplot(GO2,showCategory = 20, color = "p.adjust", layout = "kk")#通路间关联网络图method of calculating the similarity between nodes
dev.off()

#################################
# KEGG富集分析：
# 先转换ID：将gene名转换为ENTREZID
eg <- bitr(DEG_list, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

kegg <- enrichKEGG(
  gene = eg$ENTREZID,
  keyType = "kegg",
  organism  = 'human',
  pvalueCutoff  = 0.05,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.05
)

# 保存富集分析结果：


# 转换为data.frame,方便查看：
kegg_results<-as.data.frame(kegg)

pdf("KEGG_barplot.pdf",height = 8, width = 9)
## 最简单的方法绘制柱状图：直接barplot就可以搞定：
barplot(kegg, showCategory=10, x = "GeneRatio")
dev.off()


write.table(kegg, file = "kegg_Enrichment.txt",sep="\t", row.names =F, quote = F)



pdf("KEGG_dotplot.pdf",height = 8, width = 9)
## 最简单的方法绘制柱状图：直接dotplot就可以搞定：
dotplot(kegg,showCategory=10)
dev.off()



###########################GSEA################

tumor.dia.all<-read.table("~/result/normal_tumor/DEPs/2.1_tumor_dia_all.csv",sep=",",header=T)
tumor.phos.all<-read.table("~/result/normal_tumor/DEPs/2.1_tumor_phos_all.csv",sep=",",header=T)
normal.dia.all<-read.table("~/result/normal_tumor/DEPs/2.1_normal_dia_all.csv",sep=",",header=T)
normal.phos.all<-read.table("~/result/normal_tumor/DEPs/2.1_normal_phos_all.csv",sep=",",header=T)

exp<-tumor.phos.all
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
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "none" )
#########可视化

pdf("2.2.pdf",height = 8, width = 9)
ridgeplot(kkgsea,10)
dev.off()

    # 画图图
pdf("1.pdf",height = 8, width = 9)
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


gseaplot2(kkgsea,geneSetID = head(kkgsea@result$ID,12))
dotplot(kkgsea, showCategory = 15, color = "p.adjust")
