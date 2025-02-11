# 载入R包：
rm(list=ls())
library(ComplexHeatmap)
library(tidyverse)
library(circlize)
library(paletteer)


exp<-data[,3:24]
exp <- as.matrix(exp)
exp<- scale(exp, center = TRUE, scale = TRUE)

col_fun <- colorRamp2(c(-0.3, 0, 0.3), c("#5296cc", "#ffffff","#c7462b"))

#行分类
tissue<-data$tissue
row_anno <- rowAnnotation(
  Group = tissue,
  col = list(Group = c( "NAT" = "#B2DF8A", "Tumor"= "#FB9A99"))  # 为每个分类指定颜色
)


# colorbar的颜色：
d_palettes <- palettes_d_names
col_anno <- as.character(paletteer_d("ggsci::default_igv", n=22)) #随机选一个查看
names(col_anno) <- colnames(exp)

 # colorbar:列分类
top_anno <- HeatmapAnnotation(Module = colnames(exp),
                              col = list(Module = col_anno),
                              border = T,show_annotation_name = F,
                              # 去除图例：
                              show_legend = F)

# 绘图：
pdf("plots.pdf", height = 5, width = 5.5)
Heatmap(exp,
        # 行不聚类：
        cluster_rows = F,
        # 列聚类树高度：
        column_dend_height = unit(15, "mm"),
        # 设置颜色：
        col = col_fun,
        # 调整热图格子的边框颜色和粗细：
       # rect_gp = gpar(col = "white", lwd = 1.5),
        # 行名放左边：
        #row_names_side = "left",
        # 列名放上边：
        #column_names_side = "top",
        # colorbar:
        top_annotation = top_anno,
       left_annotation = row_anno,
        heatmap_legend_param = list(title = "Pearson Correlation"))

dev.off()







