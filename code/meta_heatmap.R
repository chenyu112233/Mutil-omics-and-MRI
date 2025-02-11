rm(list=ls())

#准备工作
library(ComplexHeatmap)
library(RColorBrewer)

display.brewer.all(type = "all")
#选取喜欢的颜色编码
tempColors = brewer.pal(9, "PuBu") 
scales::show_col(tempColors,labels=T)


meta<-read_excel("./data/meta/meta_hotmap.xlsx")
data<-meta[,-1]


data<-as.data.frame(t(data))
data<-as.data.frame(t(data))
data$Age<-as.numeric(data$Age)#需要转换为数值才可以下面
data$Ki.67<-as.numeric(data$Ki.67)
#连续的值需要自己设置
ha<-HeatmapAnnotation(df=data,col = list(
 "pCR"=c("0"="#D0D1E6","1"="#023858"),
 "type"=c("Luminal B"="#74A9CF","her2+"="#A6BDDB","TNBC"="#045A8D"),
 "Chemo"=c("EC-T"="#0570B0","TCb"="#74A9CF","TY"="#D0D1E6","Others"="#225EA8"), 
  "HER2"=c("Pos"="#ECE7F2","Neg"="#807DBA"),
  "ER"=c("Pos"="#ECE7F2","Neg"="#807DBA"),
  "T"=c("T1"="#F7F7F7","T1b"="#D8DAEB","T1c"="#B2ABD2","T2"="#8073AC","T3"="#542788","T4"="#49006A"),
  "N"=c("N0"="#F7F7F7","N1"="#D8DAEB","N2"="#B2ABD2","N3"="#8073AC","N3c"="#542788"),
  "Site"=c("L"="#BCBDDC","R"="#6A51A3"),
  "Grade"=c("I"="#EFEDF5","II"="#BCBDDC","III"="#9E9AC8","IIII"="#807DBA"),
  "antiHER2"=c("0"="#BCBDDC","1"="#807DBA")
                                         
                                          ))#热图注释
zero_row_mat<-matrix(nrow=0,ncol=nrow(data))
vector=data
Hm<-Heatmap(zero_row_mat, name = "mat",
           # column_split = factor(rep("TCb",68)),
            top_annotation = ha,
            #rect_gp = gpar(col = "white", lwd = 1),#单元格的边界
            #column_title = "Meta",
            cluster_columns = TRUE
            )
# 设置热图的行和列的大小


pdf(file="./result/meta_picture/meta_heatmap1.pdf",width=10,height=5)
draw(Hm,merge_legend=TRUE,heatmap_legend_side="bottom",annotation_legend_side="right")
dev.off()

