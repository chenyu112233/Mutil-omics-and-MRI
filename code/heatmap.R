rm(list=ls())
library(pheatmap)

Tumor<-meta[meta$tissue=="Tumor",]
NAT<-meta[meta$tissue=="NAT",]
Tumor.dat<-data[,colnames(data) %in% Tumor$patient_ID]
NAT.dat<-data[,colnames(data) %in% NAT$patient_ID]

dat<-cbind(NAT.dat,Tumor.dat)
rownames(dat)<-data$gene

#归一化
#dat<-as.data.frame(log2(dat+1))

new_meta<-rbind(NAT,Tumor)

annotation_col1 = data.frame(
  Group = new_meta$tissue,
  row.names = new_meta$patient_ID
)


# 归一化函数
normalize_to_minus1_1 <- function(mat) {
  # 计算每行的最大值和最小值
  max_vals <- apply(mat, 1, max)
  min_vals <- apply(mat, 1, min)
  
  # 应用归一化公式
  norm_mat <- sweep(mat, 1, min_vals, FUN = "-")
  norm_mat <- sweep(norm_mat, 1, max_vals - min_vals, FUN = "/")
  norm_mat <- 2 * norm_mat - 1
  
  return(norm_mat)
}

# 对矩阵按行归一化到 [-1, 1]
dat <- normalize_to_minus1_1(dat)



p1<-pheatmap(
  dat, 
  #scale = 'row',  #基因表达值按行标准化，以对比不同组基因表达的差异
  cluster_rows = TRUE, cluster_cols = FALSE,  #行列是否聚类
  fontsize_row = 8, fontsize_col = 8,  #字体大小设置
  color = colorRampPalette(c('green', '#F7FCF5', 'red'))(100),  #热图颜色设置，本示例由绿到红渐变表示表达值增加
  annotation_col = annotation_col1,  #定义样本分组
 # cellwidth = 25, cellheight = 8, border_color = NA,  #格子宽度高度等设置
 show_colnames = FALSE,  # 不显示列标签
 show_rownames = FALSE,
 legend_breaks = c(-1,0,1),
 #legend_labels = c("low","heigh")
 #scaleFUN = normalize_to_minus1_1
)
ggsave("Phos_heatmap.pdf",plot = p1,width = 6, height = 6)

