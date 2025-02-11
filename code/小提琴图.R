rm(list=ls())
library(readxl)
library(ggplot2)
library(openxlsx)
library(limma)
library(ggrepel)
library(ggsignif)



p2 <- ggplot(data=data,aes(x=ER_Status,y=SLC44A4,colour = ER_Status))+ 
  geom_violin(#color = 'grey',
    alpha = 0.8,#alpha = 0.8 参数控制着小提琴图的透明度。具体来说，这里的 alpha 参数用于指定填充颜色的透明度，数值范围通常在 0 到 1 之间，其中 0 表示完全透明（即不可见），1 表示完全不透明。
    scale = 'width',#小提琴宽度
    #linewidth = 1, #外轮廓粗细
    trim = TRUE)+#trim = TRUE 参数控制着小提琴图的形状。当 trim = TRUE 时，小提琴图会根据数据的分布进行修剪
  geom_boxplot(mapping=aes(x=ER_Status,y=SLC44A4,colour=ER_Status,fill=ER_Status), #箱线图
               alpha = 0.5,
               size=1.5,
               width = 0.3)+ 
  geom_jitter(mapping=aes(x=ER_Status,y=SLC44A4,colour = ER_Status), #散点
              alpha = 0.3,size=3)+
  scale_fill_manual(limits=c("Pos","Neg"), 
                    values =c( "#e59698","#abd0cd"))+
  scale_color_manual(limits=c("Pos","Neg"), 
                     values=c("#e59698","#abd0cd"))+ #颜色
  geom_signif(mapping=aes(x=ER_Status,y=SLC44A4), # 不同组别的显著性
              comparisons = list(c("Pos", "Neg")), # 哪些组进行比较
              map_signif_level=T, # T显示显著性，F显示p value
              tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0), # 修改显著性线两端的长短
              y_position = 8, # 设置显著性线的位置高度
              size=1, # 修改线的粗细
              textsize = 4, # 修改显著性标记的大小
              test = "t.test")+ # 检验的类型,可以更改
  theme_bw()+#设置白色背景
  labs(x="ER Status",y="SLC44A4 log2(TPM+1)") +
  theme(axis.title.x = element_text(size = 17),
        axis.text = element_text(size = 15),# 设置x轴标签字体大小
        axis.title.y = element_text(size = 17))# 添加标题，x轴，y轴标签

p3 <- ggplot(data=data,aes(x=T,y=SLC44A4,colour = T))+ 
  geom_violin(#color = 'grey',
    alpha = 0.8,#alpha = 0.8 参数控制着小提琴图的透明度。具体来说，这里的 alpha 参数用于指定填充颜色的透明度，数值范围通常在 0 到 1 之间，其中 0 表示完全透明（即不可见），1 表示完全不透明。
    scale = 'width',#小提琴宽度
    #linewidth = 1, #外轮廓粗细
    trim = TRUE)+#trim = TRUE 参数控制着小提琴图的形状。当 trim = TRUE 时，小提琴图会根据数据的分布进行修剪
  geom_boxplot(mapping=aes(x=T,y=SLC44A4,colour=T,fill=T), #箱线图
               alpha = 0.5,
               size=1.5,
               width = 0.3)+ 
  geom_jitter(mapping=aes(x=T,y=SLC44A4,colour = T), #散点
              alpha = 0.3,size=3)+
  scale_fill_manual(limits=c("T1","T2","T3","T4"), 
                    values =c( "#e59698","#abd0cd","#f0b87f","#A6CEE3"))+
  scale_color_manual(limits=c("T1","T2","T3",'T4'), 
                     values=c("#e59698","#abd0cd","#f0b87f","#A6CEE3"))+ #颜色
  geom_signif(mapping=aes(x=T,y=SLC44A4), # 不同组别的显著性
              comparisons = list(c("T1", "T2"),
                                 c("T2","T3"),
                                 c("T3","T4")
                              ), # 哪些组进行比较
              map_signif_level=T, # T显示显著性，F显示p value
              tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0), # 修改显著性线两端的长短
              y_position = c(8,9,10), # 设置显著性线的位置高度
              size=1, # 修改线的粗细
              textsize = 4, # 修改显著性标记的大小
              test = "t.test")+ # 检验的类型,可以更改
  theme_bw()+#设置白色背景
  labs(x="T",y="SLC44A4 log2(TPM+1)") +
  theme(axis.title.x = element_text(size = 17),
        axis.text = element_text(size = 15),# 设置x轴标签字体大小
        axis.title.y = element_text(size = 17))# 添加标题，x轴，y轴标签
# 使用 plot_grid() 拼接图形

p4 <- ggplot(data=data,aes(x=N,y=SLC44A4,colour = N))+ 
  geom_violin(#color = 'grey',
    alpha = 0.8,#alpha = 0.8 参数控制着小提琴图的透明度。具体来说，这里的 alpha 参数用于指定填充颜色的透明度，数值范围通常在 0 到 1 之间，其中 0 表示完全透明（即不可见），1 表示完全不透明。
    scale = 'width',#小提琴宽度
    #linewidth = 1, #外轮廓粗细
    trim = TRUE)+#trim = TRUE 参数控制着小提琴图的形状。当 trim = TRUE 时，小提琴图会根据数据的分布进行修剪
  geom_boxplot(mapping=aes(x=N,y=SLC44A4,colour=N,fill=N), #箱线图
               alpha = 0.5,
               size=1.5,
               width = 0.3)+ 
  geom_jitter(mapping=aes(x=N,y=SLC44A4,colour = N), #散点
              alpha = 0.3,size=3)+
  scale_fill_manual(limits=c("N0","N1","N2","N3"), 
                    values =c( "#e59698","#abd0cd","#f0b87f","#A6CEE3"))+
  scale_color_manual(limits=c("N0","N1","N2",'N3'), 
                     values=c("#e59698","#abd0cd","#f0b87f","#A6CEE3"))+ #颜色
  geom_signif(mapping=aes(x=N,y=SLC44A4), # 不同组别的显著性
              comparisons = list(c("N0", "N1"),
                                 c("N1","N2"),
                                 c("N2","N3")
              ), # 哪些组进行比较
              map_signif_level=T, # T显示显著性，F显示p value
              tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0), # 修改显著性线两端的长短
              y_position = c(8,9,10), # 设置显著性线的位置高度
              size=1, # 修改线的粗细
              textsize = 4, # 修改显著性标记的大小
              test = "t.test")+ # 检验的类型,可以更改
  theme_bw()+#设置白色背景
  labs(x="N",y="SLC44A4 log2(TPM+1)") +
  theme(axis.title.x = element_text(size = 17),
        axis.text = element_text(size = 15),# 设置x轴标签字体大小
        axis.title.y = element_text(size = 17))# 添加标题，x轴，y轴标签

p5 <- ggplot(data=data,aes(x=Histologic_Grade,y=SLC44A4,colour = Histologic_Grade))+ 
  geom_violin(#color = 'grey',
    alpha = 0.8,#alpha = 0.8 参数控制着小提琴图的透明度。具体来说，这里的 alpha 参数用于指定填充颜色的透明度，数值范围通常在 0 到 1 之间，其中 0 表示完全透明（即不可见），1 表示完全不透明。
    scale = 'width',#小提琴宽度
    #linewidth = 1, #外轮廓粗细
    trim = TRUE)+#trim = TRUE 参数控制着小提琴图的形状。当 trim = TRUE 时，小提琴图会根据数据的分布进行修剪
  geom_boxplot(mapping=aes(x=Histologic_Grade,y=SLC44A4,colour=Histologic_Grade,fill=Histologic_Grade), #箱线图
               alpha = 0.5,
               size=1.5,
               width = 0.3)+ 
  geom_jitter(mapping=aes(x=Histologic_Grade,y=SLC44A4,colour = Histologic_Grade), #散点
              alpha = 0.3,size=3)+
  scale_fill_manual(limits=c("II","III"), 
                    values =c( "#e59698","#abd0cd"))+
  scale_color_manual(limits=c("II","III"), 
                     values=c("#e59698","#abd0cd"))+ #颜色
  geom_signif(mapping=aes(x=Histologic_Grade,y=SLC44A4), # 不同组别的显著性
              comparisons = list(c("II", "III")), # 哪些组进行比较
              map_signif_level=T, # T显示显著性，F显示p value
              tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0), # 修改显著性线两端的长短
              y_position = 10, # 设置显著性线的位置高度
              size=1, # 修改线的粗细
              textsize = 4, # 修改显著性标记的大小
              test = "t.test")+ # 检验的类型,可以更改
  theme_bw()+#设置白色背景
  labs(x="Histologic Grade",y="SLC44A4 log2(TPM+1)") +
  theme(axis.title.x = element_text(size = 17),
        axis.text = element_text(size = 15),# 设置x轴标签字体大小
        axis.title.y = element_text(size = 17))# 添加标题，x轴，y轴标签


p6 <- ggplot(data=data,aes(x=type,y=SLC44A4,colour = type))+ 
  geom_violin(#color = 'grey',
    alpha = 0.8,#alpha = 0.8 参数控制着小提琴图的透明度。具体来说，这里的 alpha 参数用于指定填充颜色的透明度，数值范围通常在 0 到 1 之间，其中 0 表示完全透明（即不可见），1 表示完全不透明。
    scale = 'width',#小提琴宽度
    #linewidth = 1, #外轮廓粗细
    trim = TRUE)+#trim = TRUE 参数控制着小提琴图的形状。当 trim = TRUE 时，小提琴图会根据数据的分布进行修剪
  geom_boxplot(mapping=aes(x=type,y=SLC44A4,colour=type,fill=type), #箱线图
               alpha = 0.5,
               size=1.5,
               width = 0.3)+ 
  geom_jitter(mapping=aes(x=type,y=SLC44A4,colour = type), #散点
              alpha = 0.3,size=3)+
  scale_fill_manual(limits=c("her2+","Luminal B","TNBC"), 
                    values =c( "#e59698","#abd0cd","#f0b87f"))+
  scale_color_manual(limits=c("her2+","Luminal B","TNBC"), 
                     values=c("#e59698","#abd0cd","#f0b87f"))+ #颜色
  geom_signif(mapping=aes(x=type,y=SLC44A4), # 不同组别的显著性
              comparisons = list(c("her2+", "Luminal B"),
                                 c("her2+","TNBC"),
                                 c("TNBC","Luminal B")
              ), # 哪些组进行比较
              map_signif_level=T, # T显示显著性，F显示p value
              tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0), # 修改显著性线两端的长短
              y_position = c(8,9,10), # 设置显著性线的位置高度
              size=1, # 修改线的粗细
              textsize = 4, # 修改显著性标记的大小
              test = "t.test")+ # 检验的类型,可以更改
  theme_bw()+#设置白色背景
  labs(x="type",y="SLC44A4 log2(TPM+1)") +
  theme(axis.title.x = element_text(size = 17),
        axis.text = element_text(size = 15),# 设置x轴标签字体大小
        axis.title.y = element_text(size = 17))# 添加标题，x轴，y轴标签

p7 <- ggplot(data=data,aes(x=pCR,y=SLC44A4,colour = pCR))+ 
  geom_violin(#color = 'grey',
    alpha = 0.8,#alpha = 0.8 参数控制着小提琴图的透明度。具体来说，这里的 alpha 参数用于指定填充颜色的透明度，数值范围通常在 0 到 1 之间，其中 0 表示完全透明（即不可见），1 表示完全不透明。
    scale = 'width',#小提琴宽度
    #linewidth = 1, #外轮廓粗细
    trim = TRUE)+#trim = TRUE 参数控制着小提琴图的形状。当 trim = TRUE 时，小提琴图会根据数据的分布进行修剪
  geom_boxplot(mapping=aes(x=pCR,y=SLC44A4,colour=pCR,fill=pCR), #箱线图
               alpha = 0.5,
               size=1.5,
               width = 0.3)+ 
  geom_jitter(mapping=aes(x=pCR,y=SLC44A4,colour = pCR), #散点
              alpha = 0.3,size=3)+
  scale_fill_manual(limits=c("RD","pCR"), 
                    values =c( "#e59698","#abd0cd"))+
  scale_color_manual(limits=c("RD","pCR"), 
                     values=c("#e59698","#abd0cd"))+ #颜色
  geom_signif(mapping=aes(x=pCR,y=SLC44A4), # 不同组别的显著性
              comparisons = list(c("RD", "pCR")), # 哪些组进行比较
              map_signif_level=T, # T显示显著性，F显示p value
              tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0), # 修改显著性线两端的长短
              y_position = 10, # 设置显著性线的位置高度
              size=1, # 修改线的粗细
              textsize = 4, # 修改显著性标记的大小
              test = "t.test")+ # 检验的类型,可以更改
  theme_bw()+#设置白色背景
  labs(x="Neoadjuvant chemotherapy response",y="SLC44A4 log2(TPM+1)") +
  theme(axis.title.x = element_text(size = 17),
        axis.text = element_text(size = 15),# 设置x轴标签字体大小
        axis.title.y = element_text(size = 17))# 添加标题，x轴，y轴标签

p8 <- ggplot(data=data,aes(x=antiHER2,y=SLC44A4,colour = antiHER2))+ 
  geom_violin(#color = 'grey',
    alpha = 0.8,#alpha = 0.8 参数控制着小提琴图的透明度。具体来说，这里的 alpha 参数用于指定填充颜色的透明度，数值范围通常在 0 到 1 之间，其中 0 表示完全透明（即不可见），1 表示完全不透明。
    scale = 'width',#小提琴宽度
    #linewidth = 1, #外轮廓粗细
    trim = TRUE)+#trim = TRUE 参数控制着小提琴图的形状。当 trim = TRUE 时，小提琴图会根据数据的分布进行修剪
  geom_boxplot(mapping=aes(x=antiHER2,y=SLC44A4,colour=antiHER2,fill=antiHER2), #箱线图
               alpha = 0.5,
               size=1.5,
               width = 0.3)+ 
  geom_jitter(mapping=aes(x=antiHER2,y=SLC44A4,colour = antiHER2), #散点
              alpha = 0.3,size=3)+
  scale_fill_manual(limits=c("0","1"), 
                    values =c( "#e59698","#abd0cd"))+
  scale_color_manual(limits=c("0","1"), 
                     values=c("#e59698","#abd0cd"))+ #颜色
  geom_signif(mapping=aes(x=antiHER2,y=SLC44A4), # 不同组别的显著性
              comparisons = list(c("0", "1")), # 哪些组进行比较
              map_signif_level=T, # T显示显著性，F显示p value
              tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0), # 修改显著性线两端的长短
              y_position = 10, # 设置显著性线的位置高度
              size=1, # 修改线的粗细
              textsize = 4, # 修改显著性标记的大小
              test = "t.test")+ # 检验的类型,可以更改
  theme_bw()+#设置白色背景
  labs(x="antiHER2",y="SLC44A4 log2(TPM+1)") +
  theme(axis.title.x = element_text(size = 17),
        axis.text = element_text(size = 15),# 设置x轴标签字体大小
        axis.title.y = element_text(size = 17))# 添加标题，x轴，y轴标签


pdf("S100A2_antiHER2",width = 5,height = 4)
p8# ncol 设置列数，这里设置为2列
dev.off()
