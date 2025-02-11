rm(list=ls())
library(dplyr)
library(openxlsx)
##ggcorrplot相关性图可视化
library(ggcorrplot)
library(ggtext)
library(RColorBrewer)
library(tidyverse)
library(ComplexHeatmap)
library(tidyverse)

dat<-data[,2:196]
dat<-as.data.frame(t(dat))

list<-data$gene

list[23:44] <- paste(list[23:44], "Protein", sep = "_")
list[1:22] <- paste(list[1:22], "RNA", sep = "_")

colnames(dat)<-data$gene

#计算相关性矩阵
corr_data <- round(cor(dat), 1)
write.xlsx(corr_data,"cor.xlsx",row.names=TRUE)
#计算对应的p值
corr_data1<-corr_data[1:22,23:44]


gg<-ggcorrplot(corr_data1,
           method = "circle",
           hc.order = TRUE,
           outline.color = "grey20") +
  scale_fill_gradientn(colors =  brewer.pal(11, "Spectral"),
                     name = NULL)+
  labs(x=NULL,y=NULL,
       title = "Correlation heat map")+
#  hrbrthemes::theme_ipsum() +
  theme(  
    plot.title = element_markdown(color = "black", size=18),
    legend.key.height = unit(1, "null"),
    legend.key.width = unit(0.5, "cm"),
    legend.frame = element_rect(color="black", linewidth = 0.25),
    plot.margin = margin(10, 10, 10, 10),
    plot.background = element_rect(fill = "white", color = "white"))

ggsave("correlation_heatmap.pdf", plot = gg, device = "pdf", width = 8, height = 8)
