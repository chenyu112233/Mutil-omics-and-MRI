rm(list=ls())
library(readxl)
library(openxlsx)
library(ggplot2)
library(ggpubr)
library(ggExtra)

#
x<-as.numeric(con$S100A2)
y<-as.numeric(con$Keratinocytes)

df1=as.data.frame(cbind(x,y))
corT=cor.test(x,y,method="spearman")
cor=corT$estimate
pValue=corT$p.value
p1=ggplot(df1, aes(x, y)) + 
  xlab("S100A2_Protein")+ylab("Keratinocytes")+
  geom_point()+ geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
  stat_cor(method = 'spearman', aes(x =x, y =y))
p1 = p1 + theme(axis.title = element_text(size = 20))
p1 = p1 + theme(plot.title = element_text(size = 20))

p2=ggMarginal(p1, type = "density", xparams = list(fill = "orange"),yparams = list(fill = "blue"))


pdf(file="xcell_protein.pdf",width=5,height=5)
print(p2)
dev.off()

write.xlsx(con,"S100A2_keratinocyes.xlsx")
