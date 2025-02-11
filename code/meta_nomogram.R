rm(list=ls())
library(readxl)
library(ggplot2)

Mydata<-read_excel("./data/meta/tumor_meta1.xlsx")
Mydata<-Mydata[,-1]
names(Mydata)[11]<-"Ki.67"
names(Mydata)[10]<-"Site"

summary(Mydata)     ##初步查看我的数据

library(Hmisc) ##加载Hmisc包与rms包
library(rms)
dd<-datadist(Mydata)  #开始打包数据，注意Mydata改为自己取得名字
options(datadist="dd") #此处命令不用自己更改，直输

f_lrm <-lrm(pCR~Chemo+type+HER2_Status+ER_Status+T+N+Histologic_Grade+antiHER2+Age+Site+Ki.67
          ,data=Mydata)      #构建回归方程，使用lrm()函数构建二元LR，这里面的方程式就是结局~变量A+变量B+...，后面data来自自己命名的这个文件啦

summary(f_lrm)    #查看回归方程的结果

par(mgp=c(1.8,0.4,1),mar=c(5,5,3,1),bg="grey")     ##设置画布的命令，这里面的参数都可以调，不过我用的还算顺手，这几个数字不调也行哦，要是最后呈现的结果中有几个指标的轴线不好看，就适当更改下map里面的数字，随便调看看

nomogram <- nomogram(f_lrm,fun=function(x)1/(1+exp(-x)),   ##回归方程
                     fun.at = c(0.01,0.05,0.2,0.5,0.9,0.99), 
                     funlabel = "pCR",    ##风险轴刻度，结局是我的结局
                     conf.int = F,  #每个得分的置信度区间
                     abbrev = F)   #是否用简称代表因子变量

plot(nomogram)  #输出图片，右边plots可见列线图模型

pdf("nomogram.pdf", width = 10, height = 8)
plot(nomogram)
dev.off()
