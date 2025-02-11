rm(list=ls())
################
#MEs能够代表模块本身去跟性状进行计算相关性

##################

library(dplyr)
library(WGCNA)
library(FactoMineR)
library(factoextra)  
library(tidyverse)
library(readxl)

dia<-read.xlsx("data/clean_data/tumor_dia1.xlsx")
phos<-read.xlsx("data/clean_data/tumor_phos1.xlsx")
meta<-read.xlsx("data/meta/tumor_meta.xlsx")

dia.names<-dia[,1]
phos.names<-phos[,1]
dia<-dia[,-1]
phos<-phos[,-1]
row.names(dia)<-dia.names
row.names(phos)<-phos.names


datExpr<-phos
datExpr0=as.data.frame(t(datExpr))#行名是基因，列是样本

#设置标签,trait是标签，datExpr0是表达矩阵
pCR<-meta$pCR
Chemo<-ifelse(meta$Chemo=="TCb",0,1)
type<-ifelse(meta$type=="her2+",0,ifelse(meta$type=="TNBC",1,2))
HER2<-ifelse(meta$HER2=="Pos",1,0)
T<-ifelse(meta$T=="T1",1,ifelse(meta$T=="T2",2,3))
N<-ifelse(meta$N=="N0",0,ifelse(meta$N=="N1",1,ifelse(meta$N=="N2",2,3)))
Site<-ifelse(meta$Site=="L",1,0)
Ki.67<-meta$Ki.67
Age<-meta$Age
Grade<-ifelse(meta$Grade=="II",0,ifelse(meta$Grade=="III",1,2))
antiHER2<-ifelse(meta$antiHER2=="0",0,1)

trait<-cbind(pCR,Chemo,type,HER2,Site,T,N,Ki.67,Age,Grade,antiHER2)#标签数据
trait<-as.data.frame(trait)
rownames(trait)<-rownames(datExpr0)

#基础设置
options(stringsAsFactors = FALSE) 
enableWGCNAThreads() ## 打开多线程

#############1.1Sample cluster########### 
sampleTree = hclust(dist(datExpr0), method = "average") 
pdf(file = "1.聚类结果.pdf", width = 15, height = 8) 
par(cex = 0.6) 
par(mar = c(0,6,6,0)) 
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 2, 
     cex.axis = 1.5, cex.main = 2) 
### Plot a line to show the cut 
abline(h = 1.5e+09, col = "red")##剪切高度不确定，故无红线 
dev.off()

#######1.2聚类表达量和标签匹配##########
traitColors = numbers2colors(trait, signed = FALSE) 
pdf(file="1.2.聚类结果和标签.pdf",width=20,height=12) 
#画图
plotDendroAndColors(sampleTree, traitColors, 
                    groupLabels = names(trait), 
                    main = "Sample dendrogram and trait heatmap",
                    cex.colorLabels = 1.2, cex.dendroLabels = 1, cex.rowText = 2)
dev.off()

##########1.3聚类pca结果#######
dat.pca <- PCA(datExpr0, graph = F) 
pca <- fviz_pca_ind(dat.pca,
                    title = "Principal Component Analysis",
                    legend.title = "Groups",
                    geom.ind = c("point","text"), #"point","text"
                    pointsize = 2,
                    labelsize = 4,
                    repel = TRUE, #标签不重叠
                    col.ind = meta$pCR, # 分组上色
                    axes.linetype=NA,  # remove axeslines
                    mean.point=TRUE#去除分组中心点
) +
  theme(legend.position = "none")+  # "none" REMOVE legend
  coord_fixed(ratio = 1) #坐标轴的纵横比
pca
ggsave(pca,filename= "1.3 PCA analysis.pdf", width = 8, height = 8)



###阈值设置,选择power为后续的实验准备###
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")  #查看位于0.9以上的点，可以改变高度值
softPower=6#设置power需要自己修改                       
                        
##########2.1WGCNA构建网络和模块检测######
adjacency=adjacency(datExpr0,power=softPower)#相似矩阵
TOM=TOMsimilarity(adjacency)#两个矩阵的相似性
dissTOM=1-TOM#TOM的离散度
geneTree=hclust(as.dist(dissTOM),method="average")#层次聚类
#动态树切割和显示模块颜色
minModuleSize=30
dynamicMods=cutreeDynamic(dendro=geneTree,distM=dissTOM,
                          deepSplit=2,pamRespectsDendro=FALSE,
                          minClusterSize=minModuleSize)#切割

dynamicColors=labels2colors(dynamicMods)#转换为颜色
pdf(file="2.1 WGCNA Dynamic Tree Cut.pdf",width=8,height=6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()


#######2.2WGCNA合并模块，合并树ME：特征########
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)#计算模块的特征基因
MEs = MEList$eigengenes#提取特征基因
MEDiss = 1-cor(MEs)#离散度
METree = hclust(as.dist(MEDiss), method = "average")#层次聚类

pdf(file="2.2WGCNA合并模块树.pdf",width=7,height=6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.4  #高度
# Plot the cut line into the dendrogram
#abline(h=MEDissThres, col = "red")#画线
dev.off()


########2.3合并模块并且和绘图,tree和color#######
merge= mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs#特征基因

pdf(file="2.3 合并基因tree和color.pdf", width = 9, height = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()



#############3.1模块与表型数据关联#########
moduleColors = mergedColors
nGenes = ncol(datExpr0);
nSamples = nrow(datExpr0);
# 重新计算带有颜色标签的模块
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, trait, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
# 展示模块与表型数据的相关系数和 P值
pdf(file="3.1模块与表型数据关联热图.pdf",width=10,height=8)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(trait),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,#字体大小
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
colors = blueWhiteRed(50)
dev.off()



moduleTraitCor1<-t(moduleTraitCor)
data<-as.data.frame(moduleTraitCor1)
#表型和热图展示
library(circlize)
library(paletteer)
col_fun <- colorRamp2(c(-0.3, 0, 0.3), c("#5296cc", "#ffffff","#c7462b"))

# colorbar的颜色：
d_palettes <- palettes_d_names
col_anno <- as.character(paletteer_d("ggsci::default_igv", n=13)) #随机选一个查看
names(col_anno) <- colnames(data)

# colorbar:
top_anno <- HeatmapAnnotation(Module = colnames(data),
                              col = list(Module = col_anno),
                              border = T,show_annotation_name = F,
                              # 去除图例：
                              show_legend = F)

# 绘图：
pdf("plots.pdf", height = 5, width = 5.5)
Heatmap(data,
        # 行不聚类：
        cluster_rows = F,
        # 列聚类树高度：
        column_dend_height = unit(35, "mm"),
        # 设置颜色：
        col = col_fun,
        # 调整热图格子的边框颜色和粗细：
        rect_gp = gpar(col = "white", lwd = 1.5),
        # 行名放左边：
        row_names_side = "left",
        # 列名放上边：
        column_names_side = "top",
        # colorbar:
        top_annotation = top_anno,
        
        
        # 加星号：
        #cell_fun = function(j, i, x, y, width, height, fill) {
        # if (p_mat[i,j]) {
        # grid.text(sprintf("*"), x, y, gp = gpar(fontsize = 10))
        #  } else {
        #   grid.text(sprintf(""), x, y, gp = gpar(fontsize = 10))
        #  }
        #},
        ## 修改图例标题：
        heatmap_legend_param = list(title = "Pearson Correlation"))
dev.off()




################3.2可视化基因网络#############
nSelect = 400 
set.seed(10) 
select = sample(nGenes, size = nSelect) 
selectTOM = dissTOM[select, select] 
# 
selectTree = hclust(as.dist(selectTOM), method = "average") 
selectColors = moduleColors[select] 

plotDiss = selectTOM^7 
diag(plotDiss) = NA
library("gplots") 
pdf(file="3.2热图和原始模块树.pdf",width=9, height=9) 
mycol = colorpanel(250,'red','orange','lemonchiffon') 
TOMplot(plotDiss, selectTree, selectColors, col=mycol ,main = "Network heatmap plot, selected genes") 
dev.off()
#模块间关联性
pdf(file="3.2 合并树和热图.pdf", width=5, height=7.5) 
par(cex = 0.9) 
plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle= 90) 
dev.off()

#############4.1模块内###########
#提取模块基因
module = "green"
module_genes <- colnames(datExpr0)[moduleColors==module]
#GS：基因和性状相关性，MM：基因和模块相关性
modNames = substring(names(MEs), 3)

# 基因和模块的相关性及P值
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))#计算相关系数，模块和蛋白
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))#计算p值
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")#改变列名

# 基因和性状的相关性,需要修改合适的性状
oneTrait=pCR
oneTrait<-as.data.frame(oneTrait)
colnames(oneTrait)="pCR"
geneTraitSignificance = as.data.frame(cor(datExpr0, oneTrait, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
colnames(geneTraitSignificance) = paste("GS.", names(oneTrait), sep="")
colnames(GSPvalue) = paste("p.GS.", names(oneTrait), sep="")

column = match(module, modNames)#module在modNames中的位置
moduleGenes = moduleColors==module;

pdf(file="4.1单模块关系和蛋白.pdf", width=6, height=6) 
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in pCR", module, "module"),
                   ylab = "Protein significance for  ",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

tmp1 <- rownames(geneModuleMembership)[abs(geneModuleMembership[moduleGenes, column])>0.7]
tmp2 <- rownames(geneTraitSignificance)[abs(geneTraitSignificance[moduleGenes, 1])>0.04] 
genes <- unique(intersect(tmp1,tmp2))#tmp1和tmp2的交集
genes<-as.data.frame(genes)
colnames(genes)<-"protein_name"

write.csv(genes,"protein_greenModule.csv",quote = TRUE,row.names = FALSE)
write.csv(moduleColors,"moduleColors.csv",quote = TRUE,row.names = FALSE)
