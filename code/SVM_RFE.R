set.seed(2023)
library(e1071)
#这里填写你存放的文件路径
source("msvmRFE.R")

input<-data
input$tissue <- ifelse(input$tissue == 'NAT', 0, 1)


library(e1071)
library(Rmpi)
library(snow)
library(parallel)

nfold = 10 #10倍交叉验证
nrows = nrow(input)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))

#make a cluster
cl <- makeMPIcluster(mpi.universe.size())

clusterExport(cl, list("input","svmRFE","getWeights","svm"))
results <-parLapply(cl,folds, svmRFE.wrap, input, k=10, halve.above=100)
top.features = WriteFeatures(results, input, save=F)

clusterExport(cl, list("top.features","results", "tune","tune.control"))
featsweep = parLapply(cl,1:100, FeatSweep.wrap, results, input)
stopCluster(cl)

no.info = min(prop.table(table(input[,1])))
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))

pdf("svm_rfe.pdf", height = 8, width = 10)
PlotErrors(errors, no.info=no.info)
dev.off()
plot(top.features)
mpi.exit()

