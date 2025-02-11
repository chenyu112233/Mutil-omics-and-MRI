rm(list=ls())
library(org.Hs.eg.db)
library(openxlsx)
library(clusterProfiler)
#ENSEMBL
#ENSG00000067082.15
keytypes(org.Hs.eg.db)#查看多少种类型

old_id<-data$gene
data1<-data

tmp<-gsub("\\..*","",old_id)#将点去掉
data1$gene<-tmp



get<-c("SYMBOL","GENENAME","ENTREZID","ENSEMBL")
gene_symbols<-bitr(
  geneID=tmp,
  fromType ="ENSEMBL",
  toType = get,
  OrgDb =org.Hs.eg.db#输入类型，
)

#只有ENTREZID转换为SYMBOL需要
gene_symbols<-na.omit(gene_symbols)
#查看有多少重复的
table(duplicated(gene_symbols$SYMBOL))
dup_row<-duplicated(gene_symbols$SYMBOL)
gene_symbols<-gene_symbols[!dup_row,]

colnames(gene_symbols)[1]<-"gene"#准备合并
merged_data<-merge(gene_symbols,data1,by="gene",all.x=TRUE)



write.xlsx(merged_data,file="exp_RNA_4ID.xlsx",row.names = TRUE)
