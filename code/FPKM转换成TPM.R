data<- read_excel("1_RNA_protein_MRI/data/tumor_RNA.xlsx",col_names=TRUE)
names<-data$gene
row.names(data)<-names
data<-data[-1]
#将fpkm转换为tpm
fpkmToTpm<-function(fpkm){
  exp(log(fpkm)-log(sum(fpkm))+log(1e6))
}
tpms<-apply(data,2,fpkmToTpm)
colSums(tpms)

write.xlsx(tpms,"1_RNA_protein_MRI/data/tumor_RNA_tpm.xlsx")
