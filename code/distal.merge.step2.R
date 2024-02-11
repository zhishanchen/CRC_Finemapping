args <- commandArgs(trailingOnly=T)
input <- args[1]

inputfile1=paste0(input,".DISTAL.matrix.DriverGene.txt.Expr.txt_with_TAD.txt")
inputfile2=paste0(input,".DISTAL.matrix")

tad <- read.table(inputfile1,header=F,sep="\t")
tad2 <- as.data.frame(cbind(tad,rep("tad",dim(tad)[1])))
colnames(tad2)[12:13] <- c("gene","ccv")

df <- read.table(inputfile2,header=T,sep="\t")
df2 <- merge(df,tad2,by=c("ccv","gene"),all=T)
df3 <- df2[,c(1:11,46,47,48)]
colnames(df3)[12:14] <- c("driver_gene","expression","TAD")

write.table(df3,paste0(input,".DISTAL.matrix.DriverGene.txt.Expr.txt.Gene"),row.names=F,quote=F,sep="\t")
