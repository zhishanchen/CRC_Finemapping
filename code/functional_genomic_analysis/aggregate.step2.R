args <- commandArgs(trailingOnly=T)
input <- args[1]
type <- args[2] # DISTAL, PROMOTER, or CODING

inputfile1=paste0(input,".",type,".SCORE.txt.region.txt_HighConfidence_Gene.txt_CCV_num.txt_v2")
df <- read.table(inputfile1,header=T,sep="\t",stringsAsFactors=F)

if [[ type == "DISTAL" ]]; then
   colnames(df)[c(24,26,27)] <- c("count","indep","region")
fi 

if [[ type == "PROMOTER" ]]; then
   colnames(df)[c(4,5,17,18)] <- c("ccv","gene","indep","region")
fi

if [[ type == "CODING" ]]; then
   colnames(df)[c(1,5,19,20)] <- c("ccv","gene","indep","region")
fi

df2 <- aggregate(df[,1]~gene+indep+region, data = df,paste0,collapse=",")
colnames(df2) <- c("gene","indep","region","ccv")
write.table(df2,paste0(input,".",type,".SCORE.txt.region.txt_HighConfidence_Gene.txt_CCV_num.txt_v2_aggregate.txt"),row.names=F,quote=F,sep="\t")
