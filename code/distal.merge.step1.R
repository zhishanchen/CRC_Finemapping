library(dplyr)
library(stringr)
args <- commandArgs(trailingOnly=T)
input <- args[1]
    
files <- list.files(pattern=".intersection$")
data<-list()
for (i in 1:length(files))
{
    data[[i]]<-read.table(files[i],header=F,sep="\t")
}

endb <-  data[[1]][,c(4,12)]
colnames(endb) <- c("ccv","gene")
endb <- endb %>% distinct()
endb <- as.data.frame(cbind(endb,rep("endb",dim(endb)[1])))
colnames(endb) <- c("ccv","gene","endb")

enatlas <- data[[2]][,c(4,12)]
colnames(enatlas) <- c("ccv","gene")
enatlas <- enatlas %>% distinct()
enatlas <- as.data.frame(cbind(enatlas,rep("enatlas",dim(enatlas)[1])))
colnames(enatlas) <- c("ccv","gene","enatlas")

fantom <- data[[3]][,c(4,12)]
colnames(fantom) <- c("ccv","gene")
fantom <- fantom %>% distinct()
fantom <- as.data.frame(cbind(fantom,rep("fantom5",dim(fantom)[1])))
colnames(fantom) <- c("ccv","gene","FANTOM5")

impet  <- data[[4]][,c(4,16)]
colnames(impet) <- c("ccv","gene")
impet <- impet[impet$gene %in% str_subset(impet$gene,"ENSG"),]
impet$gene <- gsub(",.*","",impet$gene,perl=TRUE)
impet <- impet %>% distinct()
impet <- as.data.frame(cbind(impet,rep("impet",dim(impet)[1])))
colnames(impet) <- c("ccv","gene","IMPET")

pmid1  <- data[[5]][,c(4,12)]
colnames(pmid1) <- c("ccv","gene")
pmid1$gene <- gsub(",.*","",pmid1$gene,perl=T)
pmid1$gene <- gsub("(\\S+)-\\d+","\\1",pmid1$gene,perl=T)
pmid1 <- pmid1 %>% distinct()
pmid1 <- as.data.frame(cbind(pmid1,rep("pmid1",dim(pmid1)[1])))
colnames(pmid1) <- c("ccv","gene","PMID1")


pmid2 <- data[[6]][,c(4,12)]
colnames(pmid2) <- c("ccv","gene")
pmid2 <- pmid2 %>% distinct()
pmid2 <- as.data.frame(cbind(pmid2,rep("pmid2",dim(pmid2)[1])))
colnames(pmid2) <- c("ccv","gene","PMID2")

pmid3 <- data[[7]][,c(4,12)]
colnames(pmid3) <- c("ccv","gene")
pmid3 <- pmid3 %>% distinct()
pmid3 <- as.data.frame(cbind(pmid3,rep("pmid3",dim(pmid3)[1])))
colnames(pmid3) <- c("ccv","gene","PMID3")

superen  <- data[[8]][,c(4,12)]
colnames(superen) <- c("ccv","gene")
superen <- superen %>% distinct()
superen <- as.data.frame(cbind(superen,rep("superen",dim(superen)[1])))
colnames(superen) <- c("ccv","gene","SUPEREN")

df <-  Reduce(function(d1, d2) merge(d1, d2, by=c("ccv","gene"),all=T), list(pmid1,pmid2,pmid3,endb,enatlas,fantom,impet,superen))
df <- df[,c(2,1,3,4,5,6,7,8,9,10)]

gf <- read.table(paste0(input,".genomic_feature_with_score.txt"),header=T,stringsAsFactors=F,sep="\t")
n <- dim(gf)[1]
data2 <- list()
sink(paste0(input,".genomic_feature_with_score.txt.v2"))
for(i in 1:n)
{
    sum <- (gf[i,5] + gf[i,7] + gf[i,9]) # 5: accessible_chromatin, 7:distal_gene_regulatory_elements,9:TFs
    if(sum == "0"){cat(gf[i,1],gf[i,2],"0","\n",sep="\t")}
    if(sum == "1"){cat(gf[i,1],gf[i,2],"1","\n",sep="\t")}
    if(sum > "1"){cat(gf[i,1],gf[i,2],"2","\n",sep="\t")}
}
sink()
gf <- read.table(paste0(input,".genomic_feature_with_score.txt.v2"),header=F,stringsAsFactors=F,sep="\t")[,1:3]
colnames(gf) <- c("ccv","position","genomic_feature")

df2 <- merge(df,gf,by="ccv",all.x=T)
df2 <- df2[,c(1:10,12)]

#anno <- read.table("gencode.v41lift37.annotation.gtf_gene.bed",header=F,stringsAsFactors=F,sep="\t")
#colnames(anno) <- c("chr","start","end","gene_id","gene","strand","type")

anno <- read.table("gencode.v19.annotation.gtf_gene.bed",header=F,stringsAsFactors=F,sep="\t")
colnames(anno) <- c("gene_id","gene","chr","start","end","anno","type","strand")


df2 <- merge(df2,anno,by="gene",all.x=T)
ids <- as.numeric(rownames(df2[duplicated(df2[,1:2]),])) ## remove genes which are annotated to multiple chromsomes,based on gencode 19
df3 <- df2[-ids,]
df3 <- df3 %>% filter(type=="protein_coding") ## only focus on protein coding genes

write.table(df3,paste0(input,".DISTAL.matrix"),row.names=F,quote=F,sep="\t")
