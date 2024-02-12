cat("extract mQTL for trans CCVs from GTEx")
cat(".........................................")

args <- commandArgs(trailingOnly=T)
chrom <- args[1]
data <- paste0("/path/to/ColonTransverse.mQTLs.regular.txt_chr",chrom)
ccv_file <- "/path/to/mQTL/trans-mQTL/ALL_INDEP_CCV_Trans.txt.bed.hg38.v2"
output <- paste0("/path/to/mQTL/trans-mQTL/result-gtex/Trans_in_GTEx_mQTL_chr",chrom,"_output.txt")

cat(" loading data file")
mQTLs <- read.table(data,header=F,sep="\t",stringsAsFactors=F)

cat(" loading trans CCVs")
ccv <- read.table(ccv_file,header=F,sep="\t",stringsAsFactors=F)
ccv <- ccv[,c(7,11)]

ccv_mQTLs <- mQTLs[mQTLs[,10] %in% ccv[,2],]

write.table(ccv_mQTLs,output,quote=F,sep="\t",row.names=F)