

cat("extract eQTL for trans CCVs from Colomics")
cat(".........................................")

## args
#args <- commandArgs(trailingOnly=T)
#rdata_file <- args[1]
#ccv_file <- args[2]

rdata_file <- "/path/to/eQTL-mQTL-Colomics/CCV_mQTLs_eQTLs/mQTLs_CCVSNPs_MN.Rdata"
ccv_file <- "/path/to/eQTL/trans-eQTL/ALL_INDEP_CCV_Trans.txt"

cat(" loading Rdata file rdata_file")
load(rdata_file)

cat(" loading trans CCVs")
ccv <- read.table(ccv_file,header=F,sep="\t",stringsAsFactors=F)
ccv <- ccv[,7]

ccv_mQTLs <- mQTLs[mQTLs$SNP %in% ccv,]
ccv_mQTLs$se <- ccv_mQTLs$beta/ccv_mQTLs$t.stat

ccv_mQTLs <- ccv_mQTLs[,c("SNP","CpG","GeneName","beta","se","t.stat","p.value")]
write.table(ccv_mQTLs,"/path/to/mQTL/trans-mQTL/result-colomics/Trans_in_Colonomics_MN_mQTL_ALL_output.txt",quote=F,sep="\t",row.names=F)