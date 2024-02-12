

cat("extract eQTL for trans CCVs from Colomics")
cat(".........................................")

## args
#args <- commandArgs(trailingOnly=T)
#rdata_file <- args[1]
#ccv_file <- args[2]

rdata_file <- "/nobackup/sbcs/chenz27/fine-mapping6/cond.07302022/eQTL/eQTL-mQTL-Colomics/CCV_mQTLs_eQTLs/eQTLs_CCVSNPs_MN.Rdata"
ccv_file <- "/nobackup/sbcs/chenz27/fine-mapping6/cond.07302022/eQTL/eur-eQTL/ALL_INDEP_CCV_Eur.txt"

cat(" loading Rdata file rdata_file")
load(rdata_file)

cat(" loading Trans CCVs")
ccv <- read.table(ccv_file,header=F,sep="\t",stringsAsFactors=F)
ccv <- ccv[,7]

ccv_eQTLs <- eQTLs[eQTLs$SNP %in% ccv,]
ccv_eQTLs$se <- ccv_eQTLs$beta/ccv_eQTLs$t.stat

ccv_eQTLs <- ccv_eQTLs[,c("SNP","ensembl","Gene","beta","se","t.stat","p.value")]
write.table(ccv_eQTLs,"/nobackup/sbcs/chenz27/fine-mapping6/cond.07302022/eQTL/eur-eQTL/result-colomics/Eur_in_Colonomics_MN_eQTL_ALL_output.txt",quote=F,sep="\t",row.names=F)