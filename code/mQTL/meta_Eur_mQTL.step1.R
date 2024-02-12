library(dplyr)
colomics <- read.table("result-colomics/Eur_in_Colonomics_MN_mQTL_ALL_output.txt",header=T,sep="\t",row.names=NULL)
gtex <- read.table("result-gtex/Eur_in_GTEx_mQTL_ALL_output.txt.rsid",header=F,sep="\t",row.names=NULL)

colnames(colomics)[1] <- "rsid"
colnames(colomics)[2] <- "cpg_id"
colomics <- colomics[,c(1:5,7)]

colnames(gtex) <- c("cpg_id","variant_id","tss_distance","minor_allele_samples","mc","maf","pvalue","beta","se","var_pos","alt","rsid")
gtex <- gtex[,c(12,1,8:9,7,11)]

df_list <- list(colomics,gtex)

df <- Reduce(function(x, y) merge(x, y, by=c("rsid","cpg_id"),all=TRUE), df_list)  

colnames(df) <- c("rsid","cpg_id","GeneName","beta_col","se_col","p_col","beta_gtex","se_gtex","p_gtex","alt_gtex")

write.table(df,"Eur_ALL_mQTLs.txt",quote=F,sep="\t",row.names=F)



