library(dplyr)

accc <- read.table("result-accc/Trans_in_ACCC_eQTL_ALL_output.txt",header=T,sep="\t",row.names=NULL)
colomics <- read.table("result-colomics/Trans_in_Colonomics_MN_eQTL_ALL_output.txt",header=T,sep="\t",row.names=NULL)
gtex <- read.table("result-gtex/Trans_in_GTEx_eQTL_ALL_output.txt",header=T,sep="\t",row.names=NULL)
uva <- read.table("result-uva/Trans_in_UVA_eQTL_ALL_output.txt",header=T,sep="\t",row.names=NULL)

colnames(accc)[1] <- "rsid"
colnames(accc)[3] <- "Gene"
accc[,3] <- gsub("(ENSG\\d+).\\d+","\\1",accc[,3],perl=T)
accc <- accc[,c(1,3,7,8,10,11)]

colnames(colomics)[1] <- "rsid"
colnames(colomics)[2] <- "Gene"
colnames(colomics)[3] <- "Genename"
colomics <- colomics[,c(1:2,4:5,7)]

colnames(gtex)[1] <- "rsid"
colnames(gtex)[2] <- "Gene"
colnames(gtex)[3] <- "Genename"
gtex <- gtex[,c(1:2,6:9)]

colnames(uva)[1] <- "rsid"
colnames(uva)[3] <- "Gene"
colnames(uva)[4] <- "Genename"
uva[,3] <- gsub("(ENSG\\d+).\\d+","\\1",uva[,3],perl=T) 
uva <- uva[,c(1,3,4,7,8,10,11)]


df_list <- list(uva,accc,gtex, colomics)

df <- Reduce(function(x, y) merge(x, y, by=c("rsid","Gene"),all=TRUE), df_list)  

colnames(df) <- c("SNP","Gene","Genename","beta_uva","se_uva","p_uva","alt_uva","beta_accc","se_accc","p_accc","alt_accc","beta_gtex","se_gtex","p_gtex","alt_gtex","beta_col","se_col","p_col")

write.table(df,"Trans_ALL_eQTLs.txt",quote=F,sep="\t",row.names=F)



