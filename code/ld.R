library("LDlinkR")
args <- commandArgs(trailingOnly=T)

input <- args[1] # fine-mapping region

isnp_trans <- read.table("Trans_indepSNP.txt",header=F,stringsAsFactors=F,sep="\t")
isnp_eur <- read.table("Eur_indepSNP.txt",header=F,stringsAsFactors=F,sep="\t")
isnp_asian <- read.table("Asian_indepSNP.txt",header=F,stringsAsFactors=F,sep="\t")

reportedSNP <- read.table("reported_SNP.txt",header=F,stringsAsFactors=F,sep="\t")

isnp_trans  <- as.data.frame(cbind(isnp_trans,rep("indep",dim(isnp_trans)[1])))
colnames(isnp_trans) <- c("region","signal","snp","population","tag")

isnp_eur  <- as.data.frame(cbind(isnp_eur,rep("indep",dim(isnp_eur)[1])))
colnames(isnp_eur) <- c("region","signal","snp","population","tag")

isnp_asian  <- as.data.frame(cbind(isnp_asian,rep("indep",dim(isnp_asian)[1])))
colnames(isnp_asian) <- c("region","signal","snp","population","tag")

reportedSNP <- as.data.frame(cbind(reportedSNP,rep("reported",dim(reportedSNP)[1])))
colnames(reportedSNP) <- c("region","snp","tag")


isnp_trans <- isnp_trans[isnp_trans$region %in% input,]
isnp_eur <- isnp_eur[isnp_eur$region %in% input,]
isnp_asian <- isnp_asian[isnp_asian$region %in% input,]
reportedsnp <- reportedSNP[reportedSNP$region %in% input,]

out_file_trans_vs_eur <- paste0("/path/to/LD/",input,"_indepSNP_trans_VS_eur.ld.txt")
out_file_trans_vs_asian <- paste0("/path/to/LD/",input,"_indepSNP_trans_VS_asian.ld.txt")
out_file_trans_vs_reported <- paste0("/path/to/LD/",input,"_indepSNP_trans_VS_reportedSNP.ld.txt")
out_file_eur_vs_reported <- paste0("/path/to/LD/",input,"_indepSNP_eur_VS_reportedSNP.ld.txt")
out_file_asian_vs_reported <- paste0("/path/to/LD/",input,"_indepSNP_asian_VS_reportedSNP.ld.txt")


len1 <- dim(isnp_trans)[1]
len2 <- dim(isnp_eur)[1]
len3 <- dim(isnp_asian)[1]
len4 <- dim(reportedsnp)[1]

sink(out_file_trans_vs_eur)
cat("region","\t","signal","\t","indep","\t","population","\t","singal","\t","indep","\t","population","\t","European_r2","\t","European_Dprime","\n")
for(i in 1:len1)
{
    var1 <- isnp_trans[i,3]
    tag1 <- as.character(isnp_trans[i,4])
    sig1 <- as.character(isnp_trans[i,2])
    for(j in 1:len2)
    {
        var2 <- isnp_eur[j,3]
        tag2 <- as.character(isnp_eur[j,4])
        sig2 <- as.character(isnp_eur[j,2])
        tryCatch({
            sum_eur <- LDpair(var1, var2, pop = "EUR", token = "your_token", output = "table")
            d_eur <- sum_eur[1,14]
            r2_eur <- sum_eur[1,15]
            cat(input,"\t", sig1,"\t",var1,"\t",tag1,"\t",sig2,"\t",var2,"\t",tag2,"\t",r2_eur,"\t",d_eur,"\n")
         },error=function(e){})
    }
}
sink()

sink(out_file_trans_vs_asian)
cat("region","\t","signal","\t","indep","\t","population","\t","singal","indep","\t","population","\t","Asian_r2","\t","Asian_Dprime","\n")
for(i in 1:len1)
{
    var1 <- isnp_trans[i,3]
    tag1 <- as.character(isnp_trans[i,4])
    sig1 <- as.character(isnp_trans[i,2])
    for(j in 1:len3)
    {
        var2 <- isnp_asian[j,3]
        tag2 <- as.character(isnp_asian[j,4])
        sig2 <- as.character(isnp_asian[j,2])
        tryCatch({
            sum_as <- LDpair(var1, var2, pop = "EAS", token = "your_token", output = "table")
            d_as <- sum_as[1,14]
            r2_as <- sum_as[1,15]
            pop_as <- sum_as[1,3]
            cat(input,"\t", sig1 ,"\t",var1,"\t",tag1,"\t",sig2,"\t",var2,"\t",tag2,"\t",r2_as,"\t",d_as,"\n")
            },error=function(e){})
    }
}
sink()

sink(out_file_trans_vs_reported)
cat("region","\t","signal","\t","indep","\t","population","\t","indep","\t","population","\t","European_r2","\t","European_Dprime","\t","Asian_r2","\t","Asian_Dprime","\n")
for(i in 1:len1)
{
    var1 <- isnp_trans[i,3]
    tag1 <- as.character(isnp_trans[i,4])
    sig <- as.character(isnp_trans[i,2])
    for(j in 1:len4)
    {
        var2 <- reportedsnp[j,2]
        tag2 <- as.character(reportedsnp[j,3])
        tryCatch({
            sum_eur <- LDpair(var1, var2, pop = "EUR", token = "your_token", output = "table")
            d_eur <- sum_eur[1,14]
            r2_eur <- sum_eur[1,15]
            pop_eur <- sum_eur[1,3]

            sum_as <- LDpair(var1, var2, pop = "EAS", token = "your_token", output = "table")
            d_as <- sum_as[1,14]
            r2_as <- sum_as[1,15]
            pop_as <- sum_as[1,3]

            cat(input,"\t",sig,"\t",var1,"\t",tag1,"\t",var2,"\t",tag2,"\t",r2_eur,"\t",d_eur,"\t",r2_as,"\t",d_as,"\n")
         },error=function(e){})
    }
}
sink()

sink(out_file_eur_vs_reported)
cat("region","\t","signal","\t","indep","\t","population","\t","indep","\t","population","\t","European_r2","\t","European_Dprime","\n")
for(i in 1:len2)
{
    var1 <- isnp_eur[i,3]
    tag1 <- as.character(isnp_eur[i,4])
    sig <- as.character(isnp_eur[i,2])
    for(j in 1:len4)
    {
        var2 <- reportedsnp[j,2]
        tag2 <- as.character(reportedsnp[j,3])
        tryCatch({
            sum_eur <- LDpair(var1, var2, pop = "EUR", token = "your_token", output = "table")
            d_eur <- sum_eur[1,14]
            r2_eur <- sum_eur[1,15]
            pop_eur <- sum_eur[1,3]
            cat(input,"\t",sig,"\t",var1,"\t",tag1,"\t",var2,"\t",tag2,"\t",r2_eur,"\t",d_eur,"\t","\n")
        },error=function(e){})
    }
}
sink()

sink(out_file_asian_vs_reported)
cat("region","\t","signal","\t","indep","\t","population","\t","indep","\t","population","\t","Asian_r2","\t","Asian_Dprime","\n")
for(i in 1:len3)
{
    var1 <- isnp_asian[i,3]
    tag1 <- as.character(isnp_asian[i,4])
    sig <- as.character(isnp_asian[i,2])
    for(j in 1:len4)
    {
        var2 <- reportedsnp[j,2]
        tag2 <- as.character(reportedsnp[j,3])
        tryCatch({
            sum_as <- LDpair(var1, var2, pop = "EAS", token = "your_token", output = "table")
            d_as <- sum_as[1,14]
            r2_as <- sum_as[1,15]
            pop_as <- sum_as[1,3]
            cat(input,"\t",sig,"\t",var1,"\t",tag1,"\t",var2,"\t",tag2,"\t","\t",r2_as,"\t",d_as,"\n")
        },error=function(e){})
    }
}
sink()
