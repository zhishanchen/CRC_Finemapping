library(dplyr)
eur_eqtl <- read.table("Eur_ALL_eQTLs.txt",header=T,stringsAsFactors=F,sep="\t")
df <- eur_eqtl
n <- dim(df)[1]

sink("Eur_ALL_eQTLs.txt_Meta.txt")
for(i in 1:n)
{

    # 1. all datasets
    if(!is.na(df[i,4]) && !is.na(df[i,8])  && !is.na(df[i,12]))
    {
        beta<-as.numeric(c(df[i,4],df[i,8],df[i,12]))
        p<-as.numeric(c(df[i,6],df[i,10],df[i,14]))
        N<-as.numeric(c(423,368,144))
        z<-qnorm(p/2)*sign(beta)
        tot_z<-sum(z*sqrt(N))/sqrt(sum(N))
        p_value<-2-2*pnorm(abs(tot_z))
	
        if(isTRUE(p_value==0)==TRUE){
            if(tot_z < 0){p_value <- 2*pnorm(q=(-tot_z), lower.tail=F)}
            else{
	    p_value <- 2*pnorm(q=tot_z, lower.tail=F)
	    }
	 }
        cat(cat(c(unlist(df[i,]),tot_z,p_value,"meta"),sep="\t"),"\n")
     }


  # 2. in gtex, and colonomics
  if(is.na(df[i,4]) && !is.na(df[i,8]) && !is.na(df[i,12]))
    {
        beta<-as.numeric(c(df[i,8],df[i,12]))
        p<-as.numeric(c(df[i,10],df[i,14]))
        N<-as.numeric(c(368,144))
        z<-qnorm(p/2)*sign(beta)
        tot_z<-sum(z*sqrt(N))/sqrt(sum(N))
        p_value<-2-2*pnorm(abs(tot_z))

        if(isTRUE(p_value==0)==TRUE){
            if(tot_z < 0){p_value <- 2*pnorm(q=(-tot_z), lower.tail=F)}
            else{
            p_value <- 2*pnorm(q=tot_z, lower.tail=F)
            }
         }
        cat(cat(c(unlist(df[i,]),tot_z,p_value,"gtex-col"),sep="\t"),"\n")
     }


  # 3. in uva, gtex
  if(!is.na(df[i,4]) && !is.na(df[i,8]) && is.na(df[i,12]))
    {
        beta<-as.numeric(c(df[i,4],df[i,8]))
        p<-as.numeric(c(df[i,6],df[i,10]))
        N<-as.numeric(c(423,368))
        z<-qnorm(p/2)*sign(beta)
        tot_z<-sum(z*sqrt(N))/sqrt(sum(N))
        p_value<-2-2*pnorm(abs(tot_z))

        if(isTRUE(p_value==0)==TRUE){
            if(tot_z < 0){p_value <- 2*pnorm(q=(-tot_z), lower.tail=F)}
            else{
            p_value <- 2*pnorm(q=tot_z, lower.tail=F)
            }
         }
        cat(cat(c(unlist(df[i,]),tot_z,p_value,"uva-gtex"),sep="\t"),"\n")
     }

  # 4. in uva,colonomics
  if(!is.na(df[i,4]) && is.na(df[i,8]) && !is.na(df[i,12]))
     {
	   beta<-as.numeric(c(df[i,4],df[i,12]))
       	    p<-as.numeric(c(df[i,6],df[i,14]))
            N<-as.numeric(c(423,144))
            z<-qnorm(p/2)*sign(beta)
       	     tot_z<-sum(z*sqrt(N))/sqrt(sum(N))
       	      p_value<-2-2*pnorm(abs(tot_z))

        if(isTRUE(p_value==0)==TRUE){
            if(tot_z < 0){p_value <- 2*pnorm(q=(-tot_z), lower.tail=F)}
            else{
            p_value <- 2*pnorm(q=tot_z, lower.tail=F)
            }
         }
        cat(cat(c(unlist(df[i,]),tot_z,p_value,"uva-col"),sep="\t"),"\n")
     }

    # uva 
    if(!is.na(df[i,4]) && is.na(df[i,8]) && is.na(df[i,12]))
    {
        cat(cat(c(unlist(df[i,]),df[i,4],df[i,6],"uva"),sep="\t"),"\n")
    }

   # gtex
   if(is.na(df[i,4]) && !is.na(df[i,8]) && is.na(df[i,12]))
    {
        cat(cat(c(unlist(df[i,]),df[i,8],df[i,10],"gtex"),sep="\t"),"\n")
    }

   # col

   if(is.na(df[i,4]) && is.na(df[i,8]) && !is.na(df[i,12]))
    {
        cat(cat(c(unlist(df[i,]),df[i,12],df[i,14],"col"),sep="\t"),"\n")
    }
}
sink()


############## determine significant gene ###########
df <- read.table("Eur_ALL_eQTLs.txt_Meta.txt",header=F,sep="\t")

n <- (dim(df)[2]-1)
y1 <- p.adjust(df[,n],method="bonferroni")
df2 <- as.data.frame(cbind(df,y1))
write.table(df2,"Eur_ALL_eQTLs.txt_Meta.txt.padjusted.txt",sep='\t', quote = F,row.names = F)







