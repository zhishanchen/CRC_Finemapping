library(dplyr)

trans_mqtl <- read.table("Eur_ALL_mQTLs.txt",header=T,stringsAsFactors=F,sep="\t")
df <- trans_mqtl
n <- dim(df)[1]

sink("Eur_ALL_mQTLs.txt_Meta.txt")
for(i in 1:n)
{
  
  if(!is.na(df[i,4]) && !is.na(df[i,7]))
     {
	   beta<-as.numeric(c(df[i,4],df[i,7]))
       	    p<-as.numeric(c(df[i,6],df[i,9]))
            N<-as.numeric(c(144,189))
            z<-qnorm(p/2)*sign(beta)
       	     tot_z<-sum(z*sqrt(N))/sqrt(sum(N))
       	      p_value<-2-2*pnorm(abs(tot_z))

        if(isTRUE(p_value==0)==TRUE){
            if(tot_z < 0){p_value <- 2*pnorm(q=(-tot_z), lower.tail=F)}
            else{
            p_value <- 2*pnorm(q=tot_z, lower.tail=F)
            }
         }
        cat(cat(c(unlist(df[i,]),tot_z,p_value,"col-gtex"),sep="\t"),"\n")
     }

    # col
    if(!is.na(df[i,4]) && is.na(df[i,7]))
    {
        cat(cat(c(unlist(df[i,]),df[i,4],df[i,6],"col"),sep="\t"),"\n")
    }

   # gtex
   if(is.na(df[i,4]) && !is.na(df[i,7]))
    {
        cat(cat(c(unlist(df[i,]),df[i,7],df[i,9],"gtex"),sep="\t"),"\n")
    }

}
sink()










