library(dplyr)


trans_eqtl <- read.table("Trans_ALL_eQTLs.txt.v2",header=T,stringsAsFactors=F,sep="\t")
df <- trans_eqtl
n <- dim(df)[1]

sink("Trans_ALL_eQTLs.txt.v2_Meta.txt")
for(i in 1:n)
{

    # 1. all datasets
    if(!is.na(df[i,4]) && !is.na(df[i,8])  && !is.na(df[i,12]) && !is.na(df[i,16]))
    {
        beta<-as.numeric(c(df[i,4],df[i,8],df[i,12],df[i,16]))
        p<-as.numeric(c(df[i,6],df[i,10],df[i,14],df[i,18]))
        N<-as.numeric(c(423,364,368,144))
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


  # 2. in accc, gtex, and colonomics
  if(is.na(df[i,4]) && !is.na(df[i,8]) && !is.na(df[i,12]) && !is.na(df[i,16]))
    {
        beta<-as.numeric(c(df[i,8],df[i,12],df[i,16]))
        p<-as.numeric(c(df[i,10],df[i,14],df[i,18]))
        N<-as.numeric(c(364,368,144))
        z<-qnorm(p/2)*sign(beta)
        tot_z<-sum(z*sqrt(N))/sqrt(sum(N))
        p_value<-2-2*pnorm(abs(tot_z))

        if(isTRUE(p_value==0)==TRUE){
            if(tot_z < 0){p_value <- 2*pnorm(q=(-tot_z), lower.tail=F)}
            else{
            p_value <- 2*pnorm(q=tot_z, lower.tail=F)
            }
         }
        cat(cat(c(unlist(df[i,]),tot_z,p_value,"accc-gtex-col"),sep="\t"),"\n")
     }


  # 3. in uva, gtex and colonomics
  if(!is.na(df[i,4]) &&  is.na(df[i,8]) && !is.na(df[i,12]) && !is.na(df[i,16]))
    {
        beta<-as.numeric(c(df[i,4],df[i,12],df[i,16]))
        p<-as.numeric(c(df[i,6],df[i,14],df[i,18]))
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
        cat(cat(c(unlist(df[i,]),tot_z,p_value,"uva-gtex-col"),sep="\t"),"\n")
     }

  # 4. in uva, accc, and colonomics
  if(!is.na(df[i,4]) && !is.na(df[i,8]) && is.na(df[i,12]) && !is.na(df[i,16]))
     {
	   beta<-as.numeric(c(df[i,4],df[i,8],df[i,16]))
       	    p<-as.numeric(c(df[i,6],df[i,10],df[i,18]))
            N<-as.numeric(c(423,364,144))
            z<-qnorm(p/2)*sign(beta)
       	     tot_z<-sum(z*sqrt(N))/sqrt(sum(N))
       	      p_value<-2-2*pnorm(abs(tot_z))

        if(isTRUE(p_value==0)==TRUE){
            if(tot_z < 0){p_value <- 2*pnorm(q=(-tot_z), lower.tail=F)}
            else{
            p_value <- 2*pnorm(q=tot_z, lower.tail=F)
            }
         }
        cat(cat(c(unlist(df[i,]),tot_z,p_value,"uva-accc-col"),sep="\t"),"\n")
     }

   # 5. in uva, accc and gtex
    if(!is.na(df[i,4]) && !is.na(df[i,8]) && !is.na(df[i,12]) && is.na(df[i,16]))
    {
        beta<-as.numeric(c(df[i,4],df[i,8],df[i,12]))
        p<-as.numeric(c(df[i,6],df[i,10],df[i,14]))
        N<-as.numeric(c(423,364,368))
        z<-qnorm(p/2)*sign(beta)
        tot_z<-sum(z*sqrt(N))/sqrt(sum(N))
        p_value<-2-2*pnorm(abs(tot_z))

        if(isTRUE(p_value==0)==TRUE){
            if(tot_z < 0){p_value <- 2*pnorm(q=(-tot_z), lower.tail=F)}
            else{
            p_value <- 2*pnorm(q=tot_z, lower.tail=F)
            }
         }
        cat(cat(c(unlist(df[i,]),tot_z,p_value,"uva-accc-gtex"),sep="\t"),"\n")
     }


    # 6. in uva and accc
    if(!is.na(df[i,4] != "NA") && !is.na(df[i,8] != "NA") && is.na(df[i,12]) && is.na(df[i,16]))
    {
        beta<-as.numeric(c(df[i,4],df[i,8]))
        p<-as.numeric(c(df[i,6],df[i,10]))
        N<-as.numeric(c(423,364))
        z<-qnorm(p/2)*sign(beta)
        tot_z<-sum(z*sqrt(N))/sqrt(sum(N))
        p_value<-2-2*pnorm(abs(tot_z))

        if(isTRUE(p_value==0)==TRUE){
            if(tot_z < 0){p_value <- 2*pnorm(q=(-tot_z), lower.tail=F)}
            else{
            p_value <- 2*pnorm(q=tot_z, lower.tail=F)
            }
         }
        cat(cat(c(unlist(df[i,]),tot_z,p_value,"uva-accc"),sep="\t"),"\n")
     }    


    # 7. in uva and gtex
    if(!is.na(df[i,4]) && is.na(df[i,8]) && !is.na(df[i,12]) && is.na(df[i,16]))
    {
        beta<-as.numeric(c(df[i,4],df[i,12]))
        p<-as.numeric(c(df[i,6],df[i,14]))
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


    # 8. in uva and colonomics

    if(!is.na(df[i,4]) && is.na(df[i,8]) && is.na(df[i,12]) && !is.na(df[i,16]))
    {
        beta<-as.numeric(c(df[i,4],df[i,16]))
        p<-as.numeric(c(df[i,6],df[i,18]))
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


     # 9. in accc and gtex

    if(is.na(df[i,4]) && !is.na(df[i,8]) && !is.na(df[i,12]) && is.na(df[i,16]))
    {
        beta<-as.numeric(c(df[i,8],df[i,12]))
        p<-as.numeric(c(df[i,10],df[i,14]))
        N<-as.numeric(c(364,368))
        z<-qnorm(p/2)*sign(beta)
        tot_z<-sum(z*sqrt(N))/sqrt(sum(N))
        p_value<-2-2*pnorm(abs(tot_z))

        if(isTRUE(p_value==0)==TRUE){
            if(tot_z < 0){p_value <- 2*pnorm(q=(-tot_z), lower.tail=F)}
            else{
            p_value <- 2*pnorm(q=tot_z, lower.tail=F)
            }
         }
        cat(cat(c(unlist(df[i,]),tot_z,p_value,"accc-gtex"),sep="\t"),"\n")
     }

     # 10. in accc and colonomics

    if(is.na(df[i,4]) && !is.na(df[i,8]) && is.na(df[i,12]) && !is.na(df[i,16]))
    {
        beta<-as.numeric(c(df[i,8],df[i,16]))
        p<-as.numeric(c(df[i,10],df[i,18]))
        N<-as.numeric(c(364,144))
        z<-qnorm(p/2)*sign(beta)
        tot_z<-sum(z*sqrt(N))/sqrt(sum(N))
        p_value<-2-2*pnorm(abs(tot_z))

        if(isTRUE(p_value==0)==TRUE){
            if(tot_z < 0){p_value <- 2*pnorm(q=(-tot_z), lower.tail=F)}
            else{
            p_value <- 2*pnorm(q=tot_z, lower.tail=F)
            }
         }
        cat(cat(c(unlist(df[i,]),tot_z,p_value,"accc-col"),sep="\t"),"\n")
     }


     # 11. in gtex and col

    if(is.na(df[i,4]) && is.na(df[i,8]) && !is.na(df[i,12]) && !is.na(df[i,16]))
    {
        beta<-as.numeric(c(df[i,12],df[i,16]))
        p<-as.numeric(c(df[i,14],df[i,18]))
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

    # uva 
    if(!is.na(df[i,4]) && is.na(df[i,8]) && is.na(df[i,12]) && is.na(df[i,16]))
    {
        cat(cat(c(unlist(df[i,]),df[i,4],df[i,6],"uva"),sep="\t"),"\n")
    }

   # accc
   if(is.na(df[i,4]) && !is.na(df[i,8]) && is.na(df[i,12]) && is.na(df[i,16]))
    {
        cat(cat(c(unlist(df[i,]),df[i,8],df[i,10],"accc"),sep="\t"),"\n")
    }

   # gtex

   if(is.na(df[i,4]) && is.na(df[i,8]) && !is.na(df[i,12]) && is.na(df[i,16]))
    {
        cat(cat(c(unlist(df[i,]),df[i,12],df[i,14],"gtex"),sep="\t"),"\n")
    }

   # col

   if(is.na(df[i,4]) && is.na(df[i,8]) && is.na(df[i,12]) && !is.na(df[i,16]))
    {
        cat(cat(c(unlist(df[i,]),df[i,16],df[i,18],"col"),sep="\t"),"\n")
    }

}
sink()










