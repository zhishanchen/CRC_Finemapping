# CRC_Finemapping
Code for the CRC finemapping paper


## eQTL analysis 
### trans-ancestry
```
bash Trans_eQTL.ACCC.sh
Rscript Trans_eQTL.Coloc.R
bash Trans_eQTL_GTEx.sh
Rscript meta_Trans.indep.step1.R
perl -F"\t" -lane 'if(($F[14] ne $F[10]) && ($F[14] ne "NA") && ($F[10] ne "NA")){$f=join"\t",@F[0..6];$e=join"\t",@F[8..$#F];$F[7]=($F[7]*-1);print "$f\t$F[7]\t$e"}else{print $_}' Trans_ALL_eQTLs.txt > Trans_ALL_eQTLs.txt.v2
Rscript meta_Trans.indep.step2.R

```
    