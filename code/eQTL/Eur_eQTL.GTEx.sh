#!/bin/bash

echo "To perform eQTL analysis for CCVs from fine-mapping"

DATE=`date +%Y-%m-%d`
echo $DATE

###################
# setting
###################

WORK_DIR=`pwd|perl -F"/" -lane '$p=join"/",@F[0..($#F)];print $p'`
COND_DIR=`pwd|perl -F"/" -lane '$p=join"/",@F[0..($#F-2)];print $p'`
CCV_DIR="$COND_DIR/ccv" # the folder with ccv
OUT_DIR="$WORK_DIR/result-gtex"

echo "........................"
echo "run analysis in $WORK_DIR"
echo "........................"

# GWAS
GWAS_SUM_A="/path/to/CRC_consortium_ALL_risk_SNPs.txt_latest.bed.merge.region.gwas_ACCC.rsid"
GWAS_SUM_E="/path/to/CRC_consortium_ALL_risk_SNPs.txt_latest.bed.merge.region.gwas_Eur.rsid"
GWAS_SUM_T="/path/to/CRC_consortium_ALL_risk_SNPs.txt_latest.bed.merge.region.gwas_Trans.rsid"

# input - European CCVs

perl -F"\t" -lane '$keep=join"\t",@F[3..$#F];{print "$F[3]\t$F[0]\t$F[9]\t$keep"}' $CCV_DIR/All_ccv.txt_keep |grep "European" > $WORK_DIR/ALL_INDEP_CCV_Eur.txt

comp2line.hash.pl -c 7 -q $WORK_DIR/ALL_INDEP_CCV_Eur.txt -d 23 -db $GWAS_SUM_E -e |cut -f 1-11,19 |perl -F"\t" -lane 'if($F[11]=~/(\d+):(\d+)/){print "chr$1\t$2\t$2\t$_"}' > $WORK_DIR/ALL_INDEP_CCV_Eur.txt.bed.tmp

cut -f 1-3,7-10,13,14,15 $WORK_DIR/ALL_INDEP_CCV_Eur.txt.bed.tmp > $WORK_DIR/ALL_INDEP_CCV_Eur.txt.bed

CrossMap.py bed $WORK_DIR/hg19ToHg38.over.chain.gz $WORK_DIR/ALL_INDEP_CCV_Eur.txt.bed  $WORK_DIR/ALL_INDEP_CCV_Eur.txt.bed.hg38

rm $WORK_DIR/ALL_INDEP_CCV_Eur.txt.bed.tmp

echo "........................"
echo " extract eQTL from GTEx   "
echo "........................"

perl $WORK_DIR/eQTL.GTExV8.pl $WORK_DIR/ALL_INDEP_CCV_Eur.txt.bed.hg38 $OUT_DIR/ALL_INDEP_CCV_Eur.txt.bed.hg38_eQTL_GTEx8.txt.tmp

perl -F"\t" -lane '{print "$F[6]\t$F[10]\t$F[11]\t$F[12]\t$F[13]\t$F[15]\t$F[16]\t$F[14]\t$F[17]"}' $OUT_DIR/ALL_INDEP_CCV_Eur.txt.bed.hg38_eQTL_GTEx8.txt.tmp > $OUT_DIR/Eur_in_GTEx_eQTL_ALL_output.txt 
