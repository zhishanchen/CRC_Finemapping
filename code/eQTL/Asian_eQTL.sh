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
OUT_DIR="$WORK_DIR/result"

echo "........................"
echo "run analysis in $WORK_DIR"
echo "........................"

# GWAS
GWAS_SUM_A="/path/to/CRC_consortium_ALL_risk_SNPs.txt_latest.bed.merge.region.gwas_ACCC.rsid"

# Genotype

VCF="/path/to/genotype.vcf.gz"

# input - Asian CCVs
perl -F"\t" -lane '$keep=join"\t",@F[3..$#F];{print "$F[3]\t$F[0]\t$F[9]\t$keep"}' $CCV_DIR/All_ccv.txt_keep |grep "Asian" > $WORK_DIR/ALL_INDEP_CCV_Asian.txt

comp2line.hash.pl -c 7 -q $WORK_DIR/ALL_INDEP_CCV_Asian.txt -d 23 -db $GWAS_SUM_A -e > $WORK_DIR/ALL_INDEP_CCV_Asian.txt.tmp

perl -F"\t" -lane 'if($F[18]=~/(\d+):\d+/){print $1}' $WORK_DIR/ALL_INDEP_CCV_Asian.txt.tmp |sort |uniq  >  $WORK_DIR/chr.list

chrlist="$WORK_DIR/chr.list"

input="$WORK_DIR/ALL_INDEP_CCV_Asian.txt.tmp"

echo "........................"
echo "run eQTL analysis       "
echo "........................"

module load GCC/8.2.0  CUDA/10.1.105  OpenMPI/3.1.4 R/3.6.0

parallel -q echo Rscript Asian_eQTL.R {} $OUT_DIR $input  :::: $chrlist |parallel -j 8 bash -c   

cat $OUT_DIR/ACCC_eQTL_chr*_output.txt |awk '!seen[$1,$2,$3,$4,$5]++'   >  $OUT_DIR/ACCC_eQTL_ALL_output.txt


