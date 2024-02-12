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
OUT_DIR="$WORK_DIR/result-accc"

echo "........................"
echo "run analysis in $WORK_DIR"
echo "........................"

# GWAS
GWAS_SUM_T="/nobackup/sbcs/chenz27/fine-mapping6/gwas/07302022/CRC_consortium_ALL_risk_SNPs.txt_latest.bed.merge.region.gwas_Trans.rsid"

# Genotype

VCF="/scratch/sbcs/chenzs/CRC_TWAS/genotype4/merge/TWAS.R03.vcf.gz"


# input - Trans CCVs

perl -F"\t" -lane '$keep=join"\t",@F[3..$#F];{print "$F[3]\t$F[0]\t$F[9]\t$keep"}' $CCV_DIR/All_ccv.txt_keep |grep "Trans-ancestry" > $WORK_DIR/ALL_INDEP_CCV_Trans.txt

comp2line.hash.pl -c 7 -q $WORK_DIR/ALL_INDEP_CCV_Trans.txt  -d 23 -db $GWAS_SUM_T -e > $WORK_DIR/ALL_INDEP_CCV_Trans.txt.tmp

perl -F"\t" -lane 'if($F[18]=~/(\d+):\d+/){print $1}' $WORK_DIR/ALL_INDEP_CCV_Trans.txt.tmp |sort |uniq  >  $WORK_DIR/chr.list

chrlist="$WORK_DIR/chr.list"

input="$WORK_DIR/ALL_INDEP_CCV_Trans.txt.tmp"

# eQTL
# extract genotype for SNPs
#parallel -q echo perl -F\"\\t\" -lane \''if($F[18]=~/(\d+):(\d+)/){if($1=={}){print "$1\t$2"}}'\' $input \> $OUT_DIR/{}.snp :::: $chrlist | bash

#parallel -q echo tabix $VCF -R $OUT_DIR/{}.snp \> $OUT_DIR/{}.vcf.tmp :::: $chrlist |bash

#tabix -h $VCF $OUT_DIR/1.snp > $OUT_DIR/vcf.head

#parallel -q echo cat $OUT_DIR/vcf.head $OUT_DIR/{}.vcf.tmp  \> $OUT_DIR/{}.vcf :::: $chrlist |bash

## recode
#module load PLINK/1.9b_5.2
#parallel -q echo plink --vcf $OUT_DIR/{}.vcf --make-bed --const-fid --out $OUT_DIR/{} :::: $chrlist |bash
#parallel -q echo plink --bfile $OUT_DIR/{} --recode A --out $OUT_DIR/{} :::: $chrlist |bash
#parallel -q echo sed \''s/ /\t/g'\' $OUT_DIR/{}.raw \|sed \''s/NA//g'\' \> $OUT_DIR/{}.txt :::: $chrlist |bash

echo "........................"
echo "run eQTL analysis       "
echo "........................"

#module load GCC/8.2.0  CUDA/10.1.105  OpenMPI/3.1.4 R/3.6.0

parallel -q echo Rscript Trans_eQTL.ACCC.R {} $OUT_DIR $input :::: $chrlist |parallel -j 8 bash -c

cat $OUT_DIR/Trans_in_ACCC_eQTL_chr*_output.txt | awk '!seen[$1,$2,$3,$4,$5]++' > $OUT_DIR/Trans_in_ACCC_eQTL_ALL_output.txt

## format output 
perl -F"\t" -lane 'if($F[0]=~/(\d+:\d+):\S+/){print "$_$1"}' $OUT_DIR/Trans_in_ACCC_eQTL_ALL_output.txt > $OUT_DIR/Trans_in_ACCC_eQTL_ALL_output.txt.tmp

comp2line.hash.pl -c 11 -q $OUT_DIR/Trans_in_ACCC_eQTL_ALL_output.txt.tmp -d 19 -db $input -e |perl -F"\t" -lane '$keep=join"\t",@F[1..9];{print "$F[44]\t$keep"}' > $OUT_DIR/Trans_in_ACCC_eQTL_ALL_output.txt

rm $OUT_DIR/Trans_in_ACCC_eQTL_ALL_output.txt.tmp

sed -i '1i SNP\tgene\tgene_name\tmc\tmaf\tbeta\tse\ttvalue\tp-value\talt_allele' $OUT_DIR/Trans_in_ACCC_eQTL_ALL_output.txt
