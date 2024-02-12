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
GWAS_SUM_A="/nobackup/sbcs/chenz27/fine-mapping6/gwas/07302022/CRC_consortium_ALL_risk_SNPs.txt_latest.bed.merge.region.gwas_ACCC.rsid"

# Genotype

VCF="/scratch/sbcs/chenzs/CRC_TWAS/genotype4/merge/TWAS.R03.vcf.gz"

# input - Asian CCVs
perl -F"\t" -lane '$keep=join"\t",@F[3..$#F];{print "$F[3]\t$F[0]\t$F[9]\t$keep"}' $CCV_DIR/All_ccv.txt_keep |grep "Asian" > $WORK_DIR/ALL_INDEP_CCV_Asian.txt

comp2line.hash.pl -c 7 -q $WORK_DIR/ALL_INDEP_CCV_Asian.txt -d 23 -db $GWAS_SUM_A -e > $WORK_DIR/ALL_INDEP_CCV_Asian.txt.tmp

perl -F"\t" -lane 'if($F[18]=~/(\d+):\d+/){print $1}' $WORK_DIR/ALL_INDEP_CCV_Asian.txt.tmp |sort |uniq  >  $WORK_DIR/chr.list

chrlist="$WORK_DIR/chr.list"

input="$WORK_DIR/ALL_INDEP_CCV_Asian.txt.tmp"

#rm $WORK_DIR/ALL_INDEP_CCV_Asian.txt.tmp

# eQTL
# extract genotype for SNPs

#echo "........................"
#echo " extract SNP genotype   "
#echo "........................"

#parallel -q echo perl -F\"\\t\" -lane \''if($F[18]=~/(\d+):(\d+)/){if($1=={}){print "$1\t$2"}}'\' $input \> $OUT_DIR/{}.snp :::: $chrlist | bash

#parallel -q echo tabix $VCF -R $OUT_DIR/{}.snp \> $OUT_DIR/{}.vcf.tmp :::: $chrlist |bash

#tabix -h $VCF $OUT_DIR/1.snp > $OUT_DIR/vcf.head

#parallel -q echo cat $OUT_DIR/vcf.head $OUT_DIR/{}.vcf.tmp  \> $OUT_DIR/{}.vcf :::: $chrlist |bash
#rm $OUT_DIR/*.vcf.tmp

## recode
#module load PLINK/1.9b_5.2
#parallel -q echo plink --vcf $OUT_DIR/{}.vcf --make-bed --const-fid --out $OUT_DIR/{} :::: $chrlist |bash
#parallel -q echo plink --bfile $OUT_DIR/{} --recode A --out $OUT_DIR/{} :::: $chrlist |bash
#parallel -q echo sed \''s/ /\t/g'\' $OUT_DIR/{}.raw \|sed \''s/NA//g'\' \> $OUT_DIR/{}.txt :::: $chrlist |bash

## for dosage
#parallel -q echo bcftools query -f \''%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%DS\t]\n'\' $OUT_DIR/{}.vcf \> $OUT_DIR/{}.genotype.tmp :::: $chrlist | bash
#sed -n '53p' $OUT_DIR/vcf.head  |sed 's/#CHROM/CHROM/' | sed 's/QUAL\tFILTER\tINFO\tFORMAT\t//' > $OUT_DIR/vcf.head.tmp
#mv $OUT_DIR/vcf.head.tmp $OUT_DIR/vcf.head
#parallel -q echo cat $OUT_DIR/vcf.head  $OUT_DIR/{}.genotype.tmp  \> $OUT_DIR/{}.genotype  :::: $chrlist |bash

echo "........................"
echo "run eQTL analysis       "
echo "........................"

module load GCC/8.2.0  CUDA/10.1.105  OpenMPI/3.1.4 R/3.6.0

parallel -q echo Rscript Asian_eQTL.R {} $OUT_DIR $input  :::: $chrlist |parallel -j 8 bash -c   

cat $OUT_DIR/ACCC_eQTL_chr*_output.txt |awk '!seen[$1,$2,$3,$4,$5]++'   >  $OUT_DIR/ACCC_eQTL_ALL_output.txt


