# input file
file="ALL_INDEP_CCV.txt_Uniq_CCV.txt.pos.txt"
input="/path/to/$file"

## ccv bed file
perl -F"\t"  -lane 'if($F[4]=~/(\d+):(\d+)/){print "chr$1\t$2\t$2\t$_"}' $input | sort -k1,1n -k2,2n  > ccv/gene_score/distal/$file.bed
# example line in input file: rs149685934     rs149685934     region_11       Trans-ancestry  1:222148219
bed="/gpfs52/nobackup/sbcs/chenz27/fine-mapping6/cond.07302022/ccv/gene_score/distal/$file.bed" 
# example line in bed:chr17   809425  809425  rs4968126       rs4968127       region_115      Trans-ancestry  17:809425


# functional genomic data
## experimental interaction
pmid1="/path/to/CRC_functional_genomic_data/interaction/PMID30224643/PMID30224643_E_P.bed"
pmid2="/path/to/PMID32841603/GSE133928_All_filt.intra.loop_counts.txt_filtered.bed_E_P.bed.final.v2"
pmid3="/path/to/CRC_functional_genomic_data/interaction/PMID34059508.interactions.fdr0.01_used.bed_E_P.bed.final.v3"
endb_experiment="/path/to/CRC_functional_genomic_data/Endb/ENdb_target_gene_colorectal.bed.v2"

## compuataional interaction
fantom5="/path/to/CRC_functional_genomic_data/interaction/Fantom5/enhancer_used_with_associated_gene.bed"
superenhancer="/path/to/CRC_functional_genomic_data/SuperEnhancer/super_enhancer_used.bed"
impet="/path/to/CRC_functional_genomic_data/interaction/4Dgenome/4DGenome_HomoSapiens_hg19.txt_in_CRC_cell.txt_enhancer.txt"
enhancerAtlas="/path/to/CRC_functional_genomic_data/EnhancerAtlas/enhancer_altas_used.bed"

## the overlap matrix of ccvs with genomic features
matrix="/path/to/score.txt.matrix.genomic_feature.txt"

## other information
tad="/path/to/CRC_functional_genomic_data/interaction/TAD/TAD_All.bed"
driver="/path/to/Cancer_driver_gene.list_CRC"
expr="/path/to/exp_gtex_colon_transverse_used_2_log2_norm.txt_expressed_gene_name"

# overlap between ccv and funtional genomic data
bedtools intersect -a $bed -b $pmid1 -wo > $bed.pmid1.intersection
bedtools intersect -a $bed -b $pmid2 -wo > $bed.pmid2.intersection
bedtools intersect -a $bed -b $pmid3 -wo > $bed.pmid3.intersection
bedtools intersect -a $bed -b $endb_experiment -wo > $bed.endb_experiment.intersection

bedtools intersect -a $bed -b $fantom5 -wo > $bed.FANTOM5.intersection
bedtools intersect -a $bed -b $superenhancer -wo > $bed.Superenhancer.intersection
bedtools intersect -a $bed -b $impet -wo > $bed.IM-PET.bed.intersection
bedtools intersect -a $bed -b $enhancerAtlas -wo > $bed.enhancerAtlas.bed.intersection

cut -f 1-8,11 $matrix > $bed.genomic_feature_with_score.txt

Rscript distal.merge.step1.R $bed

comp2line.hash.pl -c 1 -q $bed.DISTAL.matrix -d 1 -db $driver -e | cut -f 1-19 > $bed.DISTAL.matrix.DriverGene.txt
comp2line.hash.pl -c 1 -q $bed.DISTAL.matrix.DriverGene.txt -d 1 -db $expr -last -e > $bed.DISTAL.matrix.DriverGene.txt.Expr.txt

comp2line.hash.pl -c 4 -q $bed -d 2 -db $bed.DISTAL.matrix.DriverGene.txt.Expr.txt -e | awk '{if($9!~/nocommon/ && $20!~/NA/)print $0}' |perl -F"\t" -lane 'if($F[1] < $F[21]){print "$F[0]\t$F[1]\t$F[21]\t$_"};if($F[1] > $F[22]){print "$F[0]\t$F[22]\t$F[1]\t$_"}; if(($F[1] > $F[21]) && ($F[1] < $F[22])){print "$F[0]\t$F[21]\t$F[1]\t$_"}'  > $bed.DISTAL.matrix.DriverGene.txt.Expr.txt_for_TAD.txt

bedtools intersect -a ${bed}.DISTAL.matrix.DriverGene.txt.Expr.txt_for_TAD.txt -b $tad -f 1 -wo | cut -f 1-31 |awk '!seen[$12,$13]++' >${bed}.DISTAL.matrix.DriverGene.txt.Expr.txt_with_TAD.txt

Rscript distal.merge.step2.R $bed


perl -F"\t" -lane 'if($F[0]=~/ccv/){print "$_\tscore"}else{if($F[2]=~/NA/ && $F[3]=~/NA/ && $F[4]=~/NA/ && $F[5]=~/NA/){$expr=0}else{$expr=2};if($F[6]=~/NA/ && $F[7]=~/NA/ && $F[8]=~/NA/ && $F[9]=~/NA/){$comp=0}else{$comp=1};if($F[10]!~/1/ && $F[10]!~/2/){$gf=0};if($F[10]=~/1/ && $F[10]!~/2/){$gf=1};if($F[10]=~/2/){$gf=2};if($F[11]=~/nocommon/){$dg=0}else{$dg=1};$score=$expr+$comp+$gf+$dg;if($F[13]=~/tad/ && $F[12]!~/nocommon/){print "$_\t$score"};if($F[13]=~/tad/ && $F[12]=~/nocommon/){$score2=$score*0.1;print "$_\t$score2"};if($F[13]!~/tad/ && $F[12]!~/nocommon/){$score2=$score*0.05;print "$_\t$score2"};if($F[13]!~/tad/ && $F[12]=~/nocommon/){$score2=$score*0.05*0.1;print "$_\t$score2"}}' $bed.DISTAL.matrix.DriverGene.txt.Expr.txt.Gene > $bed.DISTAL.SCORE.txt

comp2line.hash.pl -c 1 -q ${bed}.DISTAL.SCORE.txt   -d 4 -db ${bed} -e > ${bed}.DISTAL.SCORE.txt.region.txt

# explore high confidence gene with score 5 or 6
awk '{if($15 > 4)print $0}' ${bed}.DISTAL.SCORE.txt.region.txt > ${bed}.DISTAL.SCORE.txt.region.txt_HighConfidence_Gene.txt

# 
sed '1d' ${bed}.DISTAL.SCORE.txt.region.txt_HighConfidence_Gene.txt |cut -f 2,20,21,22 |sort |uniq -c |awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' > ${bed}.DISTAL.SCORE.txt.region.txt_HighConfidence_Gene.txt_CCV_num.txt

comp2line.hash.pl -mode ncol -c 2,20 -q ${bed}.DISTAL.SCORE.txt.region.txt_HighConfidence_Gene.txt -d 2,3 -db ${bed}.DISTAL.SCORE.txt.region.txt_HighConfidence_Gene.txt_CCV_num.txt -e > ${bed}.DISTAL.SCORE.txt.region.txt_HighConfidence_Gene.txt_CCV_num.txt_v2

Rscript aggregate.step2.R $bed
