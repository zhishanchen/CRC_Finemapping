
bed="/path/to/gene_score/promoter/ALL_INDEP_CCV.txt_Uniq_CCV.txt.pos.txt.bed"

promoter="/path/to/reference/Gene_anno/gencode.v19/hg19.upstream1000.downstream100.txt_with_name.txt"


driver="/path/to/Cancer_driver_gene.list_CRC"
expr="/path/to/exp_gtex_colon_transverse_used_2_log2_norm.txt_expressed_gene_name"

matrix="/path/to/score.txt.matrix.genomic_feature.txt"

bedtools intersect -a $bed -b $promoter -wo  |cut -f 1-4,14 |uniq  > $bed.CCV_in_Promoter.bed

cut -f 1-4,8,11 $matrix >  promoter_genomic_features.txt
feature="promoter_genomic_features.txt"

comp2line.hash.pl -c 4 -q $bed.CCV_in_Promoter.bed -d 1 -db $feature -e  > $bed.PROMOTER.matrix

comp2line.hash.pl -c 5 -q $bed.PROMOTER.matrix -d 1 -db $driver -e | cut -f 1-12 > $bed.PROMOTER.matrix.DriverGene.txt

comp2line.hash.pl -c 5 -q $bed.PROMOTER.matrix.DriverGene.txt -d 1 -db $expr -last -e > $bed.PROMOTER.matrix.DriverGene.txt.Expr.txt


perl -F"\t" -lane 'if($F[9]=~/1/){$ap=1}else{$ap=0}; if($F[10]=~/1/){$tf=1}else{$tf=0};if($F[11]!~/nocommon/){$dg=1}else{$dg=0};$score=$ap+$tf+$dg;if($F[12]=~/nocommon/){$score=$score*0.1;print "$_\t$score"}else{print "$_\t$score"}' $bed.PROMOTER.matrix.DriverGene.txt.Expr.txt > ${bed}.PROMOTER.SCORE.txt

awk '{if($14 > 1)print $0}' ${bed}.PROMOTER.SCORE.txt > ${bed}.PROMOTER.SCORE.txt.region.txt_HighConfidence_Gene.txt

# get the number of CCVs for each gene
cut -f 5,7,8 ${bed}.PROMOTER.SCORE.txt.region.txt_HighConfidence_Gene.txt |sort |uniq -c |awk '{print $1"\t"$2"\t"$3"\t"$4}' > ${bed}.PROMOTER.SCORE.txt.region.txt_HighConfidence_Gene.txt_CCV_num.txt


comp2line.hash.pl -model ncol -c 5,7 -q ${bed}.PROMOTER.SCORE.txt.region.txt_HighConfidence_Gene.txt  -d 2,3 -db ${bed}.PROMOTER.SCORE.txt.region.txt_HighConfidence_Gene.txt_CCV_num.txt -last -e > ${bed}.PROMOTER.SCORE.txt.region.txt_HighConfidence_Gene.txt_CCV_num.txt_v2

Rscript aggregate.step2.R ${bed}
