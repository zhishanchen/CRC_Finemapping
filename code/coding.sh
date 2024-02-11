vep="/path/to/gene_score/coding/CCVs_VEP.txt" # vep annotation for CCVs

bed="/path/to/ALL_INDEP_CCV.txt_Uniq_CCV.txt.pos.txt.bed"

driver="/path/to/Cancer_driver_gene.list_CRC"
expr="/path/to/exp_gtex_colon_transverse_used_2_log2_norm.txt_expressed_gene_name"


# get ccvs which is precited to have high or moderate impacts on genes
comp2line.hash.pl -c 1 -q High_Moderate_impact.txt -d 4 -db $vep -e |grep -v nocommon |cut -f 2,3,5,6,7,31,32,46,677-679 | awk '!seen[$1,$2,$3,$4,$5]++'  > ${vep}.CODING.matrix


comp2line.hash.pl -c 5 -q ${vep}.CODING.matrix  -d 1 -db $driver -e |cut -f 1-12 > ${vep}.CODING.matrix.DriverGene.txt
comp2line.hash.pl -c 5 -q ${vep}.CODING.matrix.DriverGene.txt -d 1 -db $expr -last -e > ${vep}.CODING.matrix.DriverGene.txt.Expr.txt


perl -F"\t" -lane 'if(($F[2]=~/missense/ &&($F[5]=~/\bdeleterious\b/ || $F[6]=~/probably_damaging/)) || (($F[2]=~/stop/ || $F[2]=~/start/))){$score=1}else{$score=0};if($F[11]!~/nocommon/){$score2=$score+1}else{$score2=$score};if($F[12]=~/nocommon/){$score3=$score2*0.1;print "$_\t$score3"}else{$score3=$score2;print "$_\t$score3"}' ${vep}.CODING.matrix.DriverGene.txt.Expr.txt > ${vep}.CODING.SCORE.txt

comp2line.hash.pl -c 1 -q ${vep}.CODING.SCORE.txt -d 4 -db $bed -e > ${vep}.CODING.SCORE.txt.region.txt

awk '{if($14>=1)print $0}' ${vep}.CODING.SCORE.txt.region.txt > ${vep}.CODING.SCORE.txt.region.txt_HighConfidence_Gene.txt

cut -f 5,19,20 ${vep}.CODING.SCORE.txt.region.txt_HighConfidence_Gene.txt|sort |uniq -c |awk '{print $1"\t"$2"\t"$3"\t"$4}' > ${vep}.CODING.SCORE.txt.region.txt_HighConfidence_Gene.txt_CCV_num.txt

comp2line.hash.pl -model ncol -c 5,19 -q  ${vep}.CODING.SCORE.txt.region.txt_HighConfidence_Gene.txt -d 2,3 -db ${vep}.CODING.SCORE.txt.region.txt_HighConfidence_Gene.txt_CCV_num.txt -last -e > ${vep}.CODING.SCORE.txt.region.txt_HighConfidence_Gene.txt_CCV_num.txt_v2

Rscript aggregate.step2.R ${vep}
