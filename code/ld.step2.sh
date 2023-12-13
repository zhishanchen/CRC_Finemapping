
## keep indepedent signals from Trans-ancestral, and explore independent signals from either European or Asian analysis which are not in LD with that in Trans-ancestral analysis

cat LD2/region_*_indepSNP_trans_VS_eur.ld.txt |grep -v "European_Dprime"  |sed 's/ //g' > LD2/Trans_indepSNP_VS_Eur.txt
cat LD2/region_*_indepSNP_trans_VS_asian.ld.txt |grep -v "Asian_Dprime" |sed 's/ //g' > LD2/Trans_indepSNP_VS_Asian.txt

perl -F"\t" -lane 'if($F[6]>0.1){print $_}' LD2/Trans_indepSNP_VS_Eur.txt > Trans_indepSNP_LD_with_Eur.txt
perl -F"\t" -lane 'if($F[6]>0.1){print $_}' LD2/Trans_indepSNP_VS_Asian.txt > Trans_indepSNP_LD_with_Asian.txt

cut -f 6 Trans_indepSNP_LD_with_Eur.txt|sort |uniq  > tmp
grep -vf tmp -w Eur_indepSNP.txt > Eur_indepSNP.txt_noLD_with_Trans.txt
cut -f 6 Trans_indepSNP_LD_with_Asian.txt|sort |uniq > tmp
grep -vf tmp -w Asian_indepSNP.txt > Asian_indepSNP.txt_noLD_with_Trans.txt
rm tmp


## Explore LD between final independent signals with reported GWAS SNPs

cat LD2/region_*_indepSNP_trans_VS_reportedSNP.ld.txt |grep -v "r2_Eur"  |sed 's/ //g' > LD2/Trans_indepSNP_VS_reportedSNP.txt
cat LD2/region_*_indepSNP_eur_VS_reportedSNP.ld.txt |grep -v "r2_Eur" |sed 's/ //g' > LD2/Eur_indepSNP_VS_reportedSNP.txt
cat LD2/region_*_indepSNP_asian_VS_reportedSNP.ld.txt |grep -v "r2_Asian"  |sed 's/ //g' > LD2/Asian_indepSNP_VS_reportedSNP.txt

# >> trans
comp2line.hash.pl -c 3 -q Trans_indepSNP.txt -d 4 -db LD2/Trans_indepSNP_VS_reportedSNP.txt -e > Trans_indepSNP_VS_reportedSNP.ld.txt
# >> european
comp2line.hash.pl -c 3 -q Eur_indepSNP.txt_noLD_with_Trans.txt -d 4 -db LD2/Eur_indepSNP_VS_reportedSNP.txt -e > Eur_indepSNP.txt_noLD_with_Trans.txt_VS_reported.ld.txt
# >> asian
comp2line.hash.pl -c 3 -q Asian_indepSNP.txt_noLD_with_Trans.txt -d 4 -db LD2/Asian_indepSNP_VS_reportedSNP.txt -e > Asian_indepSNP.txt_noLD_with_Trans.txt_VS_reported.ld.txt


