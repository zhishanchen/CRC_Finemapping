module load PLINK/2.00-alpha1

region=$1

chr=`grep "$region\b" /path/to/any_file_with_chromsome_information | cut -f 1`

trans="Trans_indepSNP.txt"
eur="Eur_indepSNP.txt"
asian="Asian_indepSNP.txt" # region_106      signal_1        rs68097734      Asian
report="reported_SNP.txt" # region_1        rs2807367

echo "---------------------------------------------"
trans_num=`grep "$region\b" $trans |wc -l |awk '{print $1}'`
echo "--$trans_num indepSNPs in $region from Trans" > LD2/${region}.log 
eur_num=`grep "$region\b" $eur |wc -l |awk '{print $1}'`
echo "--$eur_num indepSNPs in $region from Eur" >> LD2/${region}.log
asian_num=`grep "$region\b" $asian |wc -l |awk '{print $1}'`
echo "--$asian_num indepSNPs in $region from Asian" >> LD2/${region}.log 
report_num=`grep "$region\b" $report |wc -l |awk '{print $1}'`
echo "--$report_num SNPs in $region from Report" >> LD2/${region}.log 
echo "---------------------------------------------"

SNP_TRANS=(`grep "$region\b" Trans_indepSNP.txt |cut -f 3`)
SNP_EUR=(`grep "$region\b" Eur_indepSNP.txt |cut -f 3`)
SNP_ASIAN=(`grep "$region\b" Asian_indepSNP.txt |cut -f 3`)
SNP_REPORT=(`grep "$region\b" reported_SNP.txt |cut -f 2`)


echo "                         "
echo "    1. Trans_IndepSNP    "
echo "                         "
if [[ -f LD2/${region}_indepSNP_trans_VS_reportedSNP.ld.txt ]]; then
    rm LD2/${region}_indepSNP_trans_VS_reportedSNP.ld.txt
fi
if [[ -f LD2/${region}_indepSNP_trans_VS_eur.ld.txt ]]; then
    rm LD2/${region}_indepSNP_trans_VS_eur.ld.txt
fi
if [[ -f LD2/${region}_indepSNP_trans_VS_asian.ld.txt ]]; then
    rm LD2/${region}_indepSNP_trans_VS_asian.ld.txt
fi


for a in "${SNP_TRANS[@]}";
do

    echo $a
    signal=`grep "$a\b" $trans |cut -f 2`
     
    echo "                              "
    echo "  Trans VS Report   "
    echo "                              "
    for b in "${SNP_REPORT[@]}";
    do
	echo $b
	plink2 --bfile /path/to/1kg.$chr.phase3.20130502.Eur --ld $a $b | grep "\^" > tmp/$region.out
	# plink files of 1K Genome reference for panel European ancestry

	line=`wc -l tmp/$region.out |awk '{print $1}'`
	if [[ $line = 1 ]]; then
	    r2_eur=`sed -e "s/D'/D/" tmp/$region.out | perl -F"\t" -lane 'if(/r\^2 = (\S+)/){print $1}'`
	    D_eur=`sed -e "s/D'/D/" tmp/$region.out | perl -F"\t" -lane 'if(/D = (\S+)/){print $1}'`
	else
	    r2_eur=`sed -n '2p' tmp/$region.out | sed -e "s/D'/D/" | perl -F"\t" -lane 'if(/r\^2 = (\S+)/){print $1}'`
	    D_eur=`sed -n '2p' tmp/$region.out | sed -e "s/D'/D/" | perl -F"\t" -lane 'if(/D = (\S+)/){print $1}'`
	fi
	
	plink2 --bfile /path/to/1kg.$chr.phase3.20130502.asian --ld $a $b | grep "\^" > tmp/$region.out
	# plink files of 1k Genome reference panel for Asian ancestry

	line=`wc -l tmp/$region.out |awk '{print $1}'`
	if [[ $line = 1 ]]; then
	    r2_asian=`sed -e "s/D'/D/" tmp/$region.out | perl -F"\t" -lane 'if(/r\^2 = (\S+)/){print $1}'`
	    D_asian=`sed -e "s/D'/D/" tmp/$region.out | perl -F"\t" -lane 'if(/D = (\S+)/){print $1}'`
	else
	    r2_asian=`sed -n '2p' tmp/$region.out | sed -e "s/D'/D/" | perl -F"\t" -lane 'if(/r\^2 = (\S+)/){print $1}'`
	    D_asian=`sed -n '2p' tmp/$region.out | sed -e "s/D'/D/" | perl -F"\t" -lane 'if(/D = (\S+)/){print $1}'`
	fi

	echo -e "$region\t$signal\tTrans\t$a\treport\t$b\t$r2_eur\t$D_eur\t$r2_asian\t$D_asian" >> LD2/${region}_indepSNP_trans_VS_reportedSNP.ld.txt
    done

    echo "                             "
    echo "    Trans VS Eur     "
    echo "                             "

    for b in "${SNP_EUR[@]}"; 
    do
	echo $b
	plink2 --bfile /path/to/1kg.$chr.phase3.20130502.Eur --ld $a $b | grep "\^" > tmp/$region.out
	line=`wc -l tmp/$region.out |awk '{print $1}'`
        if [[ $line = 1 ]]; then
            r2_eur=`sed -e "s/D'/D/" tmp/$region.out | perl -F"\t" -lane 'if(/r\^2 = (\S+)/){print $1}'`
            D_eur=`sed -e "s/D'/D/" tmp/$region.out | perl -F"\t" -lane 'if(/D = (\S+)/){print $1}'`
        else
            r2_eur=`sed -n '2p' tmp/$region.out | sed -e "s/D'/D/" | perl -F"\t" -lane 'if(/r\^2 = (\S+)/){print $1}'`
            D_eur=`sed -n '2p' tmp/$region.out | sed -e "s/D'/D/" | perl -F"\t" -lane 'if(/D = (\S+)/){print $1}'`
        fi
	echo -e "$region\t$signal\tTrans\t$a\tEur\t$b\t$r2_eur\t$D_eur" >> LD2/${region}_indepSNP_trans_VS_eur.ld.txt
    done

    echo "                               "
    echo "    Trans VS Asian     "
    echo "                               "  

    for b in "${SNP_ASIAN[@]}";
    do
	echo $b
	plink2 --bfile /path/to/1kg.$chr.phase3.20130502.asian --ld $a $b | grep "\^" > tmp/$region.out
	line=`wc -l tmp/$region.out |awk '{print $1}'`
        if [[ $line = 1 ]]; then
            r2_asian=`sed -e "s/D'/D/" tmp/$region.out | perl -F"\t" -lane 'if(/r\^2 = (\S+)/){print $1}'`
            D_asian=`sed -e "s/D'/D/" tmp/$region.out | perl -F"\t" -lane 'if(/D = (\S+)/){print $1}'`
        else
            r2_asian=`sed -n '2p' tmp/$region.out | sed -e "s/D'/D/" | perl -F"\t" -lane 'if(/r\^2 = (\S+)/){print $1}'`
            D_asian=`sed -n '2p' tmp/$region.out | sed -e "s/D'/D/" | perl -F"\t" -lane 'if(/D = (\S+)/){print $1}'`
        fi
	echo -e "$region\t$signal\tTrans\t$a\tAsian\t$b\t$r2_asian\t$D_asian" >> LD2/${region}_indepSNP_trans_VS_asian.ld.txt
    done

done

if [[ $trans_num > 0 ]]; then
    sed -i '1i region\tsignal\tpopulation\tSNP\tpopulation\tSNP\tr2_Eur\tD_Eur\tr2_Asian\tD_Asian' LD2/${region}_indepSNP_trans_VS_reportedSNP.ld.txt
fi
if [[ $eur_num > 0 ]]; then
    sed -i '1i region\tsignal\tpopulation\tSNP\tpopulation\tSNP\tr2_Eur\tD_Eur' LD2/${region}_indepSNP_trans_VS_eur.ld.txt
fi
if [[ $asian_num > 0 ]]; then
    sed -i '1i region\tsignal\tpopulation\tSNP\tpopulation\tSNP\tr2_Asian\tD_Asian' LD2/${region}_indepSNP_trans_VS_asian.ld.txt
fi



echo "                         "
echo "    2. Eur_IndepSNP    "
echo "                         "
if [[ -f LD2/${region}_indepSNP_eur_VS_reportedSNP.ld.txt ]]; then
    rm LD2/${region}_indepSNP_eur_VS_reportedSNP.ld.txt
fi

for a in "${SNP_EUR[@]}";
do

    echo $a
    signal=`grep "$a\b" $eur |cut -f 2`

    echo "                              "
    echo "   Eur VS Report   "
    echo "                              "
    for b in "${SNP_REPORT[@]}";
    do
	echo $b
	plink2 --bfile /path/to/1kg.$chr.phase3.20130502.Eur --ld $a $b | grep "\^" > tmp/$region.out
	line=`wc -l tmp/$region.out |awk '{print $1}'`
        if [[ $line = 1 ]]; then
            r2_eur=`sed -e "s/D'/D/" tmp/$region.out | perl -F"\t" -lane 'if(/r\^2 = (\S+)/){print $1}'`
            D_eur=`sed -e "s/D'/D/" tmp/$region.out | perl -F"\t" -lane 'if(/D = (\S+)/){print $1}'`
        else
            r2_eur=`sed -n '2p' tmp/$region.out | sed -e "s/D'/D/" | perl -F"\t" -lane 'if(/r\^2 = (\S+)/){print $1}'`
            D_eur=`sed -n '2p' tmp/$region.out | sed -e "s/D'/D/" | perl -F"\t" -lane 'if(/D = (\S+)/){print $1}'`
        fi
	
	echo -e "$region\t$signal\tEur\t$a\treport\t$b\t$r2_eur\t$D_eur" >> LD2/${region}_indepSNP_eur_VS_reportedSNP.ld.txt
    done
done

if [[ $eur_num > 0 ]]; then
    sed -i '1i region\tsignal\tpopulation\tSNP\tpopulation\tSNP\tr2_Eur\tD_Eur' LD2/${region}_indepSNP_eur_VS_reportedSNP.ld.txt
fi

echo "                         "
echo "    3. Asian_IndepSNP    "
echo "                         "

if [[ -f LD2/${region}_indepSNP_asian_VS_reportedSNP.ld.txt ]]; then
    rm LD2/${region}_indepSNP_asian_VS_reportedSNP.ld.txt
fi

for a in "${SNP_ASIAN[@]}";
do
    echo $a
    signal=`grep "$a\b" $asian |cut -f 2`

    echo "                              "
    echo "   Asian VS Report   "
    echo "                              "

    for b in "${SNP_REPORT[@]}";
    do
	echo $b
	plink2 --bfile /path/to/1kg.$chr.phase3.20130502.asian --ld $a $b | grep "\^" > tmp/$region.out
	line=`wc -l tmp/$region.out |awk '{print $1}'`
        if [[ $line = 1 ]]; then
            r2_asian=`sed -e "s/D'/D/" tmp/$region.out | perl -F"\t" -lane 'if(/r\^2 = (\S+)/){print $1}'`
            D_asian=`sed -e "s/D'/D/" tmp/$region.out | perl -F"\t" -lane 'if(/D = (\S+)/){print $1}'`
        else
            r2_asian=`sed -n '2p' tmp/$region.out | sed -e "s/D'/D/" | perl -F"\t" -lane 'if(/r\^2 = (\S+)/){print $1}'`
            D_asian=`sed -n '2p' tmp/$region.out | sed -e "s/D'/D/" | perl -F"\t" -lane 'if(/D = (\S+)/){print $1}'`
        fi
	echo -e "$region\t$signal\tAsian\t$a\treport\t$b\t$r2_asian\t$D_asian" >> LD2/${region}_indepSNP_asian_VS_reportedSNP.ld.txt
    done
done
if [[ $asian_num > 0 ]]; then
    sed -i '1i region\tsignal\tpopulation\tSNP\tpopulation\tSNP\tr2_Asian\tD_Asian' LD2/${region}_indepSNP_asian_VS_reportedSNP.ld.txt
fi
