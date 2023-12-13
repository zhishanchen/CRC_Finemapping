#!/bin/bash

#
# This script is to get ccv for each indepedent signal.
# author: zhishan chen
# email: zhishan118@gmail.com

#######################################################
#  let me read in the option from command line
#######################################################

export REGION=-1
export CUTOFF=-1
export POP=-1

# now read from command line
while [[ $# -gt 0 ]]
do

    key="$1"

    case $key in
        -r|--region)
            REGION="$2"
            shift # past argument
            shift # past value
            ;;
        -c|--cutoff)
            CUTOFF="$2"
            shift # past argument
            shift # past value
            ;;
	-p|--pop)
	    POP="$2"
	    shift # past argument
	    shift # past value
	    ;;
        *) # unknown option
            echo "The option is invalid, the available options are below:"
            echo "-r(--region): the risk region;"
            echo "-c(--cutoff): the cutoff to define independent signal"
	    echo "-p(--pop): the population for conditional analysis"
	    exit 1
            ;;
    esac
done


#######################################################
####                 Setting                       ####
#######################################################

# directory
WORK_DIR=`pwd|perl -F"/" -lane '$p=join"/",@F[0..($#F-1)];print $p'`
echo "The folder for this project is $WORK_DIR"
GWAS_DIR=`echo $WORK_DIR/gwas/07302022`
echo "The folder with GWAS summaristatistic data is $GWAS_DIR"
COND_DIR=`echo $WORK_DIR/cond.07302022`
echo "The folder with conditional result is $COND_DIR"

# output dir 
if [ ! -d $COND_DIR/$POP/$REGION ]; then
    mkdir $COND_DIR/$POP
    mkdir $COND_DIR/$POP/$REGION
fi
OUT="$COND_DIR/$POP/$REGION"

# reference dir
if [[ $POP == "Trans" || $POP == "Eur" ]]; then
    REF_Eur="/scratch/sbcs/chenzs/New_Public_data/1KG/Eur/bychr"
fi
if [[ $POP == "Trans" || $POP == "Asian" ]]; then
    REF_Asian="/scratch/sbcs/chenzs/New_CRC_GWAS/genotype/Asian/forLD/6684"
fi


#######################################################
####                    GWAS                       ####
#######################################################

if [[ $POP == "Trans" ]]; then
    echo "                                                                  "
    echo "For trancestry conditional analysis, need three gwas summary files"
    echo "                                                                  "
    GWAS_SUM_T="$GWAS_DIR/CRC_consortium_ALL_risk_SNPs.txt_latest.bed.merge.region.gwas_Trans.rsid_MAF0.01_P0.05_gwas"
    GWAS_SUM_A="$GWAS_DIR/CRC_consortium_ALL_risk_SNPs.txt_latest.bed.merge.region.gwas_Trans.rsid_MAF0.01_P0.05_ACCC"
    GWAS_SUM_E="$GWAS_DIR/CRC_consortium_ALL_risk_SNPs.txt_latest.bed.merge.region.gwas_Trans.rsid_MAF0.01_P0.05_Eur"
    
    CHR=`grep "$REGION\b" $GWAS_SUM_T | cut -f 1 | sort |uniq | sed 's/chr//'`
fi

if [[ $POP == "Eur" ]]; then
    GWAS_SUM_E="$GWAS_DIR/CRC_consortium_ALL_risk_SNPs.txt_latest.bed.merge.region.gwas_Eur.rsid_MAF0.01_P1e-4"
    CHR=`grep "$REGION\b" $GWAS_SUM_E | cut -f 1 | sort |uniq |sed 's/chr//'`

fi

if [[ $POP == "Asian" ]]; then
    GWAS_SUM_A="$GWAS_DIR/CRC_consortium_ALL_risk_SNPs.txt_latest.bed.merge.region.gwas_ACCC.rsid_MAF0.01_P1e-4"
    CHR=`grep "$REGION\b" $GWAS_SUM_A | cut -f 1 | sort |uniq | sed 's/chr//'`
fi


#######################################################
####                  FUNCTION                     ####
#######################################################

# conditional analysis
function condition()
{
    if [[ $POP == "Trans" ]]; then
	local esnp=$1
	local asnp=$2
	local step=$3
	local out=$4
    fi
    
    if [[ $POP == "Eur" ]]; then
	local esnp=$1
	local step=$2
	local out=$3
    fi
    
    if [[ $POP == "Asian" ]]; then
	local asnp=$1
	local step=$2
	local out=$3
    fi
    
    if [[ $POP == "Trans" || $POP == "Asian" ]]; then
        /scratch/sbcs/chenzs/CRC_GWAS/software/gcta_1.92.3beta2/gcta64  --bfile $REF_Asian/MEGA_RsqGT03_6684_chr${CHR}.dose.recode  --chr $CHR --maf 0.001 --cojo-file $out/$REGION.a.ma  --cojo-cond $asnp --out $out/$REGION.a.$step
        comp2line.hash.pl  -c 2 -q $out/$REGION.a.$step.cma.cojo -d 1 -db $out/$REGION.a.ma -e | perl -F"\t" -lane '{print "$F[1]\t$F[14]\t$F[15]\t$F[9]\t$F[10]\t$F[11]\t$F[12]\t$F[8]"}' > $out/$REGION.a.$step.cma.cojo_2
    fi

    if [[ $POP == "Trans" || $POP == "Eur" ]]; then
	/scratch/sbcs/chenzs/CRC_GWAS/software/gcta_1.92.3beta2/gcta64 --bfile $REF_Eur/1kg.chr${CHR}.phase3.20130502.Eur --chr ${CHR} --maf 0.001 --cojo-file $out/$REGION.e.ma  --cojo-cond $esnp  --out $out/$REGION.e.$step
       comp2line.hash.pl -c 2 -q $out/$REGION.e.$step.cma.cojo -d 1 -db $out/$REGION.e.ma -e | perl -F"\t" -lane '{print "$F[0]:$F[2]\t$F[14]\t$F[15]\t$F[9]\t$F[10]\t$F[11]\t$F[12]\t$F[8]"}' > $out/$REGION.e.$step.cma.cojo_2
       sed -i 's/Chr:bp/SNP/' $out/$REGION.e.$step.cma.cojo_2    
    fi
}

# metal analysis only be used in trans-ancestry conditional analysis
function meta()
{
    local step=$1
    local out=$2	
    echo "SCHEME   STDERR" > $out/$REGION.metal
    echo "AVERAGEFREQ ON" >> $out/$REGION.metal
    echo "MINMAXFREQ ON" >> $out/$REGION.metal
    echo "############"  >> $out/$REGION.metal
    echo "MARKER SNP" >> $out/$REGION.metal
    echo "ALLELE A1 A2" >> $out/$REGION.metal
    echo "FREQ freq_geno" >> $out/$REGION.metal
    echo "EFFECT bC" >> $out/$REGION.metal
    echo "STDERR bC_se" >> $out/$REGION.metal
    echo "############" >> $out/$REGION.metal
	
    echo "PROCESS $out/$REGION.a.$step.cma.cojo_2" >> $out/$REGION.metal
    echo "PROCESS $out/$REGION.e.$step.cma.cojo_2" >> $out/$REGION.metal
    
    echo "OUTFILE $out/${REGION}_${step}_ .TBL" >> $out/$REGION.metal
    echo "ANALYZE HETEROGENEITY" >> $out/$REGION.metal
    
    /scratch/sbcs/chenzs/CRC_GWAS/software/metal $out/$REGION.metal
}


###################################################################################################
####                                           MAIN ANALYSIS                                   ####
###################################################################################################

echo "                 "
echo "Let's to get ccv for each indepedent signals in risk $REGION"
echo "                 "

if [[ -f $OUT/$REGION.indep.snplist ]] || [[ -f $OUT/$REGION.indep.snplist.eur ]]; then
    
ROUND=1
# input dir
INPUT=$OUT

while [[ $ROUND > 0 ]]
do
    # output dir
    if [ ! -d $COND_DIR/$POP/$REGION/CCV-repeat$ROUND ]; then
	mkdir $COND_DIR/$POP/$REGION/CCV-repeat$ROUND
    fi
    CCV_OUT=$COND_DIR/$POP/$REGION/CCV-repeat$ROUND

    if [[ $POP == "Trans" ]]; then
	echo "            "
	echo "copy gwas summary data to new folder for getting ccv"
	echo "            "
	cp $INPUT/$REGION.e.ma $CCV_OUT
	cp $INPUT/$REGION.a.ma $CCV_OUT
    fi
    if [[ $POP == "Eur"  ]]; then
	cp $INPUT/$REGION.e.ma $CCV_OUT
    fi
    if [[ $POP == "Asian" ]]; then
	cp $INPUT/$REGION.a.ma $CCV_OUT
    fi


    ##########################################################
    ####                          Get CCV                 ####
    ##########################################################

    #########################################
    ####         European and Asian      ####
    #########################################
    
    if [[ $POP == "Eur" || $POP == "Asian" ]]; then
	
    NUM=`wc -l $INPUT/$REGION.indep.snplist | cut -d" " -f 1`

    #################################
    ##  single independent signal  ##
    #################################
    
    if [[ $NUM == "1" ]]; then	
	echo "                                          "
	echo "There is 1 independent signals in $REGION."
	echo "                                          "

	step=1

	####### European ######
	if [[ $POP == "Eur" ]]; then
	    cp $INPUT/$REGION.indep.snplist  $CCV_OUT/snp.list
	    snp="$CCV_OUT/snp.list"
	    indep=`head -n 1 $INPUT/$REGION.indep.snplist`
	    condition $snp $step $CCV_OUT
	    # determine the order of pvalue bases on GWAS
	    PVALUE=`grep "$REGION\b" $GWAS_SUM_E |sort -k16,16g |cut -f 16 |head -n 1`
	    order=`grep "$REGION\b" $GWAS_SUM_E |sort -k16,16g |cut -f 16 |head -n 1 | perl -pe 's/[-\d.]+e(?:\+|(-))?(\d+)/$1$2/gi^V'`
	    order2=`expr $order + 3`
	    # select variants related to analyzed SNP based on conditional analysis
	    sed '1d' $CCV_OUT/$REGION.e.$step.cma.cojo | perl -F"\t" -slane '$t="1e$x";if($F[7] < $t){if($F[12]=~/NA/){print "$y\tsignal_$s\t$z\t$F[1]\t$F[12]\t$F[7]\tEuropean"}else{$fd=($F[7]/$F[12]);if($fd < 0.01){print "$y\tsignal_$s\t$z\t$F[1]\t$F[12]\t$F[7]\tEuropean"}}}' -- -x=$order2  -y=$REGION -z=$indep -s=$step | sort -t$'\t' -k5,5g > $CCV_OUT/$REGION.e.signal_$step.ccv.list.tmp
	    echo -e "$REGION\tsignal_$step\t$indep\t$indep\tNA\t$PVALUE\tEuropean" >> $CCV_OUT/$REGION.e.signal_$step.ccv.list.tmp
	    perl -F"\t" -lane 'if($F[2] eq $F[3]){print "$_\tindep"}else{print "$_\tccv"}' $CCV_OUT/$REGION.e.signal_$step.ccv.list.tmp > $CCV_OUT/$REGION.e.signal_$step.ccv.list
	    sed -i '1i regionID\tsignalID\tindepSNP\tccv\tcond_on_other_indepSNP_pvalue\toriginal_gwas_pvalue\tpopulation\ttag' $CCV_OUT/$REGION.e.signal_$step.ccv.list
	    rm $CCV_OUT/$REGION.e.signal_$step.ccv.list.tmp
	fi

	###### Asian ######
	if [[ $POP == "Asian" ]]; then
	    cp $INPUT/$REGION.indep.snplist  $CCV_OUT/snp.list
	    snp="$CCV_OUT/snp.list"
	    indep=`head -n 1 $INPUT/$REGION.indep.snplist`
	    indep_rs=`grep $indep $GWAS_SUM_A |cut -f 23`
	    condition $snp $step $CCV_OUT
	    # determine the order of pvalue bases on GWAS
	    PVALUE=`grep "$REGION\b" $GWAS_SUM_A |sort -k16,16g |cut -f 16 |head -n 1`
	    order=`grep "$REGION\b" $GWAS_SUM_A |sort -k16,16g |cut -f 16 |head -n 1 | perl -pe 's/[-\d.]+e(?:\+|(-))?(\d+)/$1$2/gi^V'`
	    order2=`expr $order + 3`
	    # select variants related to analyzed SNP based on conditional analysis
	    sed '1d' $CCV_OUT/$REGION.a.$step.cma.cojo | perl -F"\t" -slane '$t="1e$x";if($F[7] < $t){if($F[12]=~/NA/){print "$y\tsignal_$s\t$z\t$F[1]\t$F[12]\t$F[7]\tAsian"}else{$fd=($F[7]/$F[12]);if($fd < 0.01){print "$y\tsignal_$s\t$z\t$F[1]\t$F[12]\t$F[7]\tAsian"}}}' -- -x=$order2 -y=$REGION -z=$indep -s=$step | sort -t$'\t' -k5,5g  > $CCV_OUT/$REGION.a.signal_$step.ccv.list.tmp
	    echo -e "$REGION\tsignal_$step\t$indep\t$indep\tNA\t$PVALUE\tAsian" >> $CCV_OUT/$REGION.a.signal_$step.ccv.list.tmp
	    perl -F"\t" -lane 'if($F[2] eq $F[3]){print "$_\tindep"}else{print "$_\tccv"}' $CCV_OUT/$REGION.a.signal_$step.ccv.list.tmp > $CCV_OUT/$REGION.a.signal_$step.ccv.list
	    sed -i '1i regionID\tsignalID\tindepSNP\tccv\tcond_on_other_indepSNP_pvalue\toriginal_gwas_pvalue\tpopulation\ttag' $CCV_OUT/$REGION.a.signal_$step.ccv.list
	    rm $CCV_OUT/$REGION.a.signal_$step.ccv.list.tmp
	fi
    fi

    ##################################
    ## multiple independent signals ##
    ##################################   
    if [[ $NUM > "1" ]]; then
	echo "                                             "
	echo "There is $NUM independent signals in $REGION."
	echo "                                             "
	for ((i = 1; i <= $NUM; i++)); do
            step=$i
            sed "${step}d" $INPUT/$REGION.indep.snplist > $CCV_OUT/snp.1.list  #other independent signals
	    #sed -n "${step}p" $INPUT/$REGION.indep.snplist > $CCV_OUT/snp.2.list #the analyzed independent signal
            snp1="$CCV_OUT/snp.1.list" 
	    condition $snp1 $step $CCV_OUT  #conditional analysis conditioning on other independent signals

	    ###### European ######
	    if [[ $POP == "Eur" ]]; then
		# get variants with conditional p value within the most significant association after conditioning on other independent associations
		order=`awk '{if($13!~/NA/) print $0}' $CCV_OUT/$REGION.e.$step.cma.cojo |sort -k13,13g |sed -n '2p'|cut -f 13 | perl -pe 's/[-\d.]+e(?:\+|(-))?(\d+)/$1$2/gi^V'`
		order2=`expr $order + 3`
		leadsnp=`awk '{if($13!~/NA/) print $0}' $CCV_OUT/$REGION.e.$step.cma.cojo | sort -k13,13g |sed -n '2p' |cut -f 2 `
		awk '{if($13!~/NA/) print $0}' $CCV_OUT/$REGION.e.$step.cma.cojo |sed '1d' | perl -F"\t" -slane '$t="1e$x";if($F[12] < $t){print "$y\tsignal_$s\t$z\t$F[1]\t$F[12]\t$F[7]\tEuropean"}' -- -x=$order2 -y=$REGION -z=$leadsnp -s=$step > $CCV_OUT/$REGION.e.signal_$step.ccv.list.tmp
		cut -f 4 $CCV_OUT/$REGION.e.signal_$step.ccv.list.tmp > $CCV_OUT/temp

		# get variants after conditioning on the analyzed independent association, pvalue change fold > 100
		echo $leadsnp > $CCV_OUT/snp.2.list
		snp2="$CCV_OUT/snp.2.list"
		/scratch/sbcs/chenzs/CRC_GWAS/software/gcta_1.92.3beta2/gcta64  --bfile $REF_Eur/1kg.chr${CHR}.phase3.20130502.Eur --chr ${CHR} --maf 0.001 --cojo-file $CCV_OUT/$REGION.e.ma --cojo-cond $snp2 --out $CCV_OUT/$REGION.e.$step.temp
		
		sed '1d' $CCV_OUT/$REGION.e.$step.temp.cma.cojo | perl -F"\t" -slane 'if($F[12] !~/NA/){$fd=($F[7]/$F[12]); if($fd < 0.0001){print "$y\tsignal_$s\t$z\t$F[1]\t$F[12]\t$F[7]\tEuropean"}}else{if($F[12] eq "NA"){print "$y\tsignal_$s\t$z\t$F[1]\t$F[12]\t$F[7]\tEuropean"}}' -- -y=$REGION -z=$leadsnp -s=$step |cut -f 4 > $CCV_OUT/temp2
		
		# intersect to get variants which are truely correlated with analyzed varinat and within two orders
		comp2line.hash.pl -c 1 -q $CCV_OUT/temp -d 1 -db $CCV_OUT/temp2 -e |grep -v nocommon|cut -f 1 > $CCV_OUT/temp3

		tag=`wc -l $CCV_OUT/temp2 | awk '{print $1}'` # if no variants after conditioning on the analyzed SNP
		if [[ $tag == "0" ]]; then
		    echo $leadsnp > $CCV_OUT/temp3
		else
		    echo $leadsnp >> $CCV_OUT/temp3
		fi

		comp2line.hash.pl -c 1 -q $CCV_OUT/temp3 -d 4 -db $CCV_OUT/$REGION.e.signal_$step.ccv.list.tmp -e |grep -v nocommon |cut -f 2-8 > $CCV_OUT/$REGION.e.signal_$step.ccv.list.tmp2
		perl -F"\t" -lane  'if($F[2] eq $F[3]){print "$_\tindep"}else{print "$_\tccv"}' $CCV_OUT/$REGION.e.signal_$step.ccv.list.tmp2 > $CCV_OUT/$REGION.e.signal_$step.ccv.list
		sed -i '1i regionID\tsignalID\tindepSNP\tccv\tcond_on_other_indepSNP_pvalue\toriginal_gwas_pvalue\tpopulation\ttag' $CCV_OUT/$REGION.e.signal_$step.ccv.list

		rm $CCV_OUT/$REGION.e.signal_$step.ccv.list.tmp*
		rm $CCV_OUT/temp*		
	    fi

	    ###### for Asian ######
	    if  [[ $POP == "Asian" ]]; then
		# get variants with conditional p value within the most significant association after conditioning on other independent associations
		order=`awk '{if($13!~/NA/) print $0}' $CCV_OUT/$REGION.a.$step.cma.cojo |sort -k13,13g |sed -n '2p'|cut -f 13 | perl -pe 's/[-\d.]+e(?:\+|(-))?(\d+)/$1$2/gi^V'`
		order2=`expr $order + 3`
		leadsnp=`awk '{if($13!~/NA/) print $0}' $CCV_OUT/$REGION.a.$step.cma.cojo | sort -k13,13g |sed -n '2p' |cut -f 2`
		awk '{if($13!~/NA/) print $0}' $CCV_OUT/$REGION.a.$step.cma.cojo |sed '1d' | perl -F"\t" -slane '$t="1e$x";if($F[12] < $t){print "$y\tsignal$s\t$z\t$F[1]\t$F[12]\t$F[7]\tAsian"}' -- -x=$order2 -y=$REGION -z=$leadsnp -s=$step |sort -t$'\t' -k5,5g | perl -F"\t" -lane 'if($F[2] eq $F[3]){print "$_\tindep"}else{print "$_\tccv"}' > $CCV_OUT/$REGION.a.signal_$step.ccv.list.tmp
		cut -f 4 $CCV_OUT/$REGION.a.signal_$step.ccv.list.tmp > $CCV_OUT/temp
		
		# get variants after conditioning on the analyzed independent association, pvalue change fold > 100
		echo $leadsnp > $CCV_OUT/snp.2.list
		snp2="$CCV_OUT/snp.2.list"
		/scratch/sbcs/chenzs/CRC_GWAS/software/gcta_1.92.3beta2/gcta64  --bfile $REF_Asian/MEGA_RsqGT03_6684_chr${CHR}.dose.recode  --chr $CHR --maf 0.001 --cojo-file $CCV_OUT/$REGION.a.ma --cojo-cond $snp2 --out $CCV_OUT/$REGION.a.$step.temp

		sed '1d' $CCV_OUT/$REGION.a.$step.temp.cma.cojo | perl -F"\t" -slane 'if($F[12] !~/NA/){$fd=($F[7]/$F[12]); if($fd < 0.01){print "$y\tsignal_$s\t$z\t$F[1]\t$F[12]\t$F[7]\tAsian"}}else{if($F[12] eq "NA"){print "$y\tsignal_$s\t$z\t$F[1]\t$F[12]\t$F[7]\tAsian"}}' -- -y=$REGION -z=$leadsnp -s=$step |cut -f 4 > $CCV_OUT/temp2
		
		# intersect to get variants which are truely correlated with analyzed varinat and within two orders
		comp2line.hash.pl -c 1 -q $CCV_OUT/temp -d 1 -db $CCV_OUT/temp2 -e |grep -v nocommon|cut -f 1 > $CCV_OUT/temp3

		tag=`wc -l $CCV_OUT/temp2 | awk '{print $1}'` # if no variants after conditioning on the analyzed SNP
		if [[ $tag == "0" ]]; then
		    echo $leadsnp > $CCV_OUT/temp3
		else
		    echo $leadsnp >> $CCV_OUT/temp3
		fi
		
		comp2line.hash.pl -c 1 -q $CCV_OUT/temp3 -d 4 -db $CCV_OUT/$REGION.a.signal_$step.ccv.list.tmp -e |grep -v nocommon |cut -f 2-8 > $CCV_OUT/$REGION.a.signal_$step.ccv.list.tmp2
		perl -F"\t" -lane  'if($F[2] eq $F[3]){print "$_\tindep"}else{print "$_\tccv"}' $CCV_OUT/$REGION.a.signal_$step.ccv.list.tmp2 > $CCV_OUT/$REGION.a.signal_$step.ccv.list
		sed -i '1i regionID\tsignalID\tindepSNP\tccv\tcond_on_other_indepSNP_pvalue\toriginal_gwas_pvalue\tpopulation\ttag' $CCV_OUT/$REGION.a.signal_$step.ccv.list

		rm $CCV_OUT/$REGION.a.signal_$step.ccv.list.tmp*
		rm $CCV_OUT/temp*		
	    fi
	    	    
	done
    fi

         ####### to comfirm the most stable model #######
    
    grep "\bindep\b"  $CCV_OUT/$REGION.*.signal_*.ccv.list | cut -f 3 > $CCV_OUT/$REGION.indep.snplist
    a=`grep -f $INPUT/$REGION.indep.snplist -w $CCV_OUT/$REGION.indep.snplist |wc -l`
    b=`wc -l $INPUT/$REGION.indep.snplist |cut -d" " -f 1`

    if [[ $a -eq $b ]]; then
	echo "the latest indep snp is same with that from forward step"
	echo "don't need more repeat step"
	echo "only $ROUND step of get_ccv"
	ROUND="0"

	## if stop here, change variant ID to rs ID for Asian
	for ((i = 1; i <= $NUM; i++)); do
	    step=$i
	    indep_rs=`comp2line.hash.pl -c 4 -q $CCV_OUT/$REGION.a.signal_$step.ccv.list -d 8 -db $GWAS_SUM_A -e | cut -f 1-8,31 |perl -F"\t" -lane 'if($F[7]=~/indep/){print $F[8]}'`
	    comp2line.hash.pl -c 4 -q $CCV_OUT/$REGION.a.signal_$step.ccv.list -d 8 -db $GWAS_SUM_A -e | sed '1d' | cut -f 1-8,31 |perl -F"\t" -slane 'print "$F[0]\t$F[1]\t$x\t$F[8]\t$F[4]\t$F[5]\t$F[6]\t$F[7]"' -- -x=$indep_rs > $CCV_OUT/temp
	    mv $CCV_OUT/temp $CCV_OUT/$REGION.a.signal_$step.ccv.list
	    sed -i '1i regionID\tsignalID\tindepSNP\tccv\tcond_on_other_indepSNP_pvalue\toriginal_gwas_pvalue\tpopulation\ttag' $CCV_OUT/$REGION.a.signal_$step.ccv.list
	done
	
    else
	echo "the latest indep snp is a little different with that from forward step"
	echo "need additional repeat step"
	INPUT=$CCV_OUT
	ROUND=`expr $ROUND + 1`
    fi
    fi

    #########################################
    ####         Trans-ancestry          ####
    #########################################    

if [[ $POP == "Trans" ]]; then

	NUM=`wc -l $INPUT/$REGION.indep.snplist.asian | cut -d" " -f 1`

	#################################
	##  single independent signal  ##
	#################################
	
	if [[ $NUM == "1" ]]; then

	    step=1
	    
	    cp $INPUT/$REGION.indep.snplist.eur  $CCV_OUT/snp.list.eur
	    cp $INPUT/$REGION.indep.snplist.asian $CCV_OUT/snp.list.asian
	    asnp="$CCV_OUT/snp.list.asian"
	    esnp="$CCV_OUT/snp.list.eur"
	    condition $esnp $asnp $step $CCV_OUT
	    meta $step $CCV_OUT
	    # determine the order of pvalue bases on GWAS
	    indep=`grep "$REGION\b" $GWAS_SUM_T |sort -k16,16g |cut -f 8 | head -n 1`
	    order=`grep "$REGION\b" $GWAS_SUM_T |sort -k16,16g |cut -f 16 |head -n 1 | perl -pe 's/[-\d.]+e(?:\+|(-))?(\d+)/$1$2/gi^V'`
	    order2=`expr $order + 3`
	    grep "$REGION\b" $GWAS_SUM_T |perl -F"\t" -slane '$t="1e$x";if($F[15] < $t){print "$F[4]\tsignal_$s\t$z\t$F[7]\t$F[15]"}' -- -x=$order2 -z=$indep -s=$step  |sort -t$'\t' -k5,5g > $CCV_OUT/$REGION.signal_$step.ccv.list
	    # select variants related to analyzed SNP based on conditional analysis
	    comp2line.hash.pl -c 4 -q $CCV_OUT/$REGION.signal_$step.ccv.list -d 1 -db $CCV_OUT/${REGION}_${step}_1.TBL -e > $CCV_OUT/$REGION.signal_$step.ccv.list.tmp   
	    perl -F"\t" -lane 'if($F[2]eq $F[3]){print "$F[0]\t$F[1]\t$F[2]\t$F[3]\tNA\t$F[4]\tTrans-ancestry\tindep"}else{if($F[5]=~/nocommon/){print "$F[0]\t$F[1]\t$F[2]\t$F[3]\tNA\t$F[4]\tTrans-ancestry\tccv"}else{$fd=($F[4]/$F[14]);if($fd < 0.01){print "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[14]\t$F[4]\tTrans-ancestry\tccv"}}}'  $CCV_OUT/$REGION.signal_$step.ccv.list.tmp > $CCV_OUT/$REGION.signal_$step.ccv.list
	    indep_rs=`comp2line.hash.pl -c 4 -q $CCV_OUT/$REGION.signal_$step.ccv.list -d 8 -db $GWAS_SUM_T -e |perl -F"\t" -lane 'if($F[2] eq $F[3]){print $F[30]}'`
	    comp2line.hash.pl -c 4 -q $CCV_OUT/$REGION.signal_$step.ccv.list  -d 8 -db $GWAS_SUM_T -e |perl -F"\t" -slane 'if($F[2] eq $F[3]){print "$F[0]\t$F[1]\t$x\t$x\t$F[4]\t$F[5]\t$F[6]\t$F[7]"}else{print "$F[0]\t$F[1]\t$x\t$F[30]\t$F[4]\t$F[5]\t$F[6]\t$F[7]"}' -- -x=$indep_rs > $CCV_OUT/$REGION.signal_$step.ccv.list.tmp
       	    sed -i '1i regionID\tsignalID\tindepSNP\tccv\tcond_on_other_indepSNP_pvalue\toriginal_gwas_pvalue\tpopulation\ttag' $CCV_OUT/$REGION.signal_$step.ccv.list.tmp
	    mv $CCV_OUT/$REGION.signal_$step.ccv.list.tmp $CCV_OUT/$REGION.signal_$step.ccv.list
	    
	fi	 


	##################################
	## multiple independent signals ##
	##################################
	
	if [[ $NUM > "1" ]]; then

	    for ((i = 1; i <= $NUM; i++)); do

		step=$i

		sed "${step}d" $INPUT/$REGION.indep.snplist.asian > $CCV_OUT/$REGION.indep.snplist.asian.1
		sed "${step}d" $INPUT/$REGION.indep.snplist.eur > $CCV_OUT/$REGION.indep.snplist.eur.1
		asnp1="$CCV_OUT/$REGION.indep.snplist.asian.1"
		esnp1="$CCV_OUT/$REGION.indep.snplist.eur.1"
		condition $esnp1 $asnp1 $step $CCV_OUT
		meta $step $CCV_OUT

		sed -n "${step}p" $INPUT/$REGION.indep.snplist.asian > $CCV_OUT/$REGION.indep.snplist.asian.2
		sed -n "${step}p" $INPUT/$REGION.indep.snplist.eur > $CCV_OUT/$REGION.indep.snplist.eur.2
		asnp2="$CCV_OUT/$REGION.indep.snplist.asian.2"
		esnp2="$CCV_OUT/$REGION.indep.snplist.eur.2"
		
		# get variants with conditional p value within the most significant association after conditioning on other independent association
		leadsnp=`sort -k10,10g $CCV_OUT/${REGION}_${step}_1.TBL | perl -F"\t" -lane '$F[10]=~s/\?//g; $len=length($F[10]);if($len==2){print $_}'| head -n 1 | cut -f 1`
		order=`sort -k10,10g $CCV_OUT/${REGION}_${step}_1.TBL | perl -F"\t" -lane '$F[10]=~s/\?//g; $len=length($F[10]);if($len==2){print $_}' |head -n 1 |cut -f 10 | perl -pe 's/[-\d.]+e(?:\+|(-))?(\d+)/$1$2/gi'`
		order2=`expr $order + 3`
		sed '1d' $CCV_OUT/${REGION}_${step}_1.TBL | perl -F"\t" -slane '$t="1e$x";if($F[9] < $t){print "$y\tsignal_$s\t$z\t$F[0]\t$F[9]"}' -- -x=$order2 -y=$REGION -z=$leadsnp -s=$step |sort -t$'\t' -k5,5g > $CCV_OUT/$REGION.signal_$step.ccv.list.tmp
		cut -f 4 $CCV_OUT/$REGION.signal_$step.ccv.list.tmp > $CCV_OUT/temp

		# get variants after conditioning on the analyzed independent association, pvalue change fold > 100

		/scratch/sbcs/chenzs/CRC_GWAS/software/gcta_1.92.3beta2/gcta64 --bfile $REF_Eur/1kg.chr${CHR}.phase3.20130502.Eur --chr ${CHR} --maf 0.001 --cojo-file $CCV_OUT/$REGION.e.ma --cojo-cond $esnp2 --out $CCV_OUT/$REGION.e.$step.temp
		comp2line.hash.pl -c 2 -q $CCV_OUT/$REGION.e.$step.temp.cma.cojo -d 1 -db $CCV_OUT/$REGION.e.ma -e | perl -F"\t" -lane '{print "$F[0]:$F[2]\t$F[14]\t$F[15]\t$F[9]\t$F[10]\t$F[11]\t$F[12]\t$F[8]"}' > $CCV_OUT/$REGION.e.$step.temp.cma.cojo_2
		sed -i 's/Chr:bp/SNP/' $CCV_OUT/$REGION.e.$step.temp.cma.cojo_2
		
		/scratch/sbcs/chenzs/CRC_GWAS/software/gcta_1.92.3beta2/gcta64 --bfile $REF_Asian/MEGA_RsqGT03_6684_chr${CHR}.dose.recode --chr $CHR --maf 0.001 --cojo-file $CCV_OUT/$REGION.a.ma --cojo-cond $asnp2 --out $CCV_OUT/$REGION.a.$step.temp
		comp2line.hash.pl  -c 2 -q $CCV_OUT/$REGION.a.$step.temp.cma.cojo -d 1 -db $CCV_OUT/$REGION.a.ma -e | perl -F"\t" -lane '{print "$F[1]\t$F[14]\t$F[15]\t$F[9]\t$F[10]\t$F[11]\t$F[12]\t$F[8]"}' > $CCV_OUT/$REGION.a.$step.temp.cma.cojo_2

		echo "SCHEME   STDERR" > $CCV_OUT/$REGION.metal
		echo "AVERAGEFREQ ON" >> $CCV_OUT/$REGION.metal
		echo "MINMAXFREQ ON" >> $CCV_OUT/$REGION.metal
		echo "############"  >> $CCV_OUT/$REGION.metal
		echo "MARKER SNP" >> $CCV_OUT/$REGION.metal
		echo "ALLELE A1 A2" >> $CCV_OUT/$REGION.metal
		echo "FREQ freq_geno" >> $CCV_OUT/$REGION.metal
		echo "EFFECT bC" >> $CCV_OUT/$REGION.metal
		echo "STDERR bC_se" >> $CCV_OUT/$REGION.metal
		echo "############" >> $CCV_OUT/$REGION.metal

		echo "PROCESS $CCV_OUT/$REGION.a.$step.temp.cma.cojo_2" >> $CCV_OUT/$REGION.metal
		echo "PROCESS $CCV_OUT/$REGION.e.$step.temp.cma.cojo_2" >> $CCV_OUT/$REGION.metal
		
		echo "OUTFILE $CCV_OUT/${REGION}_${step}_temp_ .TBL" >> $CCV_OUT/$REGION.metal
		echo "ANALYZE HETEROGENEITY" >> $CCV_OUT/$REGION.metal

		/scratch/sbcs/chenzs/CRC_GWAS/software/metal $CCV_OUT/$REGION.metal

		comp2line.hash.pl -c 1 -q $CCV_OUT/temp -d 1 -db $CCV_OUT/${REGION}_${step}_temp_1.TBL -e |grep nocommon |cut -f 1 > $CCV_OUT/temp2
		comp2line.hash.pl -c 1 -q $CCV_OUT/temp -d 1 -db $CCV_OUT/${REGION}_${step}_temp_1.TBL -e |grep -v nocommon |cut -f 1 > $CCV_OUT/temp3
		
		comp2line.hash.pl -c 1 -q $CCV_OUT/${REGION}_${step}_temp_1.TBL -d 8 -db $GWAS_SUM_T -e |cut -f 1-15,31 |grep -v nocommon |sed '1d' |perl -F"\t" -lane '$fd=($F[9]/$F[15]);if($fd > 100){print $_}' |cut -f 1 > $CCV_OUT/temp4

		
	      	# intersect to get variants which are truely correlated with analyzed varinat and within two orders
		comp2line.hash.pl -c 1 -q $CCV_OUT/temp3 -d 1 -db $CCV_OUT/temp4 -e |grep -v nocommon|cut -f 1 >  $CCV_OUT/temp5
		cat $CCV_OUT/temp2 $CCV_OUT/temp5 > $CCV_OUT/temp6
		#echo $leadsnp >> $CCV_OUT/temp6
		comp2line.hash.pl -c 1 -q $CCV_OUT/temp6 -d 4 -db $CCV_OUT/$REGION.signal_$step.ccv.list.tmp -e  | grep -v nocommon |cut -f 2-6 > $CCV_OUT/$REGION.signal_$step.ccv.list.tmp2
		comp2line.hash.pl -c 4 -q $CCV_OUT/$REGION.signal_$step.ccv.list.tmp2 -d 8 -db $GWAS_SUM_T -e|cut -f 1-5,21,28 > $CCV_OUT/$REGION.signal_$step.ccv.list.tmp3
		index_rs=`perl -F"\t" -lane 'if($F[2] eq $F[3]){print $F[6]}' $CCV_OUT/$REGION.signal_$step.ccv.list.tmp3`
		perl -F"\t" -slane 'print "$F[0]\t$F[1]\t$x\t$F[6]\t$F[4]\t$F[5]\tTrans-ancestry"' -- -x=$index_rs $CCV_OUT/$REGION.signal_$step.ccv.list.tmp3 > $CCV_OUT/$REGION.signal_$step.ccv.list.tmp4
		perl -F"\t" -lane 'if($F[2] eq $F[3]){print "$_\tindep"}else{print "$_\tccv"}' $CCV_OUT/$REGION.signal_$step.ccv.list.tmp4 > $CCV_OUT/$REGION.signal_$step.ccv.list
		sed -i '1i regionID\tsignalID\tindepSNP\tccv\tcond_on_other_indepSNP_pvalue\toriginal_gwas_pvalue\tpopulation\ttag' $CCV_OUT/$REGION.signal_$step.ccv.list
		
		rm $CCV_OUT/temp*
		rm $CCV_OUT/$REGION.signal_$step.ccv.list.*
		
	    done
	fi

	grep "\bindep\b"  $CCV_OUT/$REGION.signal_*.ccv.list | cut -f 3 > $CCV_OUT/$REGION.indep.snplist.eur	
	comp2line.hash.pl -c 1 -q $CCV_OUT/$REGION.indep.snplist.eur -d 23 -db $GWAS_SUM_T -e |cut -f 9 > $CCV_OUT/$REGION.indep.snplist.asian
	a=`grep -f $INPUT/$REGION.indep.snplist.eur -w $CCV_OUT/$REGION.indep.snplist.eur | wc -l`
	b=`wc -l $INPUT/$REGION.indep.snplist.eur |cut -d" " -f 1`

	if [[ $a -eq $b ]]; then
            echo "the latest indep snp is same with that from forward step"
            echo "don't need more repeat step"
            echo "only $ROUND step of get_ccv"
            ROUND="0"
	else
            echo "the latest indep snp is a little different with that from forward step"
            echo "need additional repeat step"
            INPUT=$CCV_OUT
            ROUND=`expr $ROUND + 1`
	fi
    fi
done
fi 
