#!/bin/bash

#
# This script is to perform forward stepwise conditional analysis and get ccv for each indepedent signal.
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

# work directory
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

# LD reference dir
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

    tag=`grep "$REGION\b" $GWAS_SUM_T |wc -l`
    if [[ $tag == "0" ]]; then
	echo "No variants in $REGION in $POP" > $OUT/$REGION.cond.log
	exit
    fi

fi

if [[ $POP == "Eur" ]]; then
    GWAS_SUM_E="$GWAS_DIR/CRC_consortium_ALL_risk_SNPs.txt_latest.bed.merge.region.gwas_Eur.rsid_MAF0.01_P1e-4"
    CHR=`grep "$REGION\b" $GWAS_SUM_E | cut -f 1 | sort |uniq |sed 's/chr//'`

    tag=`grep "$REGION\b" $GWAS_SUM_E |wc -l`
    if [[ $tag == "0" ]]; then
	echo "No variants in $REGION in $POP" > $OUT/$REGION.cond.log
	exit
    fi
fi

if [[ $POP == "Asian" ]]; then
    GWAS_SUM_A="$GWAS_DIR/CRC_consortium_ALL_risk_SNPs.txt_latest.bed.merge.region.gwas_ACCC.rsid_MAF0.01_P1e-4"
    CHR=`grep "$REGION\b" $GWAS_SUM_A | cut -f 1 | sort |uniq | sed 's/chr//'`
    
    tag=`grep "$REGION\b" $GWAS_SUM_A |wc -l`
    if [[ $tag == "0" ]]; then
	echo "No variants in $REGION in $POP" > $OUT/$REGION.cond.log
	exit
    fi
fi


#######################################################
####                  FUNCTION                     ####
#######################################################

## conditional analysis ##
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
	/nobackup/sbcs/scratch/chenzs/CRC_GWAS/software/gcta_1.92.3beta2/gcta64 --bfile $REF_Asian/MEGA_RsqGT03_6684_chr${CHR}.dose.recode  --chr $CHR --maf 0.001 --cojo-file $out/$REGION.a.ma  --cojo-cond $asnp --out $out/$REGION.a.$step
	perl /data1/chenz27/tmp/myscript/script/comp2line.hash.pl  -c 2 -q $out/$REGION.a.$step.cma.cojo -d 1 -db $out/$REGION.a.ma -e | perl -F"\t" -lane '{print "$F[1]\t$F[14]\t$F[15]\t$F[9]\t$F[10]\t$F[11]\t$F[12]\t$F[8]"}' > $out/$REGION.a.$step.cma.cojo_2
    fi

    if [[ $POP == "Trans" || $POP == "Eur" ]]; then
	/nobackup/sbcs/scratch/chenzs/CRC_GWAS/software/gcta_1.92.3beta2/gcta64 --bfile $REF_Eur/1kg.chr${CHR}.phase3.20130502.Eur --chr ${CHR} --maf 0.001 --cojo-file $out/$REGION.e.ma  --cojo-cond $esnp  --out $out/$REGION.e.$step
	perl /data1/chenz27/tmp/myscript/script/comp2line.hash.pl -c 2 -q $out/$REGION.e.$step.cma.cojo -d 1 -db $out/$REGION.e.ma -e | perl -F"\t" -lane '{print "$F[0]:$F[2]\t$F[14]\t$F[15]\t$F[9]\t$F[10]\t$F[11]\t$F[12]\t$F[8]"}' > $out/$REGION.e.$step.cma.cojo_2
	sed -i 's/Chr:bp/SNP/' $out/$REGION.e.$step.cma.cojo_2    
    fi
}

## metal analysis ##
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
    
    /nobackup/sbcs/scratch/chenzs/CRC_GWAS/software/metal $out/$REGION.metal
}


#######################################################
####               MAIN ANALYSIS                   ####
#######################################################

####################################################
# prepare gwas summary data for conditional analysis
####################################################
echo "                                                  "
echo "extract GWAS summary data for SNPs in each region "
echo "                                                  "

if [[ $POP == "Trans" || $POP == "Asian" ]]; then
    echo "This step need strictly formatted GWAS summary data"
    ROW_marker="17"
    ROW_START=`expr $ROW_marker - 7`
    ROW_END=`expr $ROW_marker - 2`
    ROW_SNP=`expr $ROW_marker - 10`
    ROW_SAMPLE=`expr $ROW_marker + 2`
    grep "$REGION\b" $GWAS_SUM_A | perl -F"\t" -slane '$gwas=join"\t",@F[$x..$y];print "$F[$z]\t$gwas\t$F[$s]"' -- -x=$ROW_START -y=$ROW_END -z=$ROW_SNP -s=$ROW_SAMPLE > $OUT/$REGION.a.ma
    sed -i '1i SNP\tA1\tA2\tfreq\tb\tse\tp\tN' $OUT/$REGION.a.ma

fi


if [[ $POP == "Trans" || $POP == "Eur" ]]; then
    echo "This step need strictly formatted GWAS summary data"
    ROW_marker="17"
    ROW_START=`expr $ROW_marker - 7`
    ROW_END=`expr $ROW_marker - 2`
    ROW_SNP=`expr $ROW_marker + 5`
    ROW_SAMPLE=`expr $ROW_marker + 2`
    grep "$REGION\b" $GWAS_SUM_E | perl -F"\t" -slane '$gwas=join"\t",@F[$x..$y];print "$F[$z]\t$gwas\t$F[$s]"' -- -x=$ROW_START -y=$ROW_END -z=$ROW_SNP -s=$ROW_SAMPLE > $OUT/$REGION.e.ma

    if [[ $REGION == "region_31" ]]; then	
	grep -vE 'rs4421005' $OUT/$REGION.e.ma > $OUT/temp # For region_31 in European, the top1 SNP are not in 1KG, use the second one
	mv $OUT/temp $OUT/$REGION.e.ma
    fi

    if [[ $REGION == "region_11" ]]; then
	sed -i 's/rs12024666/rs149685934/g' $OUT/$REGION.e.ma # For region_11 in European, rs12024666 replaced by rs149685934 which is in 1KG ref
    fi
    
    sed -i '1i SNP\tA1\tA2\tfreq\tb\tse\tp\tN' $OUT/$REGION.e.ma
    
fi

#####################################################
# get the most significant association in each region
#####################################################
echo "                  "
echo " get the lead SNP "
echo "                                                                    "
echo " only focus on region with the most significant association < 1e-06 "
echo "                                                                    "

STEP=1

if [[ $POP == "Trans" ]]; then

    if [[ $REGION == "region_11" ]]; then
	sed -i 's/rs12024666/rs149685934/g' $GWAS_SUM_T # For region_11 in European, rs12024666 replaced by rs149685934 which is in 1KG ref
    fi
    
    ROW_SNP_a=`expr $ROW_marker - 9`
    ROW_SNP_e=`expr $ROW_marker + 6`
    ROW_SNP_PVALUE=`expr $ROW_marker - 1`
    leadsnp_a=`grep "$REGION\b" $GWAS_SUM_T | head -n 1 |cut -f $(echo $ROW_SNP_a)`
    leadsnp_e=`grep "$REGION\b" $GWAS_SUM_T | head -n 1 |cut -f $(echo $ROW_SNP_e)`
    PVALUE=`grep "$REGION\b" $GWAS_SUM_T | head -n 1 |cut -f $(echo $ROW_SNP_PVALUE)`
    echo "The most significnat association is $leadsnp_e with pvalue = $PVALUE in $POP in $REGION" > $OUT/$REGION.cond.log
    r=$(awk 'BEGIN{print ('$PVALUE'>'0.000001')?1:0}')
    if [[ $r == 1 ]]; then
	echo "The p value of most significant association is larger than 1e-06 in $POP in $REGION" >> $OUT/$REGION.cond.log
	exit
    fi
fi

if [[ $POP == "Eur" ]]; then

    if [[ $REGION == "region_11" ]]; then
	sed -i 's/rs12024666/rs149685934/g' $GWAS_SUM_E # For region_11 in European, rs12024666 replaced by rs149685934 which is in 1KG ref
    fi	
    
    ROW_SNP=`expr $ROW_marker + 6` 
    leadsnp=`grep "$REGION\b" $GWAS_SUM_E  | head -n 1 |cut -f $(echo $ROW_SNP)`
    
    if [[ $REGION == "region_31" ]]; then
	leadsnp=`grep "$REGION\b" $GWAS_SUM_E | grep -vE 'rs4421005' | head -n 1 | cut -f $(echo $ROW_SNP)` # For region_31 in European, the top1 SNP are not in 1KG, use the second one
    fi
    
    ROW_SNP_PVALUE=`expr $ROW_marker - 1`
    PVALUE=`grep "$REGION\b" $GWAS_SUM_E | head -n 1 |cut -f $(echo $ROW_SNP_PVALUE)`
    echo "The most significnat association is $leadsnp with pvalue = $PVALUE in $POP in $REGION" > $OUT/$REGION.cond.log

    if [[ $REGION == "region_16" || $REGION == "region_21" || $REGION == "region_95" || $REGION == "region_121" ]]; then
	r=$(awk 'BEGIN{print ('$PVALUE'>'0.0001')?1:0}') ## for 6 regions with most significant association with 1E-06 < P < 1E-04
	if [[ $r == 1 ]]; then
	    echo "The p value of most significant association is larger than 1e-04 in $POP in $REGION" >> $OUT/$REGION.cond.log
	    exit
	fi
    else
	r=$(awk 'BEGIN{print ('$PVALUE'>'0.000001')?1:0}')
	if [[ $r == 1 ]]; then
	    echo "The p value of most significant association is larger than 1e-06 in $POP in $REGION" >> $OUT/$REGION.cond.log
	    exit
	fi
    fi

fi

if [[ $POP == "Asian" ]]; then
    ROW_SNP=`expr $ROW_marker - 9`
    leadsnp=`grep "$REGION\b" $GWAS_SUM_A | head -n 1 |cut -f $(echo $ROW_SNP)`
    leadsnp_rsid=`grep "$REGION\b" $GWAS_SUM_A | head -n 1 |cut -f 23`
    ROW_SNP_PVALUE=`expr $ROW_marker - 1`
    PVALUE=`grep "$REGION\b" $GWAS_SUM_A | head -n 1 |cut -f $(echo $ROW_SNP_PVALUE)`
    echo "The most significnat association is $leadsnp_rsid with pvalue = $PVALUE in $POP in $REGION" > $OUT/$REGION.cond.log

    if [[ $REGION == "region_86" || $REGION == "region_106" ]]; then
        r=$(awk 'BEGIN{print ('$PVALUE'>'0.0001')?1:0}') ## for 6 regions with most significant association with 1E-06 < P < 1E-04
        if [[ $r == 1 ]]; then
            echo "The p value of most significant association is larger than 1e-04 in $POP in $REGION" >> $OUT/$REGION.cond.log
            exit
        fi
    else
        r=$(awk 'BEGIN{print ('$PVALUE'>'0.000001')?1:0}')
        if [[ $r == 1 ]]; then
            echo "The p value of most significant association is larger than 1e-06 in $POP in $REGION" >> $OUT/$REGION.cond.log
            exit
        fi
    fi
    
fi
echo "                                                                    "

##############################################
#  conditional analysis
##############################################

# Trans-ancestry

if [[ $POP == "Trans" ]]; then

    echo $leadsnp_a > $OUT/$REGION.indep.snplist.asian
    echo $leadsnp_e > $OUT/$REGION.indep.snplist.eur

    asnp="$OUT/$REGION.indep.snplist.asian"
    esnp="$OUT/$REGION.indep.snplist.eur"
    
    condition $esnp $asnp $STEP $OUT
    meta $STEP $OUT

    PVALUE=`sort -k10,10g $OUT/${REGION}_${STEP}_1.TBL | head -n 2 | sed '1d' | perl -F"\t" -lane 'if($F[10] eq "++" || $F[10] eq "--"){print $F[9]}else{print "1"}'`
    pva=`printf "%.7f" $PVALUE`
    
    tbl=`grep -v "MarkerName"  $OUT/${REGION}_${STEP}_1.TBL |wc -l`
    if [[ $tbl == "0" ]]; then
	echo "no variant is retained after round $STEP conditional analysis"
	pva="1"
    fi 

    while [[ $pva < $CUTOFF ]]
    do
        leadsnp_a=`sort -k10,10g $OUT/${REGION}_${STEP}_1.TBL | perl -F"\t" -lane 'if($F[10] eq "++" || $F[10] eq "--"){print $_}' | head -n 1 | cut -f 1`
	leadsnp_e=`grep $leadsnp_a $GWAS_SUM_T |cut -f 23`
	echo $leadsnp_a >> $OUT/$REGION.indep.snplist.asian
	echo $leadsnp_e >> $OUT/$REGION.indep.snplist.eur      
	asnp="$OUT/$REGION.indep.snplist.asian"
	esnp="$OUT/$REGION.indep.snplist.eur"
	STEP=`expr $STEP + 1`
    
	condition $esnp $asnp $STEP $OUT
	meta $STEP $OUT
    
	PVALUE=`sort -k10,10g $OUT/${REGION}_${STEP}_1.TBL | perl -F"\t" -lane 'if($F[10] eq "++" || $F[10] eq "--"){print $_}' | head -n 1 | cut -f 10`
	pva=`printf "%.7f" $PVALUE`
	tbl=`grep -v "MarkerName"  $OUT/${REGION}_${STEP}_1.TBL |wc -l`
	if [[ $tbl == "0" ]]; then
            echo "no variant is retained after round $STEP conditional analysis"
            pva="1"
	fi
    done

    echo "stepwise conditional analysis completed"
    echo "there is (are) $STEP independent signal(s) in $POP in $REGION" >> $OUT/$REGION.cond.log
    
fi

# European and Asian

if [[ $POP == "Eur" || $POP == "Asian" ]]; then

    echo $leadsnp > $OUT/$REGION.indep.snplist
    snp="$OUT/$REGION.indep.snplist"

    condition $snp $STEP $OUT

    if [[ $POP == "Eur" ]]; then
	PVALUE=`awk '{if($13!~/NA/) print $0}' $OUT/$REGION.e.$STEP.cma.cojo | sort -k13,13g |sed -n '2p' |cut -f 13`
	pva=`printf "%.7f" $PVALUE`
	cma=`awk '{if($13!~/NA/) print $0}' $OUT/$REGION.e.$STEP.cma.cojo  |grep -v "pC" |wc -l`
    fi
    
    if [[ $POP == "Asian" ]]; then
	PVALUE=`awk '{if($13!~/NA/) print $0}' $OUT/$REGION.a.$STEP.cma.cojo | sort -k13,13g |sed -n '2p' |cut -f 13`
	pva=`printf "%.7f" $PVALUE`
	cma=`awk '{if($13!~/NA/) print $0}' $OUT/$REGION.a.$STEP.cma.cojo  |grep -v "pC" |wc -l`
    fi

    if [[ $cma == "0" ]]; then
       echo "no variant is retained after round $STEP conditional analysis"
       pva="1"
    fi

    while [[ $pva < $CUTOFF ]]
    do

	if [[ $POP == "Eur" ]]; then
	    leadsnp=`awk '{if($13!~/NA/) print $0}' $OUT/$REGION.e.$STEP.cma.cojo | sort -k13,13g |sed -n '2p' |cut -f 2`
	fi
	if [[ $POP == "Asian" ]]; then
	    leadsnp=`awk '{if($13!~/NA/) print $0}' $OUT/$REGION.a.$STEP.cma.cojo | sort -k13,13g |sed -n '2p' |cut -f 2`
	fi

	echo $leadsnp >> $OUT/$REGION.indep.snplist
	snp="$OUT/$REGION.indep.snplist"
	STEP=`expr $STEP + 1`
	condition $snp $STEP $OUT
	
	if [[ $POP == "Eur" ]]; then
	    PVALUE=`awk '{if($13!~/NA/) print $0}' $OUT/$REGION.e.$STEP.cma.cojo | sort -k13,13g |sed -n '2p' |cut -f 13`
	    pva=`printf "%.7f" $PVALUE`
	    cma=`awk '{if($13!~/NA/) print $0}' $OUT/$REGION.e.$STEP.cma.cojo  |grep -v "pC" |wc -l`
	fi
	if [[ $POP == "Asian" ]]; then
	    PVALUE=`awk '{if($13!~/NA/) print $0}' $OUT/$REGION.a.$STEP.cma.cojo | sort -k13,13g |sed -n '2p' |cut -f 13`
	    pva=`printf "%.7f" $PVALUE`
	    cma=`awk '{if($13!~/NA/) print $0}' $OUT/$REGION.a.$STEP.cma.cojo  |grep -v "pC" |wc -l`
	fi

	if [[ $cma == "0" ]]; then
            echo "no variant is retained after round $STEP conditional analysis"
            pva="1"
	fi
	
    done

    echo "stepwise conditional analysis completed"
    echo "there is (are) $STEP independent signal(s) in in $POP in $REGION" >> $OUT/$REGION.cond.log
    
fi
