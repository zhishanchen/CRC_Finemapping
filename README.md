## CRC_Finemapping paper
**Title**  	       
Fine-mapping analysis including over 254,000 East Asian and European descendants identifies 136 putative colorectal cancer susceptibility genes

**Abstract**   
Genome-wide association studies (GWAS) have identified more than 200 common genetic variants independently associated with colorectal cancer (CRC) risk, but the causal variants and target genes are mostly unknown. We sought to fine-map all known CRC risk loci using GWAS data from 100,204 cases and 154,587 controls of East Asian and European ancestry. Our stepwise conditional analyses revealed 238 independent association signals of CRC risk, each with a set of credible causal variants (CCVs), of which 28 signals had a single CCV. Our cis-eQTL/mQTL and colocalization analyses using colorectal tissue-specific transcriptome and methylome data separately from 1299 and 321 individuals, along with functional genomic investigation, uncovered 136 putative CRC susceptibility genes, including 56 genes not previously reported. Analyses of single-cell RNA-seq data from colorectal tissues revealed 17 putative CRC susceptibility genes with distinct expression patterns in specific cell types. Analyses of whole exome sequencing data provided additional supports for several target genes identified in this study as CRC susceptibility genes. Enrichment analyses of the 136 genes uncover pathways not previously linked to CRC risk. Our study substantially expanded association signals for CRC and provided additional insight into the biological mechanisms underlying CRC development.

## Main analysis 	
### 1. Identify independent risk signal and CCVs
```
bash GWAS.cond.sh -r region_111 -c  0.000001 -p Eur
bash GWAS.getccv.sh -r region_111 -c 0.000001 -p Eur
```
### 2.  LD estimation between variants
```
cat region.list |parallel -q echo bash ld.plink.sh {} |bash
``` 
### 3. functional genomic data analysis
```
bash distal.sh
bash promoter.sh
bash coding.sh
``` 

### 4. eQTL analysis 
#### trans-ancestry
```
bash Trans_eQTL.ACCC.sh
Rscript Trans_eQTL.Coloc.R
bash Trans_eQTL_GTEx.sh
Rscript meta_Trans.indep.step1.R
perl -F"\t" -lane 'if(($F[14] ne $F[10]) && ($F[14] ne "NA") && ($F[10] ne "NA")){$f=join"\t",@F[0..6];$e=join"\t",@F[8..$#F];$F[7]=($F[7]*-1);print "$f\t$F[7]\t$e"}else{print $_}' Trans_ALL_eQTLs.txt > Trans_ALL_eQTLs.txt.v2
Rscript meta_Trans.indep.step2.R
```

#### European
```
bash Eur_eQTL.GTEx.sh
Rscript Eur_eQTL.Coloc.R
Rscript meta_Eur.step1.R
Rscript meta_Eur.step2.R
```    
#### East Asian
```
bash Asian_eQTL.sh
Rscript Asian_eQTL.R
```

