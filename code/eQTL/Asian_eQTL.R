
library("methods")
library("stringr")
library("data.table")
library("dplyr")

args<-commandArgs(TRUE)
chrom <- args[1]
outdir <- args[2]
ccv_file <- args[3]

"%&%" <- function(a,b) paste(a,b, sep = "")

### inputs ###
# 0/1/2
genotype_file <- paste0("/path/to/chr",chrom,".gt.gt")
# snp annotation
snp_annot_file <- paste0("/path/to/annotation/Asian_MODELsnp.chr",chrom,".annotation")
#

# expression
expression <- "/path/to/381samples_used_364sample_QN.expression_PEER_inverse.txt"

# gene annotation
gene_annot_file <- "gencode.v19.annotation.gtf"


### functions ###

get_gene_annotation <- function(gene_annot_file_name, chrom)
{
    gene_df <- read.table(gene_annot_file,header=F,stringsAsFactors =F,sep="\t",fill=T)
    gene_df1 <- filter(gene_df,V3 %in% "gene")
    geneid <- str_extract(gene_df1[,9], "ENSG\\d+.\\d+")
    genename <- gsub("gene_name (\\S+);","\\1",str_extract(gene_df1[,9], "gene_name (\\S+);"), perl=T)
    gene_used <- as.data.frame(cbind(geneid,genename,gene_df1[,c(1,4,5,3)]))
    colnames(gene_used) <- c("geneid","genename","chr","start","end","anno")
    gtf_used <- filter(gene_used,gene_used[,3] %in% ('chr' %&% chrom))
    gtf_used
}


get_gene_expression <- function(gene_expression_file_name, gene_annot) {
  expr_df <- as.data.frame(read.table(gene_expression_file_name, header = T, stringsAsFactors = F, row.names = 1))
  expr_df <- expr_df %>% select(one_of(intersect(gene_annot$geneid, colnames(expr_df))))
  row.names(expr_df) <- gsub("HCES2.(\\d+)","HCES2-\\1",row.names(expr_df),perl=T)
  expr_df <- expr_df[order(row.names(expr_df)), ]
  expr_df
}

get_flank_gene <- function(coord, gene_annot, cis_window) {
  gene_info <- (gene_annot %>% filter((start >= (coord - cis_window)  & (end <= (coord + cis_window)))))
  flankgenes <- as.character(gene_info[,1])
  flankgenes
}


############ association analysis  ######
cis_window=1000000

# genes in chrom
gene_annot <- get_gene_annotation(gene_annot_file, chrom)
expr_df <- get_gene_expression(expression, gene_annot)
dim(expr_df)
samples <- rownames(expr_df)
n_samples <- length(samples)
genes <- colnames(expr_df)
n_genes <- length(expr_df)

# genotype
#ccv
ccv <- read.table(ccv_file,header=F,stringsAsFactors=F,fill=T,sep="\t")
# snp annotation
snp_annot <- read.table(snp_annot_file, header = T, stringsAsFactors = F)
snp_annot$varID2 <- gsub("\\d+:\\d+;(\\d+:\\d+):\\S+:\\S+","\\1",snp_annot$varID,perl=T)
snp_annot_used <- snp_annot[snp_annot$varID2 %in% ccv[,19],] # snp annotation for ccv in chrom

# head of vcf file
vcfhead <- read.table("/path/to/vcf.head",header = F, stringsAsFactors = F)
vcfhead <- gsub("WG\\d+-\\S+-(CNUHHCRC_\\d+)@\\d+","\\1",vcfhead, perl=TRUE)
vcfhead <- gsub('CNUHH-CRC-(\\d+)_QC-3158-Cai_\\S+','CNUHHCRC_\\1',vcfhead, perl=T)
vcfhead <- gsub('FMMU-(CC\\d+N)_QC-3158-Cai_\\S+','FMMU_\\1',vcfhead, perl =T)
vcfhead <- gsub('(HCES2-\\d+)_QC-3158-Cai_\\S+','\\1', vcfhead, perl=T)

## for 0/1/2
gt <- fread(genotype_file,header=F,stringsAsFactors = F)
colnames(gt) <- vcfhead
gt <- as.data.frame(gt)
row.names(gt) <- gt$ID
gt <- gt [,-c(1:5)]
gt1  <- as.data.frame(t(gt))
gt_used <- gt1[,colnames(gt1) %in% snp_annot_used$varID]
gt_used <- gt_used[rownames(gt_used) %in% samples,]
gt_used <- gt_used[order(rownames(gt_used)),]

# run association

# output
fp_out<-file(paste0(outdir,"/ACCC_eQTL_chr",chrom,"_output.txt"),open="w")

sink(fp_out)

cat("SNP","gene","gene_name","mc","maf","beta","se","t-pvalue","p-value","alt_allele","\n",sep="\t")

for (i in 1:dim(gt_used)[2])
{
    snp <- colnames(gt_used)[i]
    snp2 <- gsub("\\d+:\\d+;(\\d+:\\d+:\\S+:\\S+)","\\1",snp,perl=T)
    genotype <- gt_used[,i]
    rsid <- (snp_annot_used %>% filter(snp_annot_used$varID == snp))$SNP
    coord <- (snp_annot_used %>% filter(snp_annot_used$varID == snp))$pos
    flankgenes <- get_flank_gene(coord,gene_annot,cis_window)
    flankgenes_used <- intersect(flankgenes,colnames(expr_df)) # flankgenes with expression

    for(j in 1:length(flankgenes_used))
    {

	 gene <- flankgenes_used[j]
         gene_name <- as.character(gene_annot$genename[gene_annot$geneid == gene])
         expression <- expr_df[,gene]

	 mc=sum(genotype)
    	 af=mean(genotype) / 2
  	 alt <- (snp_annot %>% filter(varID == snp))$effect
	 x1<-summary(lm(expression ~ genotype))$coefficients[2,c(1,2,3,4)]
    	  cat(rsid,snp2,gene,gene_name,mc,af,x1,alt,"\n",sep="\t")
    }
}
sink()

