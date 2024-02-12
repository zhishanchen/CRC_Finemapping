#!/usr/bin/perl
use strict;

my $eqtl="GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_eQTL_all_associations_Colon_Transverse.allpairs.txt";
my $gtf="gencode.v26.annotation.gene.gtf";
my $snp=$ARGV[0];
my $out=$ARGV[1];

my %pos;
open (F,$snp);
while(<F>)
{
	chomp;
	my @db=split/\t/,$_;
	my $chr=$db[0];
	my $start=$db[1];
	my $end=$db[2];
	#my $var=$db[9];
	if($chr=~/chr(\d+)/)
	{
	    my $chrom=$1;
	    $pos{$chrom}{$start}="$_";
	}
}
close F;

my @db;
my $gid;
my %gene;
open (F, $gtf);
while(<F>)
{
	chomp;
	@db=split/\t/;
	if(/^chr/)
	{
	    if($db[8]=~/"(ENSG\d+).\d+"; gene_type/)
		{
		    $gid=$1;
			if($db[8]=~/gene_name "(\S+)";/)
			{
			    $gene{$gid}=$1;
			}
		}	
	}
}
close F;

open OUT,">$out";
print OUT "chr\tpos\tpos\tregionID\tsignalID\tindepSNP\tccv\tpopulation\tag\tvariant\tgeneID\tgene\tmc\tmaf\tpval_nominal\tslope\tslope_se\talt_allele\n";
open F,$eqtl;
my @d;
while(<F>)
{ 
  chomp;
  @d=split/\t/;
  if(/chr(\d+)_(\d+)_\S+_(\S+)_b38/)
  {
     if(exists $pos{$1}{$2})
        {
	    my $chr=$1;
	    my $pos=$2;
	    my $alt=$3;
	    if($d[0]=~/(ENSG\d+).\d+/)
	    {
		my $g=$1;
	    	if(exists $gene{$g})
	    	{
		    print OUT "$pos{$chr}{$pos}\t$g\t$gene{$g}\t$d[4]\t$d[5]\t$d[6]\t$d[7]\t$d[8]\t$alt\n";
		}
	    }
	}	   
  }
}

