#!/usr/bin/perl

my %hash;
my %count;
my $snp_file=$ARGV[0];  # example: rs149685934     rs149685934     region_11       Trans-ancestry  1:222148219
my $anno_file=$ARGV[1]; # example: ca_ATACseq_chip-seq_used.path
my $tag;
open F,"$anno_file";
foreach(<F>)
{
    chomp;
    @input=split"\t",$_; # example line: ATACseq 66174  anno/epi_signal/ca/66174_sort_peaks.narrowPeak.bed.hg19
    $tag=$input[0]; # ATACseq
    $bed=$input[1]; # example: anno/epi_signal/ca/66174_sort_peaks.narrowPeak.bed.hg19
    open F,$bed;
    foreach(<F>)
    {
	chomp;
	my @db=split/\t/, $_; # example line: chr1    569130  569786  peak5  
	$chr=$db[0];
	foreach $pos($db[1]..$db[2])
	{
	    $pek="$chr:$pos";
	    $hash{$tag}{$pek}++;
	}                
    }
}
close(F);

my $line;
my $chr;
my $pos;
my $snp;
open F,"$snp_file";
$out=$snp_file.".$tag".".score";
open OUT,">$out";
while(<F>)
{
    chomp;
    my @db=split/\t/;
    if($db[0]=~/ccv/)
    {
	print OUT "$_\t$tag\n";
    }
    
    if($db[4]=~/(\d+):(\d+)/)
    {
	$chr="chr".$1;
	$pos=$2;
	$snp="$chr:$pos";
	if(exists $hash{$tag}{$snp})
	{
		print OUT "$_\t1\n";
	}else{
		print OUT "$_\t0\n";
	}
    }
}
