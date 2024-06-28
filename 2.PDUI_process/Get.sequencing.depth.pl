#! /usr/bin/perl -w
use strict;

my $ct = $ARGV[0];
open IN,"./1.cellranger.mapping/samples.txt" or die $!;
open OUT,">mapping_wig_location_with_depth.$ct.txt" or die $!;
#print OUT "wig\tdepth\n";

while(<IN>){
	chomp;
	my $samp=$_;
	my $depth = `samtools view -F 4 ../1.cellranger.mapping/$samp/outs/$samp.$ct.possorted_genome_bam.bam | wc -l`;
	chomp $depth;
	print OUT "../1.cellranger.mapping/$samp/outs/$samp.$ct.possorted_genome_bam.wig\t$depth\n";
}



