#! /usr/bin/perl -w
use strict;

my @cts = qw/Ast Exc Inh Mic Oli Opc Vas/;

foreach my $ct(@cts){

open OUT,">Dapars2_configure_file.$ct" or die $!;

my $wigfiles;

open IN,"./1.cellranger.mapping/samples.txt" or die $!;
while(<IN>){
	chomp;
	my $samp=$_;
	my $wig = "wig/$samp.$ct.possorted_genome_bam.wig";
	if(!defined $wigfiles){$wigfiles=$wig;next;}
	$wigfiles = "$wigfiles,$wig";
}

print OUT "# Specify the reference of 3'UTR region

Annotated_3UTR = ./Ref/GencodeV32.20230414.hg38_3UTR_annotation.uniq.bed

# A comma separated list of wig files of all samples

Aligned_Wig_files=$wigfiles

Output_directory=Dapars2_$ct

Output_result_file=Dapars2_$ct

# Specify Coverage threshold

Coverage_threshold=5

# Specify the number of threads to process the analysis

Num_Threads=12

# Provide sequencing depth file for normalization

sequencing_depth_file=mapping_wig_location_with_depth.$ct.txt\n";

close IN;
close OUT;
}

