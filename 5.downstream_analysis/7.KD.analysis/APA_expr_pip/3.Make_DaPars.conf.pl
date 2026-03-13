#! /usr/bin/perl -w
use strict; use warnings;

#my @cts = qw/Ast Exc Inh Mic Oli Opc Vas/;

#foreach my $ct(@cts){
my ($arg1, $arg2) = @ARGV;
chdir $arg2 or die $!;
open OUT,">Dapars2_configure_file" or die $!;

my $wigfiles;

open IN, $arg1 or die $!;
while(<IN>){
	chomp;
	my $samp=$_;
	my $wig = "04.bigwig/$samp.wig";
	if(!defined $wigfiles){$wigfiles=$wig;next;}
	$wigfiles = "$wigfiles,$wig";
}

print OUT "# Specify the reference of 3'UTR region

Annotated_3UTR=/data/slurm/xiongxs/Proj_APA/2.DaPars2/Ref/GencodeV32.20230414.hg38_3UTR_annotation.uniq.bed
# A comma separated list of wig files of all samples

Aligned_Wig_files=$wigfiles

Output_directory=Dapars2

Output_result_file=Dapars2

# Specify Coverage threshold

Coverage_threshold=80

# Specify the number of threads to process the analysis

Num_Threads=50

# Provide sequencing depth file for normalization

sequencing_depth_file=mapping_wig_location_with_depth.txt\n";


close IN;
close OUT;
#}

