#! /usr/bin/perl -w
use strict;

my @samp = qw/Ast Exc Inh Mic Oli Opc/;
my @start = qw/1 2/;

open OUT,">submit.sh" or die $!;

my $nbatch=20;  
my $n=0;

foreach my $samp(@samp){
 foreach my $start(@start){
	 for(my $num=1;;$num++){
		 last if $num>20;
		 next if $start > $num;
	next if -e "output.1M/$samp.chr22.pf$start-$num.out.txt.gz";
	if($n>=$nbatch){
		$n=0;
		print OUT "wait\n";
	}
	$n++;
	open OUT2,">1.$samp.pf$start-$num.runFastqtl.sh" or die $!;
	print OUT "nohup sh 1.$samp.pf$start-$num.runFastqtl.sh &\n";
	for(my $i=1;;$i++){
		last if $i==23;
		next if -e "output.1M/$samp.chr$i.pf$start-$num.out.txt.gz";
		print OUT2 "/data/slurm/xiongxs/software/fastqtl/bin/fastQTL.static --vcf /sugon/nanjh/scAPA_project/01genotype_data/liftover/10WGSfromBroad_match_samples_merged_chrall_qc_hg38_dedup.vcf.gz --bed ../${samp}_imputation_MI_3_naomit.parse.voom.sort.bed.gz --region $i --out output.1M/$samp.chr$i.pf$start-$num.out.txt.gz --cov ../4.PEER.cor/$samp.PEER_covariates.pf$start-$num.txt --window 1e6 --threshold 0.0001\n";
	}
	close OUT2;
		}
	}
}




