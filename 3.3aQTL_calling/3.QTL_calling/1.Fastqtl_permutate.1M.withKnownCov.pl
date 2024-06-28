#! /usr/bin/perl -w
use strict;

my @samp = qw/Ast Exc Inh Mic Oli Opc/;

open OUT,">submit.sh" or die $!;
my $nbatch=22;
my $n=0;

for(my $i=0;;$i++){
	last if !defined $samp[$i];
	my $samp = $samp[$i];
	for(my $chr=1;;$chr++){
		last if $chr>=23;
		next if -e "output.1M/$samp.chr$chr.allCov.out.txt.gz";
		if($n>=$nbatch){
			$n=0;
			print OUT "wait\n";
		}
		$n++;
		open OUT2,">1.$samp.allCov.chr$chr.runFastqtl.sh" or die $!;
		print OUT "nohup sh 1.$samp.allCov.chr$chr.runFastqtl.sh &\n";
		print OUT2 "/data/slurm/xiongxs/software/fastqtl/bin/fastQTL.static --vcf /sugon/nanjh/scAPA_project/01genotype_data/liftover/10WGSfromBroad_match_samples_merged_chrall_qc_hg38_dedup.vcf.gz --bed ../${samp}_imputation_MI_3_naomit.parse.voom.sort.bed.gz --region $chr --out output.1M/$samp.chr$chr.allCov.out.txt.gz --cov ../4.PEER.cor/$samp.PEER_covariates.merged.txt --window 1e6 --permute 1000\n";
		close OUT2;
	}
}




