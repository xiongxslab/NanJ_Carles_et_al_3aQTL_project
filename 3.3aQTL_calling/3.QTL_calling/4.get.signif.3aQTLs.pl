#! /usr/bin/perl -w
use strict;

my @cts = qw/Ast Exc Inh Mic Oli Opc/;

foreach my $ct(@cts){
	# read in nominal p-value files
	open IN,"../6.FastQTL.permute.1M.withCov/$ct.permute.out.qvalue.signif.txt" or die $!;
	my %cut;
	<IN>;
	while(<IN>){
		chomp;
		my @sp = split /\t/;
		$cut{$sp[0]} = $sp[-1];
	}
	close IN;
	# read in sumstat files
	open IN2,"gzip -dc $ct.apaQTL.sumstat.txt.gz|" or die $!;
	open OUT,">$ct.apaQTL.signif.txt" or die $!;
	print OUT "phenotype_id\tvariant_id\tdistance\tma_samples\tma_count\tmaf\tpval_nominal\tslope\tslope_se\n";
	while(<IN2>){
		chomp;
		my @sp = split /\t/;
		next if !exists $cut{$sp[0]};
		next if $sp[6] > $cut{$sp[0]};
		print OUT "$_\n";
	}
	close IN2;
	close OUT;
}



