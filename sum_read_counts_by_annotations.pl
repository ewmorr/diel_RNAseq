#!/usr/bin/perl
#Eric Morrison
#111919

use strict;
use warnings;

sub process_counts{
	my $coverage = $_[0];
	open(COV, "$coverage") || die "can't open coverage file\n";
	my %cov;
	#hash of arrays with seq ID index and [total bases mapped, tot. mapped/seq. length]
	while(my $cov = <COV>){
		chomp($cov);
		my @cov = split("\t", $cov);
		$cov{$cov[0]} = [$cov[1], $cov[2]];
	}
	return(\%cov);
}

sub process_anns{
	my($ann, $categoryInd) = @_;
	#hash indexed by seq name with category (e.g., phylodist string, KO number), as val
	open(ANN, "$ann") || die "Can't open annotation file\n";
	my %ann;
	while(my$ann = <ANN>){
		chomp($ann);
		my @ann = split("\t", $ann);
		$ann{$ann[0]} = $ann[$categoryInd];
	}
	return(\%ann);
}

sub add_cov_to_anns{
	my($annRef, $annRef2, $covRef) = @_;
	my %cov = %$covRef;
	#hash of annotation categories (unique names such as KO number or taxonomy
	#add cov to hash val by calling coverage hash as reading annotation file
	my %annCov;
	foreach my $covInd (keys %cov){
		
		#if either of the annotation categories is undefined skip to next
		if(defined($$annRef{$covInd}) == 0 || defined($$annRef2{$covInd}) == 0){next;}
		
		my $categoryIndex = $$annRef{$covInd}.";".$$annRef2{$covInd};
		if(defined($annCov{$categoryIndex} ) == 0){
			$annCov{$categoryIndex}{"len"} = ${ $cov{$covInd} }[0];
			$annCov{$categoryIndex}{"readCov"} = ${ $cov{$covInd} }[1];
			$annCov{$categoryIndex}{"count"} = 1;
		}elsif(defined($annCov{$categoryIndex}) == 1){
			$annCov{$categoryIndex}{"len"} += ${ $cov{$covInd} }[0];
			$annCov{$categoryIndex}{"readCov"} += ${ $cov{$covInd} }[1];
			$annCov{$categoryIndex}{"count"}++;
		}
	}
	return(\%annCov);
}

sub print_anns{
	my $annRef = $_[0];
	my %anns = %$annRef;
	
	foreach my $function (sort{$a cmp $b} keys %anns){
		print $function, "\t", $anns{$function}{"len"}, "\t", $anns{$function}{"readCov"}, "\t", $anns{$function}{"count"}, "\n";
	}
}
#MAIN
{
	my $coverage = $ARGV[0];#Mapped read counts by gene sequence from samtools idxstats
	my $annotations = $ARGV[1];#gene annotations indexed by sequence ID
	my $categoryInd = $ARGV[2]; #row number of desired category to summarise, for example KO number or taxonomy, (starting from row 0)
	my $annotations2 = $ARGV[3];#gene annotations indexed by sequence ID
	my $categoryInd2 = $ARGV[4]; #row number of desired category to summarise, for example KO number or taxonomy, (starting from row 0)

	
	my $covHashRef = process_counts($coverage);
	my $annHashRef1 = process_anns($annotations, $categoryInd);
	my $annHashRef2 = process_anns($annotations2, $categoryInd2);
	my $annCovHashRef = add_cov_to_anns($annHashRef1, $annHashRef2, $covHashRef);
	print_anns($annCovHashRef);
}
