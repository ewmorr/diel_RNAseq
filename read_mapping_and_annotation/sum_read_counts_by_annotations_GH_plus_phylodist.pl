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


sub process_anns_GH{
    my($ann) = @_;
    #hash indexed by seq name and GH category (because a gene can have multiple categories) with category as val; categories are semi-colon separated string with GH number and target name
    open(ANN, "$ann") || die "Can't open annotation file\n";
    my %ann;
    while(my$ann = <ANN>){
        chomp($ann);
        my @ann = split("\t", $ann);
        $ann{$ann[0]}{$ann[1]} = $ann[1].";".$ann[2];
    }
    return(\%ann);
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
    my %ann = %$annRef;
    my %ann2 = %$annRef2;

	#hash of annotation categories (unique names such as KO number or taxonomy
	#add cov to hash val by calling coverage hash as reading annotation file
	my %annCov;
	foreach my $covInd (keys %cov){
		
		#if either of the annotation categories is undefined or there are no reads mapped skip to next
		if(defined($ann{$covInd}) == 0 || defined($ann2{$covInd}) == 0 || ${ $cov{$covInd} }[1] == 0){next;}
		
        #foreach for GH ann categories in case there are multiple (which is few...)
        foreach my $catInd (keys %{ $ann{$covInd} }){
            my $categoryIndex = $ann{$covInd}{$catInd}.";".$ann2{$covInd};
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
	my $annotations = $ARGV[1];#gene annotations indexed by sequence ID *GH*
	my $annotations2 = $ARGV[2];#gene annotations indexed by sequence ID
	my $categoryInd2 = 4; #row number of desired category to summarise, this is harcoded for phylodist

	
	my $covHashRef = process_counts($coverage);
	my $annHashRef1 = process_anns_GH($annotations);
	my $annHashRef2 = process_anns($annotations2, $categoryInd2);
	my $annCovHashRef = add_cov_to_anns($annHashRef1, $annHashRef2, $covHashRef);
	print_anns($annCovHashRef);
}
