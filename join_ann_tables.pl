#!/usr/bin/perl
#Eric Morrison
#122118
#updated 11/22/19 to report sum avg cov instead of dividing by num genes
#
use strict;
use warnings;

sub ann_hash{
	my($dirRef, $workDir, $file) = @_;
	my @dirs = @$dirRef;
	
	my %annList;
	my %annDat;
	foreach my $dir (@dirs){
		open(TMP, "$workDir/$dir/$file") || die "Can't open annotation file\n";
		while(my $ann = <TMP>){
			chomp($ann);
			my @ann = split("\t", $ann);
			$annList{$ann[0]} = 1;
			$annDat{$ann[0]}{$dir}{"totalCov"} = $ann[1];
			$annDat{$ann[0]}{$dir}{"avgCov"} = $ann[2];
			$annDat{$ann[0]}{$dir}{"numGenes"} = $ann[3];
		}
	}
	#my @datRefs = ([%annList], [%annDat]);
	return(\%annList, \%annDat);
}

sub print_anns{
	my $annListRef = $_[0];
	my $annDatRef = $_[1];
	my $dirsRef = $_[2];
	my $file = $_[3];
	my %annList = %$annListRef;
	my %annDat = %$annDatRef;
	my @dirs = @$dirsRef;
	
	open(TOT, ">$file.totalCov.join") || die "can't open output file\n";
	open(AVG, ">$file.avgCov.join") || die "can't open output file\n";
	open(NUM, ">$file.numGenes.join") || die "can't open output file\n";

	print TOT "Category\t";
	print AVG "Category\t";
	print NUM "Category\t";

	foreach my $sample (@dirs){
		print TOT $sample, "\t";
		print AVG $sample, "\t";
		print NUM $sample, "\t";
	}
	print TOT "\n";
	print AVG "\n";
	print NUM "\n";
	
	foreach my $ann (sort{$a cmp $b} keys %annList){
		print TOT $ann, "\t";
		print AVG $ann, "\t";
		print NUM $ann, "\t";

		foreach my $sample (@dirs){
			if( defined( $annDat{$ann}{$sample} ) == 1 ){
				print TOT $annDat{$ann}{$sample}{"totalCov"}, "\t";
				print AVG $annDat{$ann}{$sample}{"avgCov"}, "\t";
				print NUM $annDat{$ann}{$sample}{"numGenes"}, "\t";
			}else{
				print TOT "0\t";
				print AVG "0\t";
				print NUM "0\t";
			}
		}
		print TOT "\n";
		print AVG "\n";
		print NUM "\n";
	}
}
#MAIN
{
	my $dir = $ARGV[0];#directory of dirs, dir names used as sample names
	my @dirs = `ls $dir`;
	chomp(@dirs);
	foreach my $dir (@dirs){
		print $dir, "\n";
	}
	my $file = $ARGV[1];#file within dirs that cintains the sample data

	my($annListRef, $annDatRef) = ann_hash(\@dirs, $dir, $file);
	print_anns($annListRef, $annDatRef, \@dirs, $file);
}
