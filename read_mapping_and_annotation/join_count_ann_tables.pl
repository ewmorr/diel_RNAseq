#!/usr/bin/perl
#Eric Morrison
#111119

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
			$annDat{$ann[0]}{$dir}{"readCountLen"} = $ann[1];
			$annDat{$ann[0]}{$dir}{"readCount"} = $ann[2];
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
	
	open(LEN, ">$file.readCountRefLen.join") || die "can't open output file\n";
	open(COUNT, ">$file.readCount.join") || die "can't open output file\n";
	open(NUM, ">$file.readCountNumGenes.join") || die "can't open output file\n";

	print LEN "Category\t";
	print COUNT "Category\t";
	print NUM "Category\t";

	foreach my $sample (@dirs){
		print LEN $sample, "\t";
		print COUNT $sample, "\t";
		print NUM $sample, "\t";

	}
	print LEN "\n";
	print COUNT "\n";
	print NUM "\n";
	
	foreach my $ann (sort{$a cmp $b} keys %annList){
		print LEN $ann, "\t";
		print COUNT $ann, "\t";
		print NUM $ann, "\t";

		foreach my $sample (@dirs){
			if( defined( $annDat{$ann}{$sample} ) == 1 ){
				print LEN $annDat{$ann}{$sample}{"readCountLen"}, "\t";
				print COUNT $annDat{$ann}{$sample}{"readCount"}, "\t";
				print NUM $annDat{$ann}{$sample}{"numGenes"}, "\t";
			}else{
				print LEN "0\t";
				print COUNT "0\t";
				print NUM "0\t";
			}
		}
		print LEN "\n";
		print COUNT "\n";
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
	my $file = $ARGV[1];#file within dirs that contains the sample data

	my($annListRef, $annDatRef) = ann_hash(\@dirs, $dir, $file);
	print_anns($annListRef, $annDatRef, \@dirs, $file);
}
