#!/usr/bin/perl
#Eric Morrison
#112119

use strict;
use warnings;

sub count_hash{
	my($dirRef, $workDir, $file) = @_;
	my @dirs = @$dirRef;
	
	my %counts;
	foreach my $dir (@dirs){
		open(TMP, "$workDir/$dir/$file") || die "Can't open file\n";
		while(my $line = <TMP>){
			chomp($line);
			my @line = split("\t", $line);
			$counts{$file} += $line[2];
		}
	}
	return(\%counts);
}

sub print_counts{
	my $countRef = $_[0];
    my %counts = %$countRef;
    
    foreach my $sample (sort{$a cmp $b} keys %counts){
        print $sample, "\t", $counts{$sample}, "\n";
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
 
	my($countRef) = count_hash(\@dirs, $dir, $file);
	print_counts($countRef);
}
