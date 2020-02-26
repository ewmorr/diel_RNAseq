#!/usr/bin/perl

use strict;
use warnings;

my $in = $ARGV[0];

open(IN, "$in") || die "can't open input\n";

chomp(my @in = <IN>);
my @header = split("\t", $in[0]);
shift @header;
print "Edge_num\tTax_name\tProtein\tTax_string\t", join("\t", @header), "\n";
shift @in;

foreach my $line (@in){
#    print $line, "\n";
    my @line = split("\t", $line);
    my @tax = split(";", $line[0]);
    shift @line;

    my $edge = join("\t", @tax[0 .. 2]);
    my $taxString = join(";", @tax[3 .. scalar(@tax)-1]);
    $edge =~ s/\'//g;
    $edge =~ s/\"//g;
    $edge =~ s/\(//g;
    $edge =~ s/\)//g;
    $edge =~ s/\#//g;
    $edge =~ s/=//g;
    print $edge, "\t", $taxString, "\t", join("\t", @line), "\n";
}

