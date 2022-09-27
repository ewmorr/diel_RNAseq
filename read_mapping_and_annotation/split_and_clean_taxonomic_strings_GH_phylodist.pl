#!/usr/bin/perl

use strict;
use warnings;

my $in = $ARGV[0];

open(IN, "$in") || die "can't open input\n";

chomp(my @in = <IN>);
my @header = split("\t", $in[0]);
shift @header;
print "GH\tKindom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\t", join("\t", @header), "\n";
shift @in;

foreach my $line (@in){
#    print $line, "\n";
    my @line = split("\t", $line);
    my @tax = split(";", $line[0]);
    shift @line;
    my $gh = join(";", @tax[0 ..1]);
    my $tax = join("\t", @tax[2 .. 8]);
    $tax =~ s/\'//g;
    $tax =~ s/\"//g;
    $tax =~ s/\(//g;
    $tax =~ s/\)//g;
    $tax =~ s/\#//g;
    print $gh, "\t", $tax, "\t", join("\t", @line), "\n";
}

