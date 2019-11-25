#!/usr/bin/perl

use strict;
use warnings;

my $in = $ARGV[0];

open(IN, "$in") || die "can't open input\n";

chomp(my @in = <IN>);
my @header = split("\t", $in[0]);
shift @header;
print "Kindom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\t", join("\t", @header), "\n";
shift @in;

foreach my $line (@in){
#    print $line, "\n";
    my @line = split("\t", $line);
    my @tax = split(";", $line[0]);
    shift @line;

    my $tax = join("\t", @tax[0 .. 6]);
    $tax =~ s/\'//g;
    $tax =~ s/\"//g;
    $tax =~ s/\(//g;
    $tax =~ s/\)//g;
    $tax =~ s/\#//g;
    print $tax, "\t", join("\t", @line), "\n";
}

