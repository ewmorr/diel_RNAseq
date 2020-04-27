#!/usr/bin/perl
#Eric Morrison
#This script searches the KO reference hierarchy (file retrieved from https://www.genome.jp/kegg-bin/get_htext?ko00001.keg) for K number and returns the associated names (quoted for upload to R).

use strict;
use warnings;

sub hash_kref{
    my $Kref = $_[0];
    open(KREF, "$Kref") || die "can't open KO refrence hierarchy file";
    chomp(my @Kref = <KREF>);
    my %Kref;
    foreach my $ko (@Kref){
        if($ko =~ /^D/){
            $ko =~ s/^D\s{6}//;
            my @ko = split('\s\s', $ko);
            $Kref{$ko[0]} = $ko[1];
        }
    }
    return(\%Kref);
}

sub print_names{
    my($Kref, $numList) = @_;
    my %Kref = %$Kref;
    open(NUM, "$numList") || die "Can't open search file\n";
    chomp(my @numList = <NUM>);
    foreach my $knum (@numList){
        my $num = $knum;
        $num =~ s/KO://;
        if(defined($Kref{$num}) == 1){
            print $knum, "\t\"", $Kref{$num}, "\"\n";
        }else{
            print $knum, "\t", "NA", "\n";
        }
     }
}


#MAIN
{
    my $Kref = $ARGV[0];
    my $numList = $ARGV[1];
    my $KrefHash = hash_kref($Kref);
    print_names($KrefHash, $numList);
}
