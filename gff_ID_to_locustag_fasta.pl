#!/usr/bin/perl
#Eric Morrison
#12/21/19
#changes "ID" from GFF to "locustag"

use strict;
use warnings;

sub hash_gff_names{
    my $file = $_[0];
    open(IN, "$file") || die "Can't open gff\n";
    chomp(my @file = <IN>);
    my %gff;
    foreach my $line (@file){
        my @line = split("\t", $line);
        $line[8] =~ /ID=(.*?);/;
        my $id = $1;
        $line[8] =~ /locus_tag=(.*?);/;
        my $ltag = $1;
        $id =~ s/\.//;
        $gff{$id} = $ltag;
    }
    return(\%gff);
}

sub replace_fasta_headers{
    my($gffRef, $fastaFile) = @_;
    my %gff = %$gffRef;
    open(IN, "$fastaFile") || die "Can't open fasta\n";
    chomp(my @fasta = <IN>);
    foreach my $line (@fasta){
        if($line =~ /^>/){
            $line =~ s/^>//;
            print ">", $gff{$line}, "\n";
        }else{
            print $line, "\n"
        }
    }
}
#MAIN
{
    my $gffFile = $ARGV[0];
    my $fastaFile = $ARGV[1];
    #my $out = $ARGV[2];
    
    my $gffRef = hash_gff_names($gffFile);
    replace_fasta_headers($gffRef, $fastaFile);
}
