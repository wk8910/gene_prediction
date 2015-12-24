#!/usr/bin/perl
use strict;
use warnings;

my $outdir=shift or die "perl $0 outdir\n";
my @gff=<$outdir/temp/evm/dataOFscaffolds/*/evm.out.gff>;
open (O,">$outdir/temp/evm/Merge.evm.out.gff");
for my $gff (@gff){
    open (F,"$gff")||die"$!";
    while (<F>) {
        chomp;
        next if /^\s*$/;
        print O "$_\n";
    }
    close F;
}
close O;
