#! /usr/bin/env perl
use strict;
use warnings;

my $outdir=shift;

my $outfile="$outdir/prediction.gff";

my @gff=<$outdir/temp/evm/dataOFscaffolds/*/evm.out.new.gff>;

&merge_gff($outfile,@gff);

sub merge_gff{
    my $outfile=shift;
    my @gff=@_;

    open(O,"> $outfile");
    foreach my $gff(@gff){
        open(I,"< $gff");
        while (<I>) {
            print O "$_";
        }
        close I;
    }
    close O;
}
