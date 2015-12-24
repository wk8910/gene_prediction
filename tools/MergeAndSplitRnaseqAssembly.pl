#!/usr/bin/perl
use strict;
use warnings;

my $out=shift or die "Usage: MergeAndSplitRnaseqAssembly.pl dir\n";
my @rna=<$out/temp/rna_seq/*pasa_assemblies.gff3>;
if (scalar(@rna) > 1){
    die "Too many pasaGFF find. Please check!\n";
}
my $rna=$rna[0];
my $outdir="02.prediction_output/gff/rna_seq";
my %list;
my $line=0;
open (F,"$rna")||die"$!";
while (<F>) {
    chomp;
    $line++;
    my @a=split(/\t/,$_);
    $a[0]=~s/\|\S+//;
    $list{$a[0]}{$line}=join("\t",@a);
}
close F;

`mkdir $outdir` if (! -e "$outdir");

for my $k1 (sort keys %list){
    open (O,">$outdir/$k1.gff");
    for my $k2 (sort{$a<=>$b} keys %{$list{$k1}}){
        print O "$list{$k1}{$k2}\n";
    }
    close O;
}
