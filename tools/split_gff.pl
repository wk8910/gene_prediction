#!/usr/bin/perl
use strict;
use warnings;

my ($indir,$outdir)=@ARGV;
die "Usage:\nperl split_gff.pl <containing pasa assembled gff files dir> <output dir>\n" if (! $outdir);
`mkdir -p $outdir` if (! -e "$outdir");
my %h;
my @in=<$indir/*pasa_assemblies.gff3>;
if (scalar(@in) != 1){
    die "None or more than one pass assembled gff files, please cotain one file\n";
}
my $in=$in[0];
my $i=0;
open (F,"$in")||die"$!";
while (<F>) {
    chomp;
    my @a=split(/\t/,$_);
    $a[1]='rna_sq';
    $i++;
    $h{$a[0]}{$i}=join("\t",@a);
}
close F;

for my $k1 (sort keys %h){
    open (O,">$outdir/$k1.gff");
    for my $k2 (sort{$a<=>$b} keys %{$h{$k1}}){
        print O "$h{$k1}{$k2}\n";
    }
    close O;
}
