#!/usr/bin/perl
use strict;
use warnings;

my $in=shift or die "Usage: perl ConverttRNAscan-SE2Gff.pl <tRNA.secondary.out> <out.gff>\n";
my $out=shift or die "Usage: perl ConverttRNAscan-SE2Gff.pl <tRNA.secondary.out> <out.gff>\n";

open (F,"$in");
open (O,">$out");
print O "##gff-version 3\n";
my $num=0;
while (<F>) {
    my $l1=$_; chomp $l1;
    my $l2=<F>;chomp $l2;
    while (<F>) {
        chomp;
        last if $_=~/^Seq\:\s+\S+$/;
    }
    my $l5=<F>;chomp $l5;
    <F>;
    $num++;
    my $outnum=sprintf("%04d","$num");
    $l1=~/^(\S+)\.trna\d+\s+\((\d+)\-(\d+)\)/ or die "wrong type: $l1\n";
    my ($chr,$start,$end)=($1,$2,$3);
    my $strand="+";
    $strand="-" if $start>$end;
    $l2=~/^Type\:\s+(\S+)\s+Anticodon\:\s+(\S+).*Score\:\s+(\S+)/ or die "wrong type: $l2\n";
    my ($type,$anticodon,$score)=($1,$2,$3);
    $l5=~/^Str\:\s+(\S+)$/ or die "wrong type: $l5\n";
    my $str=$1;
    
    print O "$chr\ttRNAscan-SE\ttRNA\t$start\t$end\t$score\t$strand\t.\tID=tRNA_$outnum;Type=$type;Anti-codon=$anticodon;structure=\"$str\";\n";
}
close O;
