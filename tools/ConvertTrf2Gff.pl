#!/usr/bin/perl
use strict;
use warnings;

my $in=shift or die "Usage: perl ConvertTrf2Gff.pl <the trf output with the -d -h parameters> <ouput file>\n";
my $out=shift or die "Usage: perl ConvertTrf2Gff.pl <the trf output with the -d -h parameters> <ouput file>\n";

open (F,"$in") or die "$!";
open (O,">$out") or die "$!";
print O "##gff-version 3\n";
my $sequence;
my $num=0;
while (<F>) {
    chomp;
    next if /^\s*$/;
    my @a=split(/\s+/,$_);
    if (/^Sequence:\s+(\S+)/){
        $sequence=$1;
    }
    if (scalar(@a)>=15 && $a[0]=~/^\d+$/ && $a[1]=~/^\d+$/ && $a[2]=~/^\d+$/){
        my $stand="+";
        $stand="-" if ($a[0]<$a[1]);
        $num++;
        my $outnum=sprintf("%05d","$num");
        print O "$sequence\tTRF\tTandemRepeat\t$a[0]\t$a[1]\t$a[7]\t$stand\t.\tID=TR$outnum;PeriodSize=$a[2];CopyNumber=$a[3];PercentMatches=$a[5];PercentIndels=$a[6];Consensus=$a[13];\n";
    }
}
close F;
close O;
