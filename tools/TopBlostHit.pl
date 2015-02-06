#!/usr/bin/perl
use strict;
use warnings;

my $blast2geneout=shift or die "perl $0 *bl2g\n";
my $out="$blast2geneout.tophit";
my %blast2gene;
my $pwd=`pwd`;chomp $pwd;
my $iterm=0;
open (F,"$blast2geneout")||die"$!";
while (<F>) {
    chomp;
    $iterm++;
    if (/^#GENE/){
        my $stand;
        if (/\~/){
            $stand='trev';
        }else{
            $stand='tfor';
        }
        s/\~\(//;
        s/\)//;
        /^#GENE\S+\s+coverage=\s+(\S+)\s+score=\s+(\S+)\s+identity=\s+(\S+)\s+(\d+)\.\.(\d+)$/;
        my $coverage=$1;
        my $score=$2;
        my $identity=$3;
        my $start=$4;
        my $end=$5;
        next if $coverage < 0.5;
        next if $identity < 0.4;
        my $line=<F>;
        chomp $line;
        die "wrong: $_\n$line\n" if $line=~/^#/;
        die "wrong: end<start\n" if $start>$end;
        my @line=split(/\s+/,$line);
        my $seq=$line[0];
        my $chr=$line[1];

        $blast2gene{$chr}{$seq}{$iterm}{start}=$start;
        $blast2gene{$chr}{$seq}{$iterm}{end}=$end;
        $blast2gene{$chr}{$seq}{$iterm}{score}=$score;
        $blast2gene{$chr}{$seq}{$iterm}{stand}=$stand;
    }
}
close F;

### one seq with multiple hits choose the top one, the priority was score. When scores were same, we will keep all hits for this seq ###
open (O,">$out")||die"$!";
my $ID=0;
print O "#ID\tChr\tSeq\tChr_start\tChr_end\tStand\n";
for my $chromosome (sort keys %blast2gene){
    for my $homolog (sort keys %{$blast2gene{$chromosome}}){
        my @mark=sort{$blast2gene{$chromosome}{$homolog}{$b}{score} <=> $blast2gene{$chromosome}{$homolog}{$a}{score}} keys %{$blast2gene{$chromosome}{$homolog}};
        my %besthit;
        $besthit{$mark[0]}++;
        for (my $i=0;$i<@mark;$i++){
            next if exists $besthit{$mark[$i]};
            for my $besthit (keys %besthit){
	if ($blast2gene{$chromosome}{$homolog}{$besthit}{score} < $blast2gene{$chromosome}{$homolog}{$mark[$i]}{score}){
	    delete $besthit{$besthit};
	    $besthit{$mark[$i]}++;
	}elsif($blast2gene{$chromosome}{$homolog}{$besthit}{score} == $blast2gene{$chromosome}{$homolog}{$mark[$i]}{score}){
	    $besthit{$mark[$i]}++;
	}
            }
        }
        for my $fltbest (sort keys %besthit){
            $ID++;
            print O "$ID\t$chromosome\t$homolog\t$blast2gene{$chromosome}{$homolog}{$fltbest}{start}\t$blast2gene{$chromosome}{$homolog}{$fltbest}{end}\t$blast2gene{$chromosome}{$homolog}{$fltbest}{stand}\n";
        }
    }
}
close O;
### done ###
