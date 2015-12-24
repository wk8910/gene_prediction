#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;

my ($in,$genome,$out)=@ARGV;
die "Usage:\n$0 <gffINputFile> <genomeSequence> <outfile>\n" if (! $out);

my %rep;
my $line=0;
print "begin reading gff\n";
open (F,"$in")||die"$!";
while (<F>) {
    chomp;
    next if /^#/;
    $line++;
    my @a=split(/\t/,$_);
    if ($a[2] =~/TEprotein|Transposon/i){
        my @pos=sort{$a<=>$b} ($a[3],$a[4]);
        $rep{$a[0]}{$line}{start}=$pos[0];
        $rep{$a[0]}{$line}{end}=$pos[1];
    }
}
close F;
print "begin reading fasta\n";
my %seq;
my $fa=Bio::SeqIO->new(-file=>"$genome",-format=>"fasta");
while (my $seq=$fa->next_seq) {
    my $id=$seq->id;
    my $seq=$seq->seq;
    $seq{$id}=$seq;
}

open (O,">$out");
for my $k1 (sort keys %seq){
    print O ">$k1\n";
    my $seq=$seq{$k1};
    my @seq=split(//,$seq);
    if (exists $rep{$k1}){
        my %site;
        for my $k2 (sort keys %{$rep{$k1}}){
            my ($start,$end)=($rep{$k1}{$k2}{start},$rep{$k1}{$k2}{end});
            for (my $i=$start;$i<$end;$i++){
	$site{$i}++;
            }
        }
        for (my $j=0;$j<@seq;$j++){
            my $k=$j+1;
            $seq[$j]='N' if exists $site{$k};
        }
        print O join("",@seq),"\n";
    }else{
        print O "$seq\n";
    }
}
close O;
