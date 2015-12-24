#! /usr/bin/env perl
use strict;
use warnings;

my ($input,$output)=@ARGV;

my %gff;
open(I,"< $input") or die "Cannot open $input\n";
my $exon_id=0;
while (<I>) {
    chomp;
    next if(/^#/);
    next unless(/Exon/);
    my @a=split(/\t/);
    my ($chr,$source,$type,$start,$end,$score,$strand,$phase,$info)=@a;
    $source="CDS";

    if($start>$end){
        my $tmp=$end;
        $end=$start;
        $start=$tmp;
    }

    $exon_id++;
    my $transcript_id=$info;
    $gff{$chr}{$strand}{$transcript_id}{$exon_id}{source}=$source;
    $gff{$chr}{$strand}{$transcript_id}{$exon_id}{start} =$start;
    $gff{$chr}{$strand}{$transcript_id}{$exon_id}{end}   =$end;
    $gff{$chr}{$strand}{$transcript_id}{$exon_id}{type}  =$type;
    $gff{$chr}{$strand}{$transcript_id}{$exon_id}{score} =$score;
    $gff{$chr}{$strand}{$transcript_id}{$exon_id}{phase} =$phase;
}
close I;

open(O,"> $output") or die "Cannot create $output\n";
foreach my $chr(sort keys %gff){
    foreach my $strand(sort keys %{$gff{$chr}}){
        foreach my $transcript_id(sort keys %{$gff{$chr}{$strand}}){
            foreach my $exon_id(sort {$gff{$chr}{$strand}{$transcript_id}{$a}{start} <=> $gff{$chr}{$strand}{$transcript_id}{$b}{start}} keys %{$gff{$chr}{$strand}{$transcript_id}}){
	my $source=$gff{$chr}{$strand}{$transcript_id}{$exon_id}{source};
	my $start =$gff{$chr}{$strand}{$transcript_id}{$exon_id}{start};
	my $end   =$gff{$chr}{$strand}{$transcript_id}{$exon_id}{end};
	my $type  =$gff{$chr}{$strand}{$transcript_id}{$exon_id}{type};
	my $score =$gff{$chr}{$strand}{$transcript_id}{$exon_id}{score};
	my $phase =$gff{$chr}{$strand}{$transcript_id}{$exon_id}{phase};
	my $info="ID=$transcript_id;Parent=$transcript_id";
	my @line=($chr,$source,$type,$start,$end,$score,$strand,$phase,$info);
	print O join "\t",@line,"\n"
            }
        }
    }
}
close O;
