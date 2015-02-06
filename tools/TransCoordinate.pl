#!/usr/bin/perl
use strict;
use warnings;

my $genewisegff=shift or die "perl $0 genewise.gff\n";
die "no genewise output\n" if (! -e "$genewisegff");
my $subseqtxt=$genewisegff;
$subseqtxt=~s/genewise.gff$/subseq.txt/;
&Coordinate("$subseqtxt","$genewisegff","$genewisegff.transCoordinate.gff");

sub Coordinate{
    my ($extract,$rawgff,$newgff)=@_;
    ### check ###
    for my $check ($extract,$rawgff){
        die "$check doesn't exists\n" if (! -e "$check");
    }
    my $lines=`wc -l $extract`;
    chomp $lines;
    $lines=~s/^(\d+)\s+.*/$1/;
    die "wrong format: lines large than 2 in the $extract\n" if ($lines != 2);
    ### check done ###
    my @out;
    my $lastline=`tail -1 $extract`;
    chomp $lastline;
    $lastline=~/^\S+\t(\S+)\t(\S+)\t\d+\t\d+\t(\d+)\t(\d+)$/ or die "wrong format: $extract\n";
    my ($chr,$seq,$chrstart,$chrend)=($1,$2,$3,$4);
    open (F,"<$rawgff")||die"$!";
    open (O,">$newgff")||die"$!";
    while (<F>) {
        chomp;
        my @a=split(/\t/,$_);
        if (/^\/\//){
            print O "\n";
            next;
        }
        next if /^\s*$/;
        next if $a[2] eq 'intron';
        $a[3]=$a[3]+$chrstart-1;
        $a[4]=$a[4]+$chrstart-1;
        if ($a[2] eq 'match'){
            $a[2]='mRNA' if $a[2] eq 'match';
            $a[8]="ID=$seq;"
        }elsif($a[2] eq 'cds'){
            $a[2]='CDS';
            $a[8]="Parent=$seq;";
        }
        print O join("\t",@a),"\n";
    }
    close F;
    close O;
}
