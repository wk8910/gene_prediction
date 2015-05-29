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
    my %gff;

    # open (F,"<$rawgff")||die"$!";
    # open (O,">$newgff")||die"$!";
    # while (<F>) {
    #     chomp;
    #     my @a=split(/\t/,$_);
    #     if (/^\/\//){
    #         print O "\n";
    #         next;
    #     }
    #     next if /^\s*$/;
    #     next if $a[2] eq 'intron';
    #     $a[3]=$a[3]+$chrstart-1;
    #     $a[4]=$a[4]+$chrstart-1;
    #     if ($a[2] eq 'match'){
    #         next;
    #         $a[2]='mRNA' if $a[2] eq 'match';
    #         $a[8]="ID=$seq;"
    #     }elsif($a[2] eq 'cds'){
    #         $a[2]='CDS';
    #         #$a[8]="Parent=$seq;";
    #         $a[8]="ID=$seq;Parent=$seq;";
    #     }
    #     print O join("\t",@a),"\n";
    # }
    # close F;
    # close O;

    open(I,"< $rawgff") or die "Cannot open $rawgff\n";
    my $exon_id=0;
    while (<I>) {
        chomp;
        next if (/^\/\//);
        next if(/^#/);
        next unless(/cds/);
        my @a=split(/\t/);
        $a[3]=$a[3]+$chrstart-1;
        $a[4]=$a[4]+$chrstart-1;
        my ($chr,$source,$type,$start,$end,$score,$strand,$phase,$info)=@a;

        my $transcript_id=$seq;
        if($start>$end){
            my $tmp=$end;
            $end=$start;
            $start=$tmp;
        }

        $exon_id++;
        #my $transcript_id=$1;
        $gff{$chr}{$strand}{$transcript_id}{$exon_id}{source}=$source;
        $gff{$chr}{$strand}{$transcript_id}{$exon_id}{start} =$start;
        $gff{$chr}{$strand}{$transcript_id}{$exon_id}{end}   =$end;
        $gff{$chr}{$strand}{$transcript_id}{$exon_id}{type}  =$type;
        $gff{$chr}{$strand}{$transcript_id}{$exon_id}{score} =$score;
        $gff{$chr}{$strand}{$transcript_id}{$exon_id}{phase} =$phase;
    }
    close I;

    open(O,"> $newgff") or die "Cannot create $newgff\n";
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
}
