#! /usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;

my $file=shift;
my $outdir=shift;
`mkdir $outdir` if(!-e "$outdir");

my $in=Bio::SeqIO->new(-file=>$file,-format=>'fasta');
while(my $s=$in->next_seq){
    my $id=$s->id;
    my $seq=$s->seq;
    $id=~s/\|/_/g;
    my $outfile="./$outdir/$id.fa";
    open(OUT,"> $outfile")||die("Cannot creat $outfile!\n");
    print OUT ">$id\n$seq\n";
    close OUT;
}
