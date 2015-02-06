#!/usr/bin/perl
use strict;
use warnings;

my $outdir=shift or die "perl $0 \$outdir\n";
my @in=<$outdir/homolog/genewise/*/*/genewise.gff.transCoordinate.gff>;

my %diffout;
for my $in (@in){
    $in=~/^prediction_test\/homolog\/genewise\/([^-]+)-\S+\/\S+\/genewise.gff.transCoordinate.gff$/;
    my $name=$1;
    my $line1=`head -1 $in`;
    chomp $line1;
    my @line1=split(/\t/,$line1);
    my @start = sort{$a<=>$b} ($line1[3],$line1[4]);
    $diffout{$name}{$line1[0]}{$start[0]}{$line1[8]}=$in;
}

for my $k1 (sort keys %diffout){
    my $line=0;
    open (O,">$outdir/homolog/$k1.homolog.gff") || die "$!";
    for my $k2 (sort keys %{$diffout{$k1}}){
        for my $k3 (sort{$a<=>$b} keys %{$diffout{$k1}{$k2}}){
            for my $k4 (sort keys %{$diffout{$k1}{$k2}{$k3}}){
	$line++;
	open (F,"$diffout{$k1}{$k2}{$k3}{$k4}")||die"$!";
	while (<F>) {
	    print O "$_";
	}
	close F;
            }
        }
    }
    close O;
}

