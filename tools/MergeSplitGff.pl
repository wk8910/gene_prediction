#!/usr/bin/perl
use strict;
use warnings;

my $outdir=shift or die "perl $0 \$outdir\n";
my @in=<$outdir/temp/homolog/blast2gene/*/*/genewise.gff.transCoordinate.gff>;

my %diffout;
for my $in (@in){
    $in=~/temp\/homolog\/blast2gene\/([^-]+)-\S+\/\S+\/genewise.gff.transCoordinate.gff$/;
    my $name=$1;
    my $line1=`head -1 $in`;
    chomp $line1;
    my @line1=split(/\t/,$line1);
    next if(@line1 != 9);
    my @start = sort{$a<=>$b} ($line1[3],$line1[4]);
    $diffout{$name}{$line1[0]}{$start[0]}{$line1[8]}=$in;
}

for my $species (sort keys %diffout){
    my $line=0;
    `mkdir $outdir/gff/homolog/$species` if(!-e "$outdir/gff/homolog/$species");
    for my $scaffold (sort keys %{$diffout{$species}}){
        open(O,"> $outdir/gff/homolog/$species/$scaffold.gff");
        for my $position (sort{$a<=>$b} keys %{$diffout{$species}{$scaffold}}){
            for my $name (sort keys %{$diffout{$species}{$scaffold}{$position}}){
	$line++;
	open (F,"$diffout{$species}{$scaffold}{$position}{$name}") or die "$!";
	while (<F>) {
	    print O "$_";
	}
	close F;
            }
        }
        close O;
    }
}

