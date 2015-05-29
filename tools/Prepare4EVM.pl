#! /usr/bin/env perl
use strict;
use warnings;

my ($outdir,$weights)=@ARGV;

die "Usage: $0 <outdir> <weights file>\n" if(@ARGV<2);

my @scaffolds=<$outdir/scaffolds/*.fa>;

my @ab_initio=<$outdir/gff/ab_initio/*>;
my @homolog  =<$outdir/gff/homolog/*>;

&MergeGff4EVM("ab_initio",@ab_initio);
&MergeGff4EVM("homolog",@homolog);

my %input;
foreach my $ab_initio(@ab_initio){
    $ab_initio=~/([^\/]+)$/;
    my $type=$1;
    $input{ab_initio}{$type}=1;
}
foreach my $homolog(@homolog){
    $homolog=~/([^\/]+)$/;
    my $type=$1;
    $input{homolog}{$type}=1;
}

&writeWeightFile4EVM(\%input,$weights);

sub MergeGff4EVM{
    my $type=shift;
    my @dir=@_;

    foreach my $fa(@scaffolds){
        $fa=~/([^\/]+)\.fa$/;
        my $name=$1;
        `mkdir $outdir/temp/evm/dataOFscaffolds/$name` if(!-e "$outdir/temp/evm/dataOFscaffolds/$name");
        open(O,"> $outdir/temp/evm/dataOFscaffolds/$name/$type.gff");
        foreach my $dir(@dir){
            $dir=~/([^\/]+)$/;
            my $newsource=$1;
            next if(!-e "$dir/$name.gff");
            open(I,"< $dir/$name.gff") or die "Cannot open $dir/$name.gff\n";
            while (<I>) {
	chomp;
	next if(!/\S/);
	next if(/^#/);
	my @a=split(/\s+/);
	$a[1]=$newsource;
	print O join "\t",@a,"\n";
            }
            close I;
        }
        close O;
    }
}

sub writeWeightFile4EVM{
    my $input=shift;
    my $weights=shift;
    open(O,"> $weights") or die "Cannot create weights file: $weights\n";
    foreach my $type(sort keys %{$input{ab_initio}}){
        print O "ABINITIO_PREDICTION\t$type\t1\n";
    }
    foreach my $type(sort keys %{$input{homolog}}){
        print O "PROTEIN\t$type\t5\n";
    }
    close O;
}
