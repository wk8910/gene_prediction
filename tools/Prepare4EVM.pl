#! /usr/bin/env perl
use strict;
use warnings;

my ($outdir,$weights,$ab_initio_stat,$homolog_stat,$rna_seq_stat)=@ARGV;
die "perl $0 \$outdir \$weights \$ab_initio_stat \$homolog_stat \$rna_seq_stat\nexample: perl $0 02.prediction_output 02.prediction_output/temp/evm/weights.txt yes yes no" if (! $rna_seq_stat);
my %input;
my @scaffolds=<$outdir/scaffolds/*.fa>;
if ($ab_initio_stat eq 'yes'){
    my @ab_initio=<$outdir/gff/ab_initio/*>;
    &MergeGff4EVM("ab_initio",@ab_initio);
    &Fillinput("ab_initio",@ab_initio);
}
if ($homolog_stat eq 'yes'){
    my @homolog=<$outdir/gff/homolog/*>;
    &MergeGff4EVM("homolog",@homolog);
    &Fillinput("homolog",@homolog);
}
if ($rna_seq_stat eq 'yes'){
    my @rna_seq=<$outdir/gff/rna_seq/*>;
    &MergeGff4EVM("rna_seq",@rna_seq);
    &Fillinput("rna_seq",@rna_seq);
}

&writeWeightFile4EVM(\%input,$weights);

sub Fillinput{
    my $type=shift;
    my @dir=@_;
    for my $dir (@dir){
        $dir=~/([^\/]+)$/;
        my $software=$1;
        $input{$type}{$software}=1;
    }
}
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
    my %method2print;
    $method2print{ab_initio}{print}="ABINITIO_PREDICTION";
    $method2print{ab_initio}{score}="1";
    $method2print{homolog}{print}="PROTEIN";
    $method2print{homolog}{score}="5";
    $method2print{rna_seq}{print}="TRANSCRIPT";
    $method2print{rna_seq}{score}="10";
    
    open(O,"> $weights") or die "Cannot create weights file: $weights\n";
    foreach my $method (sort keys %input){
        foreach my $type(sort keys %{$input{$method}}){
            print O "$method2print{$method}{print}\t$type\t$method2print{$method}{score}\n";
        }
    }
    close O;
}
