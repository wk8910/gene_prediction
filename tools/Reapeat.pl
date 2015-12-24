#!/usr/bin/perl
use strict;
use warnings;


my ($name,$annotfa,$thread)=@ARGV;
die "Useage: perl $0 your_species_name your_species_assembly num_threads\n" if (scalar(@ARGV) != 3);
my $refspecies="all";
my ($trf,$repeatmask_dir,$repeatproteinmaske_dir,$repeatmodeler_dir)=("/home/share/users/yangyongzhi2012/tools/genome_annotation/TRF/trf407b.linux64","/home/share/users/yangyongzhi2012/tools/genome_annotation/RepeatMasker/RepeatMasker","/home/share/users/yangyongzhi2012/tools/genome_annotation/RepeatMasker-lastest/RepeatMasker","/home/share/users/yangyongzhi2012/tools/genome_annotation/RepeatMasker/RepeatModeler");

my $outdir="01.repeat";
`mkdir $outdir` if (! -e "$outdir");
my $trf_outdir="$outdir/TRF";
my $repeatmask_outdir="$outdir/RepeatMasker";
my $repeatproteinmasker="$outdir/RepeatProteinMasker";
my $repeatmodeler_outdir="$outdir/RepeatModeler";
$annotfa=~/([^\/]+)$/;
my $annotfaname=$1;
my $pwd=`pwd`;chomp $pwd;
`mkdir $trf_outdir` if (! -e "$trf_outdir");
`ln -s $pwd/$annotfa $trf_outdir/$annotfaname` if (! -e "$trf_outdir/$annotfaname");
`mkdir $repeatmask_outdir` if (! -e "$repeatmask_outdir");
`ln -s $pwd/$annotfa $repeatmask_outdir/$annotfaname` if (! -e "$repeatmask_outdir/$annotfaname");
`mkdir $repeatproteinmasker` if (! -e "$repeatproteinmasker");
`ln -s $pwd/$annotfa $repeatproteinmasker/$annotfaname` if (! -e "$repeatproteinmasker/$annotfaname");
`mkdir $repeatmodeler_outdir` if (! -e "$repeatmodeler_outdir");
`ln -s $pwd/$annotfa $repeatmodeler_outdir/$annotfaname` if (! -e "$repeatmodeler_outdir/$annotfaname");

open (O,">$0.sh");
## TRF ##
print O "cd $trf_outdir ; $trf $annotfaname 2 7 7 80 10 50 500 -d -h ; cd ../\n";
## RepeatMasker ##
print O "cd $repeatmask_outdir ; $repeatmask_dir/RepeatMasker -pa $thread -species $refspecies -nolow -norna -no_is -gff $annotfaname ; cd ../\n";
## RepeatProteinMasker ##
print O "cd $repeatproteinmasker ; $repeatproteinmaske_dir/RepeatProteinMask -engine abblast -noLowSimple -pvalue 1e-04 $annotfaname ; cd ../\n";
## RepeatModeler ##
print O "cd $repeatmodeler_outdir ; $repeatmodeler_dir/BuildDatabase -name $name $annotfaname 2>&1 | tee 01.BuildDatabase.log ; $repeatmodeler_dir/RepeatModeler -pa $thread -database $name 2>&1 | tee 02.RepeatModeler.log ; perl run.repeatmodeler.pl ; cd ../\n";
close O;
