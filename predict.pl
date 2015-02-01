#! /usr/bin/env perl
use strict;
use warnings;

my ($config)=@ARGV;

die "Usage: $0 <config file>\n" if(@ARGV<1);

### Reading Params Section Start ###

my %config=&readconfig($config);

# simple parameters
my $program_base=&get_param("program_base");
my $outdir=&get_param("outdir");
my $ref=&get_param("ref");

# programs

my $parallel=&get_param("parallel");

my $augustus=&get_param("augustus");
my $augustus_training=&get_param("augustus_training");

my $genemark=&get_param("genemark");
my $genemark_training=&get_param("genemark_training");

my $glimmerhmm=&get_param("glimmerhmm");
my $glimmerhmm_training=&get_param("glimmerhmm_training");

my $geneid=&get_param("geneid");
my $geneid_training=&get_param("geneid_training");

my $snap=&get_param("snap");
my $snap_training=&get_param("snap_training");

# complex prarameters
my %protein=&get_param("protein","multi");

### Reading Params Section End ###


### Prepare Section Start ###

if(-e $outdir){
    print STDERR "Perform clean...\n";
    `rm -r $outdir`;
}
`mkdir $outdir`;
`mkdir $outdir/scaffolds`;
`mkdir $outdir/gff`;

`$program_base/tools/split_fasta.pl $ref $outdir/scaffolds`;

my @scaffolds=<$outdir/scaffolds/*.fa>;

### Prepare Section End ###



### Denovo Prediction Section Start ###

if($augustus && $augustus_training){
    &run_augustus($augustus,$augustus_training);
}

if($genemark && $genemark_training){
    &run_genemark($genemark,$genemark_training);
}

if($glimmerhmm && $glimmerhmm_training){
    &run_glimmerhmm($glimmerhmm,$glimmerhmm_training);
}

if($geneid && $geneid_training){
    &run_geneid($geneid,$geneid_training);
}

if($snap && $snap_training){
    &run_snap($snap,$snap_training);
}

### END OF PROGRAM ###


### Sub functions ###

sub run_augustus{
    my ($bin,$species)=@_;

    `mkdir $outdir/gff/augustus` if(!-e "$outdir/gff/augustus");

    foreach my $fa(@scaffolds){
        $fa=~/([^\/]+)\.fa$/;
        my $name=$1;
        print "$bin --species=$species $fa > $outdir/gff/augustus/$name.gff\n";
    }
}

sub run_genemark{
    my ($bin,$mtx)=@_;

    `mkdir $outdir/gff/genemark` if(!-e "$outdir/gff/genemark");

    foreach my $fa(@scaffolds){
        $fa=~/([^\/]+)\.fa$/;
        my $name=$1;
        print "$bin -m $mtx -o $outdir/gff/genemark/$name.gff $fa\n";
    }
}

sub run_glimmerhmm{
    my ($bin,$dir)=@_;

    `mkdir $outdir/gff/glimmerhmm` if(!-e "$outdir/gff/glimmerhmm");

    foreach my $fa(@scaffolds){
        $fa=~/([^\/]+)\.fa$/;
        my $name=$1;
        print "$bin $fa $dir -o $outdir/gff/glimmerhmm/$name.gff -g -f\n";
    }
}

sub run_geneid{
    my ($bin,$param)=@_;

    my $gff_dir="$outdir/gff/geneid";

    `mkdir $gff_dir` if(!-e $gff_dir);

    foreach my $fa(@scaffolds){
        $fa=~/([^\/]+)\.fa$/;
        my $name=$1;
        print "$bin -3 -P $param $fa > $gff_dir/$name.gff\n";
    }
}

sub run_snap{
    my ($bin,$hmm)=@_;

    my $gff_dir="$outdir/gff/snap";

    `mkdir $gff_dir` if(!-e $gff_dir);

    foreach my $fa(@scaffolds){
        $fa=~/([^\/]+)\.fa$/;
        my $name=$1;
        print "$bin $hmm $fa -gff > $gff_dir/$name.gff\n";
    }
}

sub get_param{
    my ($param,$type)=@_;
    if(!$type){
        my $detail="";
        foreach my $key(sort keys %config){
            next unless($key eq $param);
            my @value = keys $config{$key};
            $detail   = $config{$key}{$value[0]};
            last;
        }
        return($detail);
    }
    elsif($type eq "multi") {
        my %detail;
        foreach my $key(sort keys %config){
            next unless($key eq $param);
            my @value=keys $config{$key};
            foreach my $value(sort @value){
	my $result=$config{$key}{$value};
	$detail{$result}++;
            }
            last;
        }
        return(%detail);
    }
}

sub readconfig{
    my $config_file=shift;
    open(SUBI,"< $config_file");
    my %sub_config;
    my $line_no=0;
    while (<SUBI>) {
        chomp;
        next if(/^#/);
        next if(/^$/);
        if(/#/){
            s/#.*$//;
        }
        next unless(/^(\S+)\s+=\s+(\S+)/);
        $line_no++;
        $sub_config{$1}{$line_no}=$2;
    }
    close SUBI;
    return(%sub_config);
}
