#! /usr/bin/env perl
use strict;
use warnings;

my ($config)=@ARGV;

die "Usage: $0 <config file>\n" if(@ARGV<1);

### Reading Params Section Start ###

my %config=&readconfig($config);

# simple parameters
my $program_base=&get_param("program_base","single");
my $outdir=&get_param("outdir","single");
my $ref=&get_param("ref","single");

# programs
my $augustus=&get_param("augustus","single");
my $parallel=&get_param("parallel","single");
my $genemark=&get_param("genemark","single");
my $genemark_mtx=&get_param("genemark_mtx","single");

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

if($augustus){
    my $species="human";
    &run_augustus($augustus,$species);
}

# if(!$genemark && $genemark_mtx){
#     &run_genemark($genemark,$genemark_mtx);
# }

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

sub get_param{
    my ($param,$type)=@_;
    if($type eq "single"){
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
