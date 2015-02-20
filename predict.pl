#! /usr/bin/env perl
use strict;
use warnings;

my ($config,$thread_num)=@ARGV;

die "Usage: $0 <config file> [number of threads]\n" if(@ARGV<1);

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

my $formatdb=&get_param("formatdb");
my $blastall=&get_param("blastall");

my $genewise=&get_param("genewise");

# complex prarameters
my %protein=&get_param("protein","multi");

### Reading Params Section End ###


### Prepare Section Start ###

if(-e $outdir){
    print STDERR "Perform clean...\n";
    `rm -r $outdir`;
}
`mkdir $outdir`;
`mkdir $outdir/run`;

`mkdir $outdir/scaffolds`;
`mkdir $outdir/gff`;
`mkdir $outdir/gff/ab_initio`;
`mkdir $outdir/gff/homolog`;

`mkdir $outdir/temp`;
`mkdir $outdir/temp/homolog/`;

`$program_base/tools/split_fasta.pl $ref $outdir/scaffolds`;

my @scaffolds=<$outdir/scaffolds/*.fa>;

### Prepare Section End ###



### Denovo Prediction Section Start ###

my $run_denovo="$outdir/run/01.denovo.sh";

open(R,"> $run_denovo") or die "Cannot write $run_denovo !\n";

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

close R;

if($thread_num){
    my $command="$parallel -j $thread_num < $run_denovo";
    print "$command\n";
    system($command);
}
else {
    print "sh $run_denovo\n";
}

### Homolog Prediction Section Start ###


if (scalar(keys %protein)>0){
    if ($formatdb && $blastall && $genewise){
        &run_homolog($formatdb,$blastall,$genewise);

        if($thread_num){
            my $command1="sh $outdir/run/02.homolog.01.blast.sh";
            print "$command1\n";
            system($command1);

            my $command2="$parallel -j $thread_num < $outdir/run/02.homolog.02.genewise_input.sh";
            print "$command2\n";
            system($command2);

            my $command3="$parallel -j $thread_num < $outdir/run/02.homolog.03.genewise_run.sh";
            print "$command3\n";
            system($command3);

            my $command4="sh $outdir/run/02.homolog.04.genewise_output.sh";
            print "$command4\n";
            system($command4);
        }else{
            print "sh $outdir/run/02.homolog.01.blast.sh
sh $outdir/run/02.homolog.02.genewise_input.sh
sh $outdir/run/02.homolog.03.genewise_run.sh
sh $outdir/run/02.homolog.04.genewise_output.sh
";
        }
    }
}

### END OF PROGRAM ###


### Sub functions ###

sub run_augustus{
    my ($bin,$species)=@_;
    
    my $gff_dir="$outdir/gff/ab_initio/augustus";
    
    `mkdir $gff_dir` if(!-e $gff_dir);
    
    foreach my $fa(@scaffolds){
        $fa=~/([^\/]+)\.fa$/;
        my $name=$1;
        print R "$bin --species=$species $fa > $gff_dir/$name.gff\n";
    }
}

sub run_genemark{
    my ($bin,$mtx)=@_;
    
    my $gff_dir="$outdir/gff/ab_initio/genemark";
    
    `mkdir $gff_dir` if(!-e $gff_dir);

    foreach my $fa(@scaffolds){
        $fa=~/([^\/]+)\.fa$/;
        my $name=$1;
        print R "$bin -m $mtx -o $gff_dir/$name.gff $fa\n";
    }
}

sub run_glimmerhmm{
    my ($bin,$dir)=@_;

    my $gff_dir="$outdir/gff/ab_initio/glimmerhmm";

    `mkdir $gff_dir` if(!-e $gff_dir);

    foreach my $fa(@scaffolds){
        $fa=~/([^\/]+)\.fa$/;
        my $name=$1;
        print R "$bin $fa $dir -o $gff_dir/$name.gff -g -f\n";
    }
}

sub run_geneid{
    my ($bin,$param)=@_;

    my $gff_dir="$outdir/gff/ab_initio/geneid";

    `mkdir $gff_dir` if(!-e $gff_dir);

    foreach my $fa(@scaffolds){
        $fa=~/([^\/]+)\.fa$/;
        my $name=$1;
        print R "$bin -3 -P $param $fa > $gff_dir/$name.gff\n";
    }
}

sub run_snap{
    my ($bin,$hmm)=@_;

    my $gff_dir="$outdir/gff/ab_initio/snap";

    `mkdir $gff_dir` if(!-e $gff_dir);

    foreach my $fa(@scaffolds){
        $fa=~/([^\/]+)\.fa$/;
        my $name=$1;
        print R "$bin $hmm $fa -gff > $gff_dir/$name.gff\n";
    }
}

sub run_homolog{
    my ($db,$blast,$bin)=@_;
    for my $fa (@scaffolds){
        `$db -i $fa -p F`;
    }
    my $blast_dir="$outdir/temp/homolog/blast2gene";
    `mkdir $blast_dir` if (! -e "$blast_dir");
    my $genewise_dir="$outdir/temp/homolog/genewise";
    `mkdir $genewise_dir` if (! -e "$genewise_dir");
    open (BLSH,">$outdir/run/02.homolog.01.blast.sh");
    open (INPUT,">$outdir/run/02.homolog.02.genewise_input.sh");
    #open (GENEWISE,">$outdir/run/02.homolog.03.genewise_run.sh");
    for my $protein (sort keys %protein){
        $protein=~/([^\/]+)$/;
        my $protein_name=$1;
        for my $fa (@scaffolds){
            $fa=~/([^\/]+)\.fa$/;
            my $ref_name=$1;
            my $cpu="";
            $cpu="-a $thread_num " if ($thread_num);
            print BLSH "$blastall -p tblastn -d $fa -i $protein -e 1E-5 -o $blast_dir/$protein_name-$ref_name.tblastn $cpu; $program_base/tools/blast.parse.pl $blast_dir/$protein_name-$ref_name.tblastn $blast_dir/$protein_name-$ref_name.tblastn.bp ; $program_base/tools/blast2gene.pl $blast_dir/$protein_name-$ref_name.tblastn.bp > $blast_dir/$protein_name-$ref_name.tblastn.bp.bl2g ; $program_base/tools/TopBlostHit.pl $blast_dir/$protein_name-$ref_name.tblastn.bp.bl2g\n";
            print INPUT "$program_base/tools/genewiseINPUT.pl $outdir $blast_dir $blast_dir/$protein_name-$ref_name.tblastn.bp.bl2g.tophit $protein $fa $bin $program_base/tools/TransCoordinate.pl\n";
        }
    }
    close BLSH;
    close INPUT;
    open (GENEWISEOUT,">$outdir/run/02.homolog.04.genewise_output.sh");
    print GENEWISEOUT "$program_base/tools/MergeSplitGff.pl $outdir\n";
    close GENEWISEOUT;
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
    open(SUBI,"< $config_file") or die "Cannot open configuration file: $config_file\n";
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
