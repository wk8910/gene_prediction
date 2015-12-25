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
my $transcript=&get_param("transcript");
my $species=&get_param("species");
my $gff2protein=&get_param("gff2protein");
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

my $evm=&get_param("evm");
my $evm2gff3=&get_param("evm2gff3");

my $pasa_dir=&get_param("pasa_dir");

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
`mkdir $outdir/temp/raw_gff`;
`mkdir $outdir/temp/raw_gff/ab_initio`;
`mkdir $outdir/temp/rna_seq`;

`mkdir $outdir/temp/evm`;
`mkdir $outdir/temp/evm/dataOFscaffolds`;

`$program_base/tools/split_fasta.pl $ref $outdir/scaffolds`;

my @scaffolds=<$outdir/scaffolds/*.fa>;

### Prepare Section End ###



### Denovo Prediction Section Start ###

my $run_denovo="$outdir/run/01.denovo.01.run.sh";

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

&run_job($run_denovo);

my $TransferFormat="$outdir/run/01.denovo.02.transfer.sh";
open(R,"> $TransferFormat");
&TransferFormat();
close R;

&run_job($TransferFormat);

### Denovo Prediction Section Done ###

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

### HOMOLOG PREDICTION SECTION END ###

### RNA_Seq PREDICTION SECTION START ###

if ($transcript && $pasa_dir){
    &run_rna_seq($transcript,$pasa_dir);
    &run_job("$outdir/run/03.rna_seq.sh");
}

### RNA_Seq PREDICTION SECTION DONE ###

### EVM SECTION START ###

my $merge_gff="$outdir/run/04.evm.01.prepare.sh";
open(R,"> $merge_gff");
my $weights="$outdir/temp/evm/weights.txt";
&MergeGff4EVM();
close R;

&run_job($merge_gff,"single");

my $run_evm="$outdir/run/04.evm.02.run.sh";
open(R,"> $run_evm");
&runEVM();
close R;

&run_job($run_evm);

my $collect_gff="$outdir/run/04.evm.03.collect.sh";
open (R,"> $collect_gff");
&collect_gff();
close R;

&run_job($collect_gff,"single");

### EVM SECTION END ###

### UPDATA ANNOTATION BY PASA START ###

if ($transcript && $pasa_dir){
    &unpdate_by_pasa($pasa_dir);
    &run_job("$outdir/run/05.pasa_update.sh","single");
}

### UPDATA ANNOTATION BY PASA DONE ###

### END OF PROGRAM ###


### Sub functions ###

sub run_augustus{
    my ($bin,$species)=@_;
    
    my $gff_dir="$outdir/temp/raw_gff/ab_initio/augustus";
    
    `mkdir $gff_dir` if(!-e $gff_dir);
    
    foreach my $fa(@scaffolds){
        $fa=~/([^\/]+)\.fa$/;
        my $name=$1;
        print R "$bin --species=$species $fa > $gff_dir/$name.gff\n";
    }
}

sub run_genemark{
    my ($bin,$mod)=@_;
    
    my $gff_dir="$outdir/temp/raw_gff/ab_initio/genemark";
    
    `mkdir $gff_dir` if(!-e $gff_dir);

    foreach my $fa(@scaffolds){
        $fa=~/([^\/]+)\.fa$/;
        my $name=$1;
        print R "$bin -f gff3 -m $mod -o $gff_dir/$name.gff $fa\n";
    }
}

sub run_glimmerhmm{
    my ($bin,$dir)=@_;

    my $gff_dir="$outdir/temp/raw_gff/ab_initio/glimmerhmm";

    `mkdir $gff_dir` if(!-e $gff_dir);

    foreach my $fa(@scaffolds){
        $fa=~/([^\/]+)\.fa$/;
        my $name=$1;
        print R "$bin $fa $dir -o $gff_dir/$name.gff -g -f\n";
    }
}

sub run_geneid{
    my ($bin,$param)=@_;

    my $gff_dir="$outdir/temp/raw_gff/ab_initio/geneid";

    `mkdir $gff_dir` if(!-e $gff_dir);

    foreach my $fa(@scaffolds){
        $fa=~/([^\/]+)\.fa$/;
        my $name=$1;
        print R "$bin -3 -P $param $fa > $gff_dir/$name.gff\n";
    }
}

sub run_snap{
    my ($bin,$hmm)=@_;

    my $gff_dir="$outdir/temp/raw_gff/ab_initio/snap";

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
            print BLSH "$blastall -p tblastn -d $fa -i $protein -e 1E-5 -o $blast_dir/$protein_name-$ref_name.tblastn $cpu; $program_base/tools/blast.parse.pl $blast_dir/$protein_name-$ref_name.tblastn $blast_dir/$protein_name-$ref_name.tblastn.bp ; $program_base/tools/blast2gene.pl $blast_dir/$protein_name-$ref_name.tblastn.bp > $blast_dir/$protein_name-$ref_name.tblastn.bp.bl2g ; $program_base/tools/TopBlostHit.pl $blast_dir/$protein_name-$ref_name.tblastn.bp.bl2g ; rm $blast_dir/$protein_name-$ref_name.tblastn $blast_dir/$protein_name-$ref_name.tblastn.bp $blast_dir/$protein_name-$ref_name.tblastn.bp.bl2g\n";
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
            my @value = keys %{$config{$key}};
            $detail   = $config{$key}{$value[0]};
            last;
        }
        return($detail);
    }
    elsif($type eq "multi") {
        my %detail;
        foreach my $key(sort keys %config){
            next unless($key eq $param);
            my @value=keys %{$config{$key}};
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

sub MergeGff4EVM{
    print R "perl $program_base/tools/Prepare4EVM.pl $outdir $weights yes yes ";
    if ($transcript){
        print R "yes\n";
    }
}

sub runEVM{
    foreach my $fa(@scaffolds){
        $fa=~/([^\/]+)\.fa$/;
        my $name=$1;
        my $ab_initio="$outdir/temp/evm/dataOFscaffolds/$name/ab_initio.gff";
        my $homolog="$outdir/temp/evm/dataOFscaffolds/$name/homolog.gff";
        my $transparament="";
        if ($transcript){
            $transparament="--transcript_alignments $outdir/temp/evm/dataOFscaffolds/$name/rna_seq.gff ";
        }
        print R "$evm --genome $fa --weights $weights --gene_predictions $ab_initio --protein_alignments $homolog $transparament> $outdir/temp/evm/dataOFscaffolds/$name/evm.out; $evm2gff3 $outdir/temp/evm/dataOFscaffolds/$name/evm.out $name > $outdir/temp/evm/dataOFscaffolds/$name/evm.out.gff3\n";
    }
}

sub run_job{
    my ($script,$type)=@_;
    if(!$type){
        $type="multi";
    }
    if($type eq "single"){
        if($thread_num){
            my $command="sh $script";
            print "$command\n";
            system($command);
        }
        else {
            print "sh $script\n";
        }
    }
    elsif($type eq "multi") {
        if($thread_num){
            my $command="$parallel -j $thread_num < $script";
            print "$command\n";
            system($command);
        }
        else {
            print "sh $script\n";
        }
    }
}

sub TransferFormat{
    foreach my $fa(@scaffolds){
        $fa=~/([^\/]+)\.fa$/;
        my $name=$1;
        if($augustus && $augustus_training){
            my $gff="$outdir/temp/raw_gff/ab_initio/augustus/$name.gff";
            `mkdir $outdir/gff/ab_initio/augustus` if(!-e "$outdir/gff/ab_initio/augustus");
            print R "perl $program_base/tools/TransferFormat_augustus.pl $gff $outdir/gff/ab_initio/augustus/$name.gff\n";
        }
        if($genemark && $genemark_training){
            my $gff="$outdir/temp/raw_gff/ab_initio/genemark/$name.gff";
            `mkdir $outdir/gff/ab_initio/genemark` if(!-e "$outdir/gff/ab_initio/genemark");
            print R "perl $program_base/tools/TransferFormat_genemark.pl $gff $outdir/gff/ab_initio/genemark/$name.gff\n";
        }
        if($glimmerhmm && $glimmerhmm_training){
            my $gff="$outdir/temp/raw_gff/ab_initio/glimmerhmm/$name.gff";
            `mkdir $outdir/gff/ab_initio/glimmerhmm` if(!-e "$outdir/gff/ab_initio/glimmerhmm");
            print R "perl $program_base/tools/TransferFormat_glimmerhmm.pl $gff $outdir/gff/ab_initio/glimmerhmm/$name.gff\n";
        }
        if($geneid && $geneid_training){
            my $gff="$outdir/temp/raw_gff/ab_initio/geneid/$name.gff";
            `mkdir $outdir/gff/ab_initio/geneid` if(!-e "$outdir/gff/ab_initio/geneid");
            print R "perl $program_base/tools/TransferFormat_geneid.pl $gff $outdir/gff/ab_initio/geneid/$name.gff\n";
        }
        if($snap && $snap_training){
            my $gff="$outdir/temp/raw_gff/ab_initio/snap/$name.gff";
            `mkdir $outdir/gff/ab_initio/snap` if(!-e "$outdir/gff/ab_initio/snap");
            print R "perl $program_base/tools/TransferFormat_snap.pl $gff $outdir/gff/ab_initio/snap/$name.gff\n";
        }
    }
}

sub collect_gff{
    print R "perl $program_base/tools/collect_gff.pl $outdir\n";
}

sub run_rna_seq{
    my ($trans,$pasadir)=@_;
    ## check pasa ##
    if (! "$pasadir/seqclean/seqclean/seqclean"){
        die "Not found compiled seqclean file in \"$pasadir/seqclean/seqclean/seqclean\"\n";
    }
    if (! "$pasadir/scripts/Launch_PASA_pipeline.pl"){
        die "Not found executable file in \"$pasadir/scripts/Launch_PASA_pipeline.pl\"\n";
    }
    ## check done ##
    open (RUN,">$outdir/run/03.rna_seq.sh");
    my $runpasadir="$outdir/temp/rna_seq";
    print RUN "cd $runpasadir 
$pasadir/seqclean/seqclean/seqclean $trans -v $pasadir/seqclean/seqclean/UniVec -c 1 -r Trinity.fasta.cln -o Trinity.fasta.clean | tee 01.seqclean.log
$pasadir/scripts/Launch_PASA_pipeline.pl -c  $pasadir/pasa_conf/alignAssembly.config -C -R -g $ref -t Trinity.fasta.clean --ALIGNERS blat,gmap --CPU 1 2>&1 | tee 02.pasa_assembly.log
cd ../../../
$program_base/split_gff.pl $pasadir $outdir/temp/rna_seq/rna_seq
";
}

sub unpdate_by_pasa{
    my ($pasadir)=@_;
    my $runpasadir="$outdir/temp/rna_seq";
    
    open (RUN,">$outdir/run/05.pasa_update.sh");
    print RUN "cd $runpasadir ; ln -s ../../prediction.gff orig_annotations_sample.gff3
$pasadir/scripts/Load_Current_Gene_Annotations.dbi -c $pasadir/pasa_conf/alignAssembly.config -g $ref -P orig_annotations_sample.gff3  2>&1 | tee 03.Loading_your_preexisting_protein-coding_gene_annotations.log
$pasadir/scripts/Launch_PASA_pipeline.pl -c $pasadir/pasa_conf/annotationCompare.config -A -g $ref -t Trinity.fasta.clean 2>&1 | tee 04.Compare.log
$program_base/tools/Sort.Redundancy.Rename.pl ./ $species.evmpasa.gff $species
$gff2protein $species.evmpasa.gff $ref prot > $species.evmpasa.pep
$gff2protein $species.evmpasa.gff $ref CDS > $species.evmpasa.cds
";
    close RUN;
}
