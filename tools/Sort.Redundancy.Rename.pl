#!/usr/bin/perl
use strict;
use warnings;

my ($ingffdir,$outgff,$species)=@ARGV;
die "Usage: 01.fltXandRedundancy.pl <ingff (evm out gff)> <outgff> <species name>\n" if (! $species);
my @ingff=<$ingffdir/*gene_structures_post_PASA_updates.*gff3>;
if (scalar(@ingff) != 1){
    die "please contain one *gene_structures_post_PASA_updates.*gff3 file in the dir: $ingffdir\n";
}
my $ingff=$ingff[0];
my %genetype;
my %trans2gene;
my %gene2trans;
my %fail;
open (F,"$ingff")||die"$!";
while (<F>) {
    chomp;
    next unless /^#/;
    if (/^#\s+ORIGINAL:\s+(\S+)/){
        my $t0=$1;
        $genetype{$t0}="gene";
    }elsif(/^#\s+PASA_UPDATE:\s+(\S+),/){
        my $t1=$1;
        if (/single\s+gene\s+model\s+update/){
            $genetype{$t1}="gene";
        }elsif(/alt-splice\s+addition/){
            $genetype{$t1}="alt";
        }else{
            die "$_\n";
        }
    }elsif(/^#PROT\s+(\S+)\s+(\S+)\s+(\S+)/){
        my ($t,$g,$s)=($1,$2,$3);
        $fail{$t}++;
        my $checklen=length($s);
        if ($checklen < 50){
            $fail{$t}++;
        }
        $trans2gene{$t}{$g}++;
        $gene2trans{$g}{$t}++;
    }else{
        die "$_\n";
    }
}
close F;
my %out;
for my $k1 (sort keys %trans2gene){
    next if exists $fail{$k1};
    my $type=$genetype{$k1};
    my @gene=sort keys %{$trans2gene{$k1}};
    die "one trans to many genes: $k1\n" if (@gene>1);
    my $gene=$gene[0];
    if ($type eq 'alt'){
        my @trans=sort keys %{$gene2trans{$gene}};
        if (@trans == 1){
            #print "$trans[0]\t$gene\n";
            for my $g (sort keys %gene2trans){
	my @g=split(/_/,$g);
	next if @g < 2;
	my $check=0;
	for (my $i=0;$i<@g;$i++){
	    $check++ if ($g[$i] eq $gene);
	}
	if ($check>0){
	    $out{gene2trans}{$g}{$k1}++;
	    $out{trans2gene}{$k1}{$g}++;
	}
            }
        }else{
            $out{gene2trans}{$gene}{$k1}++;
            $out{trans2gene}{$k1}{$gene}++;
        }
    }else{
        $out{gene2trans}{$gene}{$k1}++;
        $out{trans2gene}{$k1}{$gene}++;
    }
}

undef (%genetype);
undef (%trans2gene);
undef (%gene2trans);
undef (%fail);
my %gff;
my $line=0;
open (F,"$ingff")||die"$!";
while (<F>) {
    chomp;
    next if /^#/;
    next if /^\s*$/;
    my @a=split(/\t/,$_);
    $line++;
    if ($a[2] eq 'gene'){
        $a[8]=~/ID=([^;]+);/ or die "$_\n";
        $gff{gene}{$1}=$_;
    }elsif($a[2] eq 'mRNA'){
        $a[8]=~/ID=([^;]+);/ or die "$_\n";
        $gff{mrna}{$1}{0}=$_;
    }else{
        $a[8]=~/Parent=([^;]+)/ or die "$_\n";
        $gff{mrna}{$1}{$line}=$_;
    }
}
close F;
my $num=0;
open (O,">$outgff")||die"$!";
for my $k1 (sort keys %{$out{gene2trans}}){
    $num++;
    my $outnum=sprintf("%07d",$num);
    my $outline1="$gff{gene}{$k1}";
    my @outline1=split(/\t/,$outline1);
    $outline1[1]='EVM';
    my $outgene="$species"."G$outnum";
    $outline1[8]="ID=$outgene";
    print O join("\t",@outline1),"\n";
    my $mrnanum=0;
    for my $k2 (sort keys %{$out{gene2trans}{$k1}}){
        $mrnanum++;
        my $outmrna="$species"."T$outnum.$mrnanum";
        for my $k3 (sort{$a<=>$b} keys %{$gff{mrna}{$k2}}){
            my $outline2=$gff{mrna}{$k2}{$k3};
            #$outline2=~s/Name=\S+//;
            #$outline2=~s/evm.model[^;]+/$outmrna/g;
            #$outline2=~s/evm.TU[^;]+/$outgene/;
            my @outline2=split(/\t/,$outline2);
            $outline2[1]="EVM";
            if($outline2[2] eq 'mRNA'){
	$outline2[8]="ID=$outmrna;Parent=$outgene";
            }elsif($outline2[2]=~/UTR|exon/){
	$outline2[8]=~/\.(\w+);Parent=/;
	$outline2[8]="ID=$outmrna.$1;Parent=$outmrna";
            }elsif($outline2[2] eq 'CDS'){
	$outline2[8]="ID=cds.$outmrna;Parent=$outmrna";
            }else{
	die"wrong type $outline2\n";
            }
            print O join("\t",@outline2),"\n";
        }
    }
}
close O;
