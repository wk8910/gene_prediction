#! /usr/bin/env perl
use strict;
use warnings;

my ($config)=@ARGV;

die "Usage: $0 <config file>\n" if(@ARGV<1);

my %config=&readconfig($config);
foreach my $key(sort keys %config){
    my $value=$config{$key};
    print "$key\t$value\n";
}

sub readconfig{
    my $config_file=shift;
    open(SUBI,"< $config_file");
    my %sub_config;
    while (<SUBI>) {
        chomp;
        next if(/^#/);
        next if(/^$/);
        next unless(/^(\S+)\s+=\s+(\S+)/);
        $sub_config{$1}=$2;
    }
    close SUBI;
    return(%sub_config);
}
