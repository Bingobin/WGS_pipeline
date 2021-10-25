#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my ($mutect2, $strelka, $varscan, $vardict, $help);
GetOptions(
		"mutect2:s" => \$mutect2,
		"strelka:s" => \$strelka,
		"varscan:s" => \$varscan,
		"vardict:s" => \$vardict,
		"h"=>\$help,
		
);

my $usage=<<USAGE;
Discription: Combind the different resut from five snv caller(mutect2, strelka, vardict, varscan, somaticsniper & muse);
Author:Liuyabin
Data: 2018-10
usage:perl $0 -mutect2 <snv.process> -strelka <snv.process>  -varscan <snv.process> -vardict <snv.process>
USAGE

die $usage if (!$mutect2 || !$vardict || !$varscan || !$strelka || $help);

#my $mutect2 = shift or die $!;
#my $vardict = shift or die $!;
#my $varscan = shift or die $!;
#my $somaticsniper = shift or die $!;
#my $muse = shift or die $!;

my %hash;

read_input($mutect2, "Mutect2", \%hash);
read_input($strelka, "Strelka", \%hash);
read_input($varscan, "VarScan", \%hash);
read_input($vardict, "VarDict", \%hash);

for my $i (sort keys %hash){
	for my $j (sort {$a<=>$b} (keys %{$hash{$i}})){
		print "$i\t$j\t$hash{$i}{$j}{'num'}\t$hash{$i}{$j}{'lab'}\n";
	}
}

sub read_input{
	my $file = shift;
	my $label = shift;
	my $h = shift;
	open IN, "$file" or die $!;
	while(<IN>){
		chomp;
		next if (/^chrom/);
		my @tmp = split /\s+/;
		$h->{$tmp[0]}->{$tmp[1]}{'lab'} .= "$label\t$tmp[2]\t$tmp[3]\t$tmp[7]\t$tmp[8]\t$tmp[9]\t";
		$h->{$tmp[0]}->{$tmp[1]}{'num'} += 1;
	}
	close IN;
}
