#!/usr/bin/perl

use warnings;
use strict;

my $multianno = shift or die $!;
my $overlap = shift or die $!;
my $avinput = shift or die $!;

my %hash;
open IN, "$overlap" or die $!;
while(<IN>){
	chomp;
	my @tmp = split /\t/;
	my $head = "$tmp[0]_$tmp[1]";
	$hash{$head}{'num'} = $tmp[2];
	for (my $i=3; $i<=scalar(@tmp)-6; $i+=6){
		$hash{$head}{'soft'} .= "$tmp[$i],";
	}
	$hash{$head}{'soft'} =~ s/,$//;
#	print "$hash{$head}{'num'}\t$hash{$head}{'soft'}\n";
}
close IN;
open IN, "$avinput" or die $!;
while(<IN>){
	chomp;
	my @tmp = split /\t/;
	my $head = "$tmp[0]_$tmp[1]";
	$hash{$head}{'vaf'} = "$tmp[5]\t$tmp[6]\t$tmp[7]";
}
close IN;
open IN, "$multianno" or die $!;
chomp(my $title = <IN>);
$title =~ s/gff3/genehancer/g;
$title .= "\tRef_dep\tVar_dep\tVaf%\tsoft_num\tsoft_name";
print "$title\n";
while(<IN>){
	chomp;
	s/\t/\t;/g;
	my @tmp = split /\t/;
	$tmp[0] =~ s/^;//;
	$tmp[1] =~ s/^;//;
	my $head = "$tmp[0]_$tmp[1]";
	for (my $i=0; $i<scalar(@tmp); $i++){
		if ($tmp[$i] eq ";"){
			$tmp[$i] = ".";
		}else{
			$tmp[$i] =~ s/^;//;
		}
	}
	$hash{$head}{'annot'} = join("\t",@tmp);
	print "$hash{$head}{'annot'}\t$hash{$head}{'vaf'}\t$hash{$head}{'num'}\t$hash{$head}{'soft'}\n";
}
