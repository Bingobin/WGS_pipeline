#!/usr/bin/perl

use strict;
use warnings;

my $ol = shift or die $!;

##Chr     Start   End     Ref     Alt  Ref_dep Var_dep Vaf%    soft_num        soft_name

my %hash;
open IN, "$ol" or die $!;
while(<IN>){
	chomp;
	my @tmp = split /\t/;
#	next if ($tmp[0] eq "chrX");
#	next if ($tmp[0] eq "chrY");
#	next if ($tmp[0] eq "chrM");
	my $ref_dep = $tmp[6];
	my $alt_dep = $tmp[7];
	next if (($ref_dep + $alt_dep) < 5);
	my $vaf = sprintf("%.4f", $alt_dep / ($ref_dep + $alt_dep) * 100);
#	next if ($vaf > 80 || $vaf < 1);
	my $chr = $tmp[0];
	my $soft_num = $tmp[2];
	my $soft_name;
	for (my $i=3; $i<=scalar(@tmp)-6; $i+=6){
		$soft_name .= "$tmp[$i],";
	}
	$soft_name  =~ s/,$//;
	my ($start, $end, $ref, $alt);
	if ($tmp[5] =~ /^-/){
		$start = $tmp[1] + 1;
		$end = $tmp[1] + length($tmp[5]) - 1;
		$tmp[5] =~ s/^-//;
		$ref = $tmp[5];
		$alt = "-";
	}elsif($tmp[5] =~ /^\+/){
		$start = $tmp[1];
		$end = $tmp[1];
		$tmp[5] =~ s/^\+//;
		$ref = "-";
		$alt = $tmp[5];
	}else{
		if(length($tmp[4]) == 1){
			$start = $tmp[1];
			$end = $tmp[1];
			$ref = "-";
			$alt = substr($tmp[5],1);
		}elsif(length($tmp[5]) == 1){
			$start = $tmp[1] + 1;
			$end = $tmp[1] + length($tmp[4]) - 1;
			$ref = substr($tmp[4],1);
			$alt = "-";
		}else{
			$start = $tmp[1] + 1;
			$end = $tmp[1] + length($tmp[4]) - 1;
			$ref = substr($tmp[4],1);
			$alt = substr($tmp[5],1);
		}
	}
	print "$chr\t$start\t$end\t$ref\t$alt\t$ref_dep\t$alt_dep\t$vaf\t$soft_num\t$soft_name\n";
}
close IN;
