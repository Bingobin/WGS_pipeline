#!/usr/bin/perl

use warnings;
use strict;

print "chrom\tposition\tref\tvar\tnormal_reads1\tnormal_reads2\tnormal_var_freq\ttumor_reads1\ttumor_reads2\ttumor_var_freq\ttumor_reads1_plus\ttumor_reads1_minus\ttumor_reads2_plus\ttumor_reads2_minus\tnormal_reads1_plus\tnormal_reads1_minus\tnormal_reads2_plus\tnormal_reads2_minus\n";

my %hash;
$hash{'A'} = 1+3;
$hash{'C'} = 2+3;
$hash{'G'} = 3+3;
$hash{'T'} = 4+3;

my $file = shift or die $!;
open IN, "$file" or die $!;
while(<IN>){
	chomp;
	next if (/^#/);
	my ($t_r1, $t_r2, $t_vaf, $n_r1, $n_r2, $n_vaf);
	my @tmp = split /\s+/;
	my @T = split (/:/, $tmp[10]);
	my @N = split (/:/, $tmp[9]);
	if(length($tmp[3])==1 && length($tmp[4])==1){
		$t_r1 = (split(/,/, $T[$hash{$tmp[3]}]))[0];
		$t_r2 = (split(/,/, $T[$hash{$tmp[4]}]))[0];
		next if (($t_r1+$t_r2) == 0);
		$t_vaf = sprintf("%.2f",$t_r2/($t_r1+$t_r2)*100);
		$n_r1 = (split(/,/, $N[$hash{$tmp[3]}]))[0];
		$n_r2 = (split(/,/, $N[$hash{$tmp[4]}]))[0];
		next if (($n_r1+$n_r2) == 0);
		$n_vaf = sprintf("%.2f",$n_r2/($n_r1+$n_r2)*100);
	}else{
		$t_r1 = (split(/,/, $T[2]))[0];
		$t_r2 = (split(/,/, $T[3]))[0];
		next if (($t_r1+$t_r2) == 0);
		$t_vaf = sprintf("%.2f",$t_r2/($t_r1+$t_r2)*100);
		$n_r1 = (split(/,/, $N[2]))[0];
		$n_r2 = (split(/,/, $N[3]))[0];
		next if (($n_r1+$n_r2) == 0);
		$n_vaf = sprintf("%.2f",$n_r2/($n_r1+$n_r2)*100);
	}
	print "$tmp[0]\t$tmp[1]\t$tmp[3]\t$tmp[4]\t$n_r1\t$n_r2\t$n_vaf\t$t_r1\t$t_r2\t$t_vaf\t5\t5\t5\t5\t5\t5\t5\t5\n";
}
close IN;
