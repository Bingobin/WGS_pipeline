#!/usr/bin/perl

use warnings;
use strict;

print "chrom\tposition\tref\tvar\tnormal_reads1\tnormal_reads2\tnormal_var_freq\ttumor_reads1\ttumor_reads2\ttumor_var_freq\ttumor_reads1_plus\ttumor_reads1_minus\ttumor_reads2_plus\ttumor_reads2_minus\tnormal_reads1_plus\tnormal_reads1_minus\tnormal_reads2_plus\tnormal_reads2_minus\n";

my $file = shift or die $!;
open IN, "$file" or die $!;
while(<IN>){
	chomp;
	next if (/^#/);
	my @tmp = split /\s+/;
	next if ($tmp[1] =~ /_/);
	my @T = split (/:/, $tmp[10]);
	my @N = split (/:/, $tmp[9]);
	next if ($T[0] =~ /2/ || $N[0] =~ /2/);
	my ($t_r1_p, $t_r1_m, $t_r2_p, $t_r2_m) = split(/,/, $T[3]);
	my $t_r1 = $t_r1_p + $t_r1_m;
	my $t_r2 = $t_r2_p + $t_r2_m;
	my $t_vaf = sprintf("%.2f",$t_r2/($t_r1+$t_r2)*100);
	my ($n_r1_p, $n_r1_m, $n_r2_p, $n_r2_m) = split(/,/, $N[3]);
	my $n_r1 = $n_r1_p + $n_r1_m;
	my $n_r2 = $n_r2_p + $n_r2_m;
	my $n_vaf = sprintf("%.2f",$n_r2/($n_r1+$n_r2)*100);
	print "$tmp[0]\t$tmp[1]\t$tmp[3]\t$tmp[4]\t$n_r1\t$n_r2\t$n_vaf\t$t_r1\t$t_r2\t$t_vaf\t$t_r1_p\t$t_r1_m\t$t_r2_p\t$t_r2_m\t$n_r1_p\t$n_r1_m\t$n_r2_p\t$n_r2_m\n";
}
close IN;
