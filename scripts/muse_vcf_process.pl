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
	next if ($tmp[6] ne "PASS");
	next if ($tmp[1] =~ /_/);
	my @T = split (/:/, $tmp[9]);
	my @N = split (/:/, $tmp[10]);
	next if ($T[0] =~ /2/ || $N[0] =~ /2/);
	my ($t_r1, $t_r2) = split(/,/, $T[2]);
	my $t_vaf = sprintf("%.2f",$t_r2/($t_r1+$t_r2)*100);
	my ($n_r1, $n_r2) = split(/,/, $N[2]);
	my $n_vaf = sprintf("%.2f",$n_r2/($n_r1+$n_r2)*100);
	print "$tmp[0]\t$tmp[1]\t$tmp[3]\t$tmp[4]\t$n_r1\t$n_r2\t$n_vaf\t$t_r1\t$t_r2\t$t_vaf\t5\t5\t5\t5\t5\t5\t5\t5\n";
}
close IN;
