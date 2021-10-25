#!/usr/bin/perl

use strict;

use warnings;

print "chrom\tposition\tref\tvar\tnormal_reads1\tnormal_reads2\tnormal_var_freq\ttumor_reads1\ttumor_reads2\ttumor_var_freq\ttumor_reads1_plus\ttumor_reads1_minus\ttumor_reads2_plus\ttumor_reads2_minus\tnormal_reads1_plus\tnormal_reads1_minus\tnormal_reads2_plus\tnormal_reads2_minus\n";

my $file = shift or die $!;
open IN, "$file" or die $!;
while(<IN>){
	chomp;
	next if (/^#/);
	my @tmp = split /\s+/;
	my @N = split (/:/, $tmp[-1]);
	my @T = split (/:/, $tmp[-2]);
	$N[1] =~ s/,/\t/;
	$N[2] = $N[2] * 100;
	$T[1] =~ s/,/\t/;
	$T[2] = $T[2] * 100;
	print "$tmp[0]\t$tmp[1]\t$tmp[3]\t$tmp[4]\t$N[1]\t$N[2]\%\t$T[1]\t$T[2]\%\t$T[-2]\t$T[-1]\t$T[3]\t$T[4]\t$N[-2]\t$N[-1]\t$N[3]\t$N[4]\n";

}
close IN;
