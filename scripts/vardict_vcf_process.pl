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
	my @info = split (/;/,$tmp[7]);
	$info[9] =~ /SSF=(.+)/;
	my $pv = $1;
	next if ($pv > 0.05);
	$info[3] =~ /DP=(.+)/;
	my $dp = $1;
	next if ($dp > 200);
	$info[6] =~ /SHIFT3=(.+)/;
	my $shift = $1;
	next if ($shift > 3);
	$info[7] =~ /MSI=(.+)/;
	my $msi = $1;
	next if ($msi  > 3);
	$info[8] =~ /MSILEN=(.+)/;
	my $msilen = $1;
	next if ($msilen  > 1);
	my @N = split (/[:,]/, $tmp[-1]);
	my @T = split (/[:,]/, $tmp[-2]);
	my $nr1 = $N[7];
	my $nr2 = $N[8];
	my $nvaf;
	if ($N[8] == 0){
		$nvaf = 0;
	}else{
		$nvaf = sprintf("%.2f",$N[8] / ($N[7] + $N[8]) * 100);
	}
	my $nr1_p = $N[5];
	my $nr1_m = $N[6];
	my $nr2_p = $N[3];
	my $nr2_m = $N[4];
	my $tr1 = $T[7];
	my $tr2 = $T[8];
	my $tvaf;
	if ($T[8] == 0){
		$tvaf = 0;
	}else{
		$tvaf = sprintf("%.2f",$T[8] / ($T[7] + $T[8]) * 100);
	}
	my $tr1_p = $T[5];
	my $tr1_m = $T[6];
	my $tr2_p = $T[3];
	my $tr2_m = $T[4];
	print "$tmp[0]\t$tmp[1]\t$tmp[3]\t$tmp[4]\t$nr1\t$nr2\t$nvaf\%\t$tr1\t$tr2\t$tvaf\%\t$tr1_p\t$tr1_m\t$tr2_p\t$tr2_m\t$nr1_p\t$nr1_m\t$nr2_p\t$nr2_m\n";

}
close IN;
