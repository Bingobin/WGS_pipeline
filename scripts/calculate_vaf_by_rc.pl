#!/usr/bin/perl

use strict;
use warnings;

my $rc = shift or die $!;
my $ol = shift or die $!;

my %hash;
open IN, "$ol" or die $!;
while(<IN>){
	chomp;
	my @tmp = split /\t/;
	my $header = "$tmp[0]\t$tmp[1]\t$tmp[4]";
	$hash{$header} = $tmp[5];
}
close IN;

my %chr;
open IN, "/lustre/home/acct-medkkw/medlyb/chr.list" or die $!;
while(<IN>){
	chomp;
	$chr{$_} = 1;
}
close IN;

open IN, "$rc" or die $!;
#print "CHR\tPOS\tREF\tDEP\tA\tC\tG\tT\n";
while(<IN>){
	chomp;
	my @tmp = split /\t/;
	next unless (defined $chr{$tmp[0]});
#	next if (scalar(@tmp) > 10);
#	next if ($tmp[0] eq "chrX");
#	next if ($tmp[0] eq "chrY");
#	next if ($tmp[0] eq "chrM");
	my @A = split (/:/, $tmp[5]);
	my @C = split (/:/, $tmp[6]);
	my @G = split (/:/, $tmp[7]);
	my @T = split (/:/, $tmp[8]);
#	print "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$A[1]\t$C[1]\t$G[1]\t$T[1]\n";
	next if ($tmp[3] <5);
	my $header = "$tmp[0]\t$tmp[1]\t$tmp[2]";
	&cal_vaf($tmp[0], $tmp[1], $tmp[3],\@A, $tmp[2], $hash{$header}) if ($tmp[2] eq "A");
	&cal_vaf($tmp[0], $tmp[1], $tmp[3],\@C, $tmp[2], $hash{$header}) if ($tmp[2] eq "C");
	&cal_vaf($tmp[0], $tmp[1], $tmp[3],\@G, $tmp[2], $hash{$header}) if ($tmp[2] eq "G");
	&cal_vaf($tmp[0], $tmp[1], $tmp[3],\@T, $tmp[2], $hash{$header}) if ($tmp[2] eq "T");


}
close IN;
sub cal_vaf{
	my $chr = shift;
	my $pos = shift;
	my $dep = shift;
	my $base = shift;
	my $ref = shift;
	my $var = shift;
	my $r1 = $base->[1];
	my $r2 = $dep - $r1;
	my $vaf = sprintf("%.4f",$r2 / $dep * 100); 
	print "$chr\t$pos\t$pos\t$ref\t$var\t$r1\t$r2\t$vaf\n";
	#if ($vaf > 1 && $vaf < 80 );
}
