#!/usr/bin/perl

use strict;
use warnings;
use File::Basename qw(fileparse basename dirname);

my $ctx = shift or die $!;
my $prefix = basename($ctx);
$prefix =~ /(\S+)\.ctx$/;
$prefix = $1;

open IN, "$ctx" or die $!;
open OUT, ">$prefix.avinput" or die $!;

while(<IN>){
	chomp;
	next if /^#/;
	my @tmp = split /\s+/;
	$tmp[10] =~ /bam|(\d+)/;
	$tmp[10] = $1;
	my $comments=join (",",@tmp);
	print OUT "$tmp[0]\t$tmp[1]\t$tmp[1]\t0\t0\tcomments:$comments\n";
	print OUT "$tmp[3]\t$tmp[4]\t$tmp[4]\t0\t0\tcomments:$comments\n";
}
close IN;
close OUT;
