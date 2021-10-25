#!/usr/bin/perl

use strict;
use warnings;

my $file = shift or die $!;

my %hash;
open IN, "$file" or die $!;
my $head = <IN>;
print $head;
while(<IN>){
	chomp;
	my @tmp = split /\t/;
	$hash{$tmp[4]}{$tmp[5]}{$tmp[12]}{$tmp[14]} = $_;
}
close IN;

my @chr;
open IN, "/lustre/home/acct-medkkw/medlyb/chr.list" or die $!;
while(<IN>){
	chomp;
	push @chr, $_;
}
close IN;

for my $i (@chr){
	for my $j (sort {$a <=> $b} keys %{$hash{$i}}){
		for my $k (keys %{$hash{$i}{$j}}){
			for my $h (keys %{$hash{$i}{$j}{$k}}){
				 print "$hash{$i}{$j}{$k}{$h}\n";
			}
		}
	}
}

close IN;
