#!/usr/bin/perl

=head1 Usage
	perl somatic_process_filter.pl [options] <prefix_somatic.snv|indel.process>
=cut

use warnings;
use strict;
use Cwd qw(abs_path);
use Getopt::Long;
use File::Basename qw(basename);

my ($nvaf, $tvaf, $min_dep, $max_dep, $vaf_dep, $vaf_sub, $help);
GetOptions
(
	 "nvaf:f" => \$nvaf,
	 "tvaf:f" => \$tvaf,
	 "min_dep:i" => \$min_dep,
	 "max_dep:i" => \$max_dep,
	 "vaf_dep:f" => \$vaf_dep,
	 "vaf_sub:f" => \$vaf_sub,
	 "h|help" => \$help
);

die `pod2text $0` if ($help || @ARGV != 1);

my $file = shift;
open IN, "$file" or die $!;

$nvaf ||= 0.05;
$tvaf ||= 0.05;
$min_dep ||= 5;
$max_dep ||= 200;
$vaf_dep ||= 1;
$vaf_sub ||= 0.01;

while(<IN>){
	chomp;
	next if (/^chrom/);
	my @tmp = split /\s+/;
	$tmp[6] =~ s/%//;
	$tmp[9] =~ s/%//;
	my $vaf_plus = 0;
	my $vaf_minus = 0;
	$vaf_plus = $tmp[12]/($tmp[10] + $tmp[12]) if ($tmp[12] > 0);
	$vaf_minus = $tmp[13]/($tmp[11] + $tmp[13]) if ($tmp[13] > 0);
	if (
		($tmp[6]/100) <= $nvaf && 
		($tmp[9]/100) >= $tvaf &&
		($tmp[7] + $tmp[8]) > $min_dep && 
		($tmp[7] + $tmp[8]) < $max_dep &&
		$tmp[8] > $vaf_dep &&
		$vaf_plus > $vaf_sub &&
		$vaf_minus > $vaf_sub 
	){
		print "$_\n";
	}
}
close IN;
