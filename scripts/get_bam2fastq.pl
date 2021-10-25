#!/usr/bin/perl

use warnings;
use strict;

my $prefix = shift or die $!;

print <<EOF
#!/bin/bash

#SBATCH -J $prefix\_bam2fq
#SBATCH -p cpu
#SBATCH --mail-type=end
#SBATCH --mail-user=ybliu1988\@163.com
#SBATCH -n 8
#SBATCH --ntasks-per-node=8
#SBATCH -o bam2fastq.slurm.o.%j
#SBATCH -e bam2fastq.slurm.e.%j

SAMTOOLS=/lustre/home/acct-medkkw/medlyb/bin/spack/opt/linux-centos7-x86_64/gcc-5.4.0/samtools-1.6-ruacwclmvjj5ekto3hopfi3ymm5l6yqx/bin/samtools
BEDTOOLS=/lustre/home/acct-medkkw/medlyb/bin/spack/opt/linux-centos7-x86_64/gcc-5.4.0/bedtools2-2.27.1-bduy2d7yx3cl2uukzsnlool7pohko7v4/bin/bedtools

uname -n &&\\
echo "Mission starts at `date`" &&\\
	 \$SAMTOOLS sort -n -@ 16 $prefix.bam -T $prefix -o $prefix.bam.qsort &&\\
	 \$BEDTOOLS bamtofastq -i $prefix.bam.qsort -fq $prefix\_1.fq -fq2 $prefix\_2.fq &&\\
	rm $prefix.bam.qsort &&\\
	gzip $prefix\_1.fq $prefix\_2.fq &&\\
echo "Mission completes at `date`"

EOF
