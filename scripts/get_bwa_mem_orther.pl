#!/usr/bin/perl

use warnings;
use strict;

my $namepre = shift or die "perl $0 XHF_MN CHG028224 \"L001 L002 L003 combined\" ../../00.data: $!";
my $datapre = shift or die "perl $0 XHF_MN CHG028224 \"L001 L002 L003 combined\" ../../00.data: $!";
my $lane = shift or die "perl $0 XHF_MN CHG028224 \"L001 L002 L003 combined\" ../../00.data: $!";
my $datadir = shift or die "perl $0 XHF_MN CHG028224 \"L001 L002 L003 combined\" ../../00.data: $!";


my @lane = split (/\s+/, $lane);
$lane = "";
for my $i (@lane){
	$lane .= " \"$i\" ";
}


print <<EOF
#!/bin/bash

#SBATCH -J $namepre\_align
#SBATCH -p small
#SBATCH -n 2
#SBATCH -o bwa_mem.slurm.o.%j
#SBATCH -e bwa_mem.slurm.e.%j

uname -n &&\\

TMP=/lustre/home/acct-medkkw/medlyb/tmp
REF=/lustre/home/acct-medkkw/medlyb/database/annotation/gatk_ann/hg38/bwaindex2/Homo_sapiens_assembly38.fasta
PICARD=/lustre/home/acct-medkkw/medlyb/bin/picard.jar
GATK=/lustre/home/acct-medkkw/medlyb/bin/GenomeAnalysisTK.jar
GATK_ann=/lustre/home/acct-medkkw/medlyb/database/annotation/gatk_ann/hg38/hg38bundle
JAVA=/usr/bin/java


##########################################################################################
DATADIR=$datadir
DATAPRE=$datapre
PREFIX=$namepre
SM=\$PREFIX
LAB=\${SM}_LIB1
LANE=( $lane )
MERGE_FILE=""

########################Step06: Other####################################################

echo "Prepared for verscan ... at `date`" &&\\
samtools mpileup -q 20 -f \$REF \${PREFIX}.sorted.merged.reorder.rmdup.realign.recal.bam > \${PREFIX}.mpileup &&\\
#		gzip \${PREFIX}.mpileup &&\\
#echo "Prepared for copycat ... at `date`" &&\\
#bam-window -w 10000 -l -r -i \${PREFIX}.sorted.merged.reorder.rmdup.realign.recal.bam -o \${PREFIX}.10k.windfile &&\\
#/lustre/home/acct-medkkw/medlyb/soft/samtools-0.1.16/samtools pileup -cvf \$REF \${PREFIX}.sorted.merged.reorder.rmdup.realign.recal.bam > \${PREFIX}.10col.mpileup &&\\
#echo "Prepared for lumpy ... at `date`" &&\\
#samtools view -b -F 1294 \${PREFIX}.sorted.merged.reorder.rmdup.realign.recal.bam > \${PREFIX}.discordants.bam 
#samtools view -h \${PREFIX}.sorted.merged.reorder.rmdup.realign.recal.bam | /lustre/home/acct-medkkw/medlyb/bin/lumpy-sv/scripts/extractSplitReads_BwaMem -i stdin | samtools view -Sb - > \${PREFIX}.splitters.bam &&\\

echo "All Mission Done at `date`"
#########################ALL DONE#########################################################
EOF
