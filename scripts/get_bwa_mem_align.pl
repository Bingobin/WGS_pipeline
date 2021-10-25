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
#SBATCH -p cpu
#SBATCH --mail-type=end
#SBATCH --mail-user=ybliu1988\@163.com
#SBATCH -n 40
#SBATCH --exclusive
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

##########################Step01: bwa align################################################
echo "Mission bwa starts at `date`" &&\\

for i in \${LANE[@]};
do
file1=\$DATADIR/\${DATAPRE}_\${i}_R1.fastq.gz
file2=\$DATADIR/\${DATAPRE}_\${i}_R2.fastq.gz
echo "BWA \${PREFIX} \$i starts at `date`" &&\ bwa mem -t 40 -M -R "\@RG\\tID:\$i\\tLB:\$LAB\\tSM:\$SM\\tPL:ILLUMINA" -Y \$REF \$file1 \$file2 | samtools view -Sb -o \${PREFIX}_\${i}.bam &&\ echo "BWA \${PREFIX} \$i DONE at `date`"
MERGE_FILE="\$MERGE_FILE I=\${PREFIX}_\${i}.bam"
done

echo "Mission bwa completes at `date`" &&\\

##########################Step02: Picard Merge############################################

#echo "Mission Sort and  Merge completes starts at `date`" &&\\
#\$JAVA -Xmx60g -jar \$PICARD MergeSamFiles SO=coordinate \$MERGE_FILE   O=\${PREFIX}.sorted.merged.bam TMP_DIR=\$TMP &&\\
#echo "Mission  sort completes at `date`" &&\\

##########################Step03: Picard Reorder##########################################

#echo "Mission reorder at `date`" &&\\
#\$JAVA -Xmx60g -jar \$PICARD ReorderSam  I=\${PREFIX}.sorted.merged.bam  O=\${PREFIX}.sorted.merged.reorder.bam R=\$REF VALIDATION_STRINGENCY=LENIENT TMP_DIR=\$TMP && rm \${PREFIX}.sorted.merged.bam &&\\
#samtools stats \${PREFIX}.sorted.merged.reorder.bam >  \${PREFIX}.sorted.merged.reorder.bam.stats &&\\
#echo "Mission reorder completes at `date`" &&\\

##########################Step04: Picard rm duplication###################################

#echo "Mission markdup starts at `date`" &&\\
#\$JAVA -Xmx60g -jar \$PICARD MarkDuplicates I=\${PREFIX}.sorted.merged.reorder.bam O=\${PREFIX}.sorted.merged.reorder.rmdup.bam METRICS_FILE=\${PREFIX}.sorted.merged.reorder.rmdup.bam.mat REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true TMP_DIR=\$TMP && rm \${PREFIX}.sorted.merged.reorder.bam  &&\\
#echo "Mission markdup completes at `date`"&&\\

#########################Step05: GATK BQSR###############################################

#echo "Mission realnRecal starts at `date`" &&\\
#\$JAVA -Djava.io.tmpdir=\$TMP -Xmx60g -jar \$GATK -T RealignerTargetCreator -R \$REF -I \${PREFIX}.sorted.merged.reorder.rmdup.bam  -known \$GATK_ann/Homo_sapiens_assembly38.known_indels.vcf -known \$GATK_ann/Mills_and_1000G_gold_standard.indels.hg38.vcf -o \${PREFIX}.sorted.merged.reorder.rmdup.intervals  &&\\
#echo "RT done at `date`" &&\\

#\$JAVA -Djava.io.tmpdir=\$TMP -Xmx25g -jar \$GATK -T IndelRealigner  -R \$REF -I \${PREFIX}.sorted.merged.reorder.rmdup.bam -known \$GATK_ann/Homo_sapiens_assembly38.known_indels.vcf -known \$GATK_ann/Mills_and_1000G_gold_standard.indels.hg38.vcf -targetIntervals \${PREFIX}.sorted.merged.reorder.rmdup.intervals -o \${PREFIX}.sorted.merged.reorder.rmdup.realign.bam  && rm \${PREFIX}.sorted.merged.reorder.rmdup.bam &&\\
#echo "IR done at `date`" &&\\

#\$JAVA -Djava.io.tmpdir=\$TMP -Xmx60g -jar \$GATK -T BaseRecalibrator -R \$REF -I \${PREFIX}.sorted.merged.reorder.rmdup.realign.bam -knownSites \$GATK_ann/dbsnp_146.hg38.vcf  -knownSites \$GATK_ann/Mills_and_1000G_gold_standard.indels.hg38.vcf  -knownSites \$GATK_ann/Homo_sapiens_assembly38.known_indels.vcf   -o \${PREFIX}.sorted.merged.reorder.rmdup.realign.recal.table -nct 40  &&\\
#echo "BR done at `date`" &&\\

#\$JAVA -Djava.io.tmpdir=\$TMP -Xmx25g -jar \$GATK -T PrintReads -R \$REF -I \${PREFIX}.sorted.merged.reorder.rmdup.realign.bam -BQSR \${PREFIX}.sorted.merged.reorder.rmdup.realign.recal.table -o \${PREFIX}.sorted.merged.reorder.rmdup.realign.recal.bam -nct 40 && rm \${PREFIX}.sorted.merged.reorder.rmdup.realign.bam &&\\
#echo "PR done at `date`" &&\\
#samtools stats \${PREFIX}.sorted.merged.reorder.rmdup.realign.recal.bam >  \${PREFIX}.sorted.merged.reorder.rmdup.realign.recal.bam.stats &&\\
#echo "stats done at `date`" &&\\
#echo "Mission realnRecal completes at `date`" &&\\

########################Step06: Other####################################################

#echo "Prepared for verscan ... at `date`" &&\\
#samtools mpileup -q 20 -f \$REF \${PREFIX}.sorted.merged.reorder.rmdup.realign.recal.bam > \${PREFIX}.mpileup &&\\
#		gzip \${PREFIX}.mpileup &&\\
#echo "Prepared for copycat ... at `date`" &&\\
#bam-window -w 10000 -l -r -i \${PREFIX}.sorted.merged.reorder.rmdup.realign.recal.bam -o \${PREFIX}.10k.windfile &&\\
#/lustre/home/acct-medkkw/medlyb/soft/samtools-0.1.16/samtools pileup -cvf \$REF \${PREFIX}.sorted.merged.reorder.rmdup.realign.recal.bam > \${PREFIX}.10col.mpileup &&\\
#echo "Prepared for lumpy ... at `date`" &&\\
#samtools view -b -F 1294 \${PREFIX}.sorted.merged.reorder.rmdup.realign.recal.bam > \${PREFIX}.discordants.bam 
#samtools view -h \${PREFIX}.sorted.merged.reorder.rmdup.realign.recal.bam | /lustre/home/acct-medkkw/medlyb/bin/lumpy-sv/scripts/extractSplitReads_BwaMem -i stdin | samtools view -Sb - > \${PREFIX}.splitters.bam &&\\

echo "All Mission Done!"
#########################ALL DONE#########################################################
EOF
