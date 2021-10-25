#!/bin/bash

PRE=${1}_D0
PREFIX=${1}_D0-MN
DIR=../../02.somatic
D0_BAM=../../01.align/$PRE/${PRE}.sorted.merged.reorder.rmdup.realign.recal.bam

BIN=/lustre/home/acct-medkkw/medlyb/script/wgs_pipeline
GENOME=/lustre/home/acct-medkkw/medlyb/database/annotation/gatk_ann/hg38/bwaindex2/Homo_sapiens_assembly38.fasta

perl $BIN/somatic_process_snv_overlap.pl  -mutect2 $DIR/mutect2/${PREFIX}.mutect2.snv.process -vardict $DIR/vardict/${PREFIX}.vardict.snv.process -varscan $DIR/varscan/${PREFIX}.varscan.snv.process -somaticsniper $DIR/somaticsniper/${PREFIX}.somaticsniper.snv.process -muse $DIR/muse/${PREFIX}.muse.snv.process -strelka $DIR/manta_strelka/${PREFIX}.strelka.snv.process  > ${PREFIX}.combind.snv.overlap &&\
awk '{OFS="\t";print $1,$2,$2}' ${PREFIX}.combind.snv.overlap  > ${PREFIX}.combind.snv.overlap.pos &&\
bam-readcount -w 10 -q 20 -b 20 -i -f $GENOME -l ${PREFIX}.combind.snv.overlap.pos $D0_BAM > ${PREFIX}.combind.snv.overlap.pos.rc  &&\
perl $BIN/calculate_vaf_by_rc.pl ${PREFIX}.combind.snv.overlap.pos.rc ${PREFIX}.combind.snv.overlap > ${PREFIX}.combind.snv.overlap.avinput 
perl /lustre/home/acct-medkkw/medlyb/soft/annovar/table_annovar.pl ${PREFIX}.combind.snv.overlap.avinput /lustre/home/acct-medkkw/medlyb/soft/annovar/humandb/ -buildver hg38 -protocol refGene,rmsk,gff3,cosmic87,avsnp150 --gff3dbfile hg38_genehancer.txt   -operation g,r,r,f,f --remove  --outfile ${PREFIX}.combind.snv.overlap  &&\
perl $BIN/annovar_add_vaf_soft.pl ${PREFIX}.combind.snv.overlap.hg38_multianno.txt ${PREFIX}.combind.snv.overlap ${PREFIX}.combind.snv.overlap.avinput   > ${PREFIX}.combind.snv.overlap.anno &&\
awk '{if($(NF-1) >= 2){print $0}}'  ${PREFIX}.combind.snv.overlap.anno  > ${PREFIX}.combind.snv.overlap.anno.m2 &&\
perl $BIN/annovar2maf.pl ${PREFIX}.combind.snv.overlap.anno $PRE >  ${PREFIX}.combind.snv.overlap.anno.maf
#rm nohup.out
