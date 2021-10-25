PRE=${1}_D0
PREFIX=${1}_D0-MN
DIR=../../02.somatic
BIN=/lustre/home/acct-medkkw/medlyb/script/wgs_pipeline


perl $BIN/somatic_process_indel_overlap.pl  -mutect2  $DIR/mutect2/${PREFIX}.mutect2.indel.process -vardict $DIR/vardict/${PREFIX}.vardict.indel.process -varscan $DIR/varscan/${PREFIX}.varscan.indel.process -strelka $DIR/manta_strelka/${PREFIX}.strelka.indel.process  > ${PREFIX}.combind.indel.overlap
perl $BIN/calculate_avinput_by_indel.pl  ${PREFIX}.combind.indel.overlap  > ${PREFIX}.combind.indel.overlap.avinput ##remove indel dep <5 
perl /lustre/home/acct-medkkw/medlyb/soft/annovar/table_annovar.pl ${PREFIX}.combind.indel.overlap.avinput /lustre/home/acct-medkkw/medlyb/soft/annovar/humandb/ -buildver hg38 -protocol refGene,rmsk,gff3,cosmic87,avsnp150 --gff3dbfile hg38_genehancer.txt   -operation g,r,r,f,f --remove  --outfile ${PREFIX}.combind.indel.overlap
perl $BIN/annovar_indel_vaf_soft.pl ${PREFIX}.combind.indel.overlap.hg38_multianno.txt ${PREFIX}.combind.indel.overlap.avinput   > ${PREFIX}.combind.indel.overlap.anno
awk '{if($(NF-1) >= 2){print $0}}'  ${PREFIX}.combind.indel.overlap.anno  > ${PREFIX}.combind.indel.overlap.anno.m2
perl $BIN/annovar2maf.pl ${PREFIX}.combind.indel.overlap.anno $PRE > ${PREFIX}.combind.indel.overlap.anno.maf
#bedtools intersect  -a ${PREFIX}.combind.indel.overlap.anno -b ~/script/wgs_pipeline/GeneHancer_version_4_4.txt  -wa -wb >  ${PREFIX}.combind.indel.genehancer.txt
