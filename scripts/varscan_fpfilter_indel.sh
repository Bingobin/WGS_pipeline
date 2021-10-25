PREFIX=${1}_D0-MN
DIR=/lustre/home/acct-medkkw/medlyb/project/$1/WGS/02.somatic_D0/varscan
BAM=/lustre/home/acct-medkkw/medlyb/project/$1/WGS/01.align/${1}_D0/${1}_D0.sorted.merged.reorder.rmdup.realign.recal.bam

awk 'BEGIN{getline}{if($4 ~ /^-/){i=$2+length($4);print $1"\t"$2"\t"i}else{j=$2+1;print $1"\t"$2"\t"j}}' $DIR/${PREFIX}.varscan.indel.Somatic.hc > $DIR/${PREFIX}.varscan.indel.Somatic.hc.pos
bam-readcount -b 20 -f /lustre/home/acct-medkkw/medlyb/database/annotation/gatk_ann/hg38/hg38bundle/Homo_sapiens_assembly38.fasta -l $DIR/${PREFIX}.varscan.indel.Somatic.hc.pos $BAM > $DIR/${PREFIX}.varscan.indel.Somatic.hc.rc
java -Djava.io.tmpdir=/lustre/home/acct-medkkw/medlyb/tmp -jar  /lustre/home/acct-medkkw/medlyb/bin/VarScan/VarScan.v2.4.3.jar fpfilter $DIR/${PREFIX}.varscan.indel.Somatic.hc $DIR/${PREFIX}.varscan.indel.Somatic.hc.rc --output-file $DIR/${PREFIX}.varscan.indel.Somatic.hc.pass
cut -f 1-7,9-11,16-23 $DIR/${PREFIX}.varscan.indel.Somatic.hc.pass > $DIR/${PREFIX}.varscan.indel.process.tmp
perl /lustre/home/acct-medkkw/medlyb/script/wgs_pipeline/somatic_process_filter.pl $DIR/${PREFIX}.varscan.indel.process.tmp > $DIR/${PREFIX}.varscan.indel.process
