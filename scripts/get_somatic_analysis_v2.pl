#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;

my($gvcf, $mutect2, $vardict, $varscan, $somaticsniper, $muse, $manta_strelka, $copycat, $cnvkit, $breakdancer, $lumpy,  $cancer, $normal, $datadir, $outdir, $prefix, $all, $help);

GetOptions(
 	 "gvcf" => \$gvcf,
	 "mutect2" => \$mutect2,
	 "vardict" => \$vardict,
	 "varscan" => \$varscan,
	 "somaticsniper" => \$somaticsniper,
	 "muse" => \$muse,
	 "manta_strelka" => \$manta_strelka,
	 "copycat" => \$copycat,
	 "cnvkit" => \$cnvkit,
	 "breakdancer" => \$breakdancer,
	 "lumpy" => \$lumpy,
	 "c:s" => \$cancer,
	 "n:s" => \$normal,
	 "d:s" => \$datadir,
	 "o:s" => \$outdir,
	 "p:s" => \$prefix,
	 "a" => \$all,
	 "h" => \$help
		
);

my $usage=<<USAGE;
Usage:perl $0 -c XHF_D0 -n XHF_MN -p XHF_D0-MN -o /lustre/home/acct-medkkw/medlyb/project/LZY/WGS/02.somatic_D0 -d /lustre/home/acct-medkkw/medlyb/project/XHF/WGS/01.align -a
Availble analysis: gvcf mutect2 vardict somaticsniper muse manta_strelka copycat cnvkit  breakdancer lumpy 

USAGE

die $usage if ($help || !$cancer ||  !$normal || !$datadir || !$prefix || !$outdir);

my @chr;
my $BIN = "/lustre/home/acct-medkkw/medlyb/script/wgs_pipeline";

open IN, "/lustre/home/acct-medkkw/medlyb/chr.list" or die $!;
while(<IN>){
	chomp;
	push @chr, $_;
}

close IN;

my $GENOME="/lustre/home/acct-medkkw/medlyb/database/annotation/gatk_ann/hg38/bwaindex2/Homo_sapiens_assembly38.fasta";
my $VarScan="/lustre/home/acct-medkkw/medlyb/bin/VarScan/VarScan.v2.4.3.jar";
my $VARDICT_HOME="/lustre/home/acct-medkkw/medlyb/bin/VarDict-1.5.1/bin";
my $GATK="/lustre/home/acct-medkkw/medlyb/bin/GenomeAnalysisTK.jar";
my $COSMIC="/lustre/home/acct-medkkw/medlyb/database/annotation/mutect_ann/hg38_v79/Cosmic.v79.hg38.sorted.vcf";
my $DBSNP="/lustre/home/acct-medkkw/medlyb/database/annotation/gatk_ann/hg38/hg38bundle/dbsnp_146.hg38.vcf";
my $normal_bam="$datadir/$normal/$normal.sorted.merged.reorder.rmdup.realign.recal.bam";
my $cancer_bam="$datadir/$cancer/$cancer.sorted.merged.reorder.rmdup.realign.recal.bam";
my $TMP="/lustre/home/acct-medkkw/medlyb/tmp";


if ($mutect2 || $all){
	my $soft = "mutect2";
	my $dir = "$outdir/$soft";
	my $outpre = "$prefix.$soft";
	if (-d $dir){
		die "ERROR:The dietectory $dir has exist!";
	}else{
		`mkdir $dir`;
	}
#	die "ERROR:The dietectory $dir has exist!" ? (-d $dir): `mkdir $dir`;
	for my $i (@chr){
		open OUT, ">$dir/$outpre.$i.slurm" or die $!;
		print OUT<<EOF;
#!/bin/bash

#SBATCH -J $outpre.$i
#SBATCH -p cpu
#SBATCH --mail-type=end
#SBATCH --mail-user=ybliu1988\@163.com
#SBATCH -n 2
#SBATCH --ntasks-per-node=2
#SBATCH -o $outpre.$i.slurm.o.%j
#SBATCH -e $outpre.$i.slurm.e.%j
echo "Mission starts at `date`" &&\
uname -n &&\

gatk Mutect2 -R $GENOME -I $cancer_bam -tumor $cancer  -I $normal_bam -normal $normal  --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter -O $outpre.$i.vcf -L $i  &&\\
gatk FilterMutectCalls -V $outpre.$i.vcf -O $outpre.$i.f.vcf &&\\

echo "Mission completes at `date`"
EOF
		close OUT;
	}
	open OUT, ">$dir/$outpre.process.sh" or die $!;
	print OUT<<EOF;
for i in `cat ~/chr.list`; do awk '{if(\$7 == "PASS"){print \$0}}'  $outpre.\$i.f.vcf;done >> $outpre.PASS.vcf&&\\
awk '{i=length(\$4);j=length(\$5);if(i==1 && j ==1 ){print \$0}}'  $outpre.PASS.vcf > $outpre.snv.vcf &&\\
awk '{i=length(\$4);j=length(\$5);if(i==1 && j ==1 ){}else{print \$0}}'  $outpre.PASS.vcf > $outpre.indel.vcf &&\\
perl $BIN/mutect2_vcf_process.pl $outpre.indel.vcf > $outpre.indel.process &&\\
perl $BIN/mutect2_vcf_process.pl $outpre.snv.vcf > $outpre.snv.process
EOF
	close OUT;
}

if ($varscan || $all){
	my $soft = "varscan";
	my $dir = "$outdir/$soft";
	my $outpre = "$prefix.$soft";
	if (-d $dir){
		die "ERROR:The dietectory $dir has exist!";
	}else{
		`mkdir $dir`;
	}
	open OUT, ">$dir/$outpre.slurm" or die $!;
	print OUT<<EOF;
#!/bin/bash
#SBATCH -J $outpre
#SBATCH -p cpu
#SBATCH --mail-type=end
#SBATCH --mail-user=ybliu1988\@163.com
#SBATCH -n 2
#SBATCH --ntasks-per-node=2
#SBATCH -o $outpre.slurm.o.%j
#SBATCH -e $outpre.slurm.e.%j
echo "Mission starts at `date`" &&\\
uname -n &&\\
echo "Unzip mpileup starts at `date`" &&\\
gzip -d $datadir/$normal/$normal.mpileup.gz && gzip -d $datadir/$cancer/$cancer.mpileup.gz &&\\
echo "Varscan call starts at `date`" &&\\
java -Djava.io.tmpdir=$TMP -jar $VarScan somatic $datadir/$normal/$normal.mpileup $datadir/$cancer/$cancer.mpileup --normal-purity 0.98 --tumor-purity 0.90 --min-var-freq 0.05 --output-snp $outpre.snv --output-indel $outpre.indel &&\\
echo "Filter starts at `date`" &&\\
java -Djava.io.tmpdir=$TMP -jar $VarScan processSomatic $outpre.snv &&\\
java -Djava.io.tmpdir=$TMP -jar $VarScan somaticFilter $outpre.snv.Somatic.hc  --indel-file $outpre.indel --min-var-freq 0.05 --output-file $outpre.snv.Somatic.hc.filter &&\\
java -Djava.io.tmpdir=$TMP -jar $VarScan processSomatic $outpre.indel &&\\
awk 'BEGIN{getline}{print \$1"\\t"\$2"\\t"\$2}' $outpre.snv.Somatic.hc.filter > $outpre.snv.Somatic.hc.filter.pos &&\\
awk 'BEGIN{getline}{if(\$4 ~ /^-/){i=\$2+length(\$4);print \$1"\t"\$2"\t"i}else{j=\$2+1;print \$1"\t"\$2"\t"j}}' $outpre.indel.Somatic.hc > $outpre.indel.Somatic.hc.pos &&\\
bam-readcount -b 20 -f $GENOME -l $outpre.snv.Somatic.hc.filter.pos $cancer_bam > $outpre.snv.Somatic.hc.filter.rc &&\\
bam-readcount -b 20 -f $GENOME -l $outpre.indel.Somatic.hc.pos $cancer_bam > $outpre.indel.Somatic.hc.rc &&\\
java -Djava.io.tmpdir=$TMP -jar  $VarScan fpfilter $outpre.snv.Somatic.hc.filter $outpre.snv.Somatic.hc.filter.rc --output-file $outpre.snv.Somatic.hc.filter.pass &&\\
java -Djava.io.tmpdir=$TMP -jar  $VarScan fpfilter $outpre.indel.Somatic.hc $outpre.indel.Somatic.hc.rc --output-file $outpre.indel.Somatic.hc.pass &&\\
#perl $BIN/somatic_process_filter.pl <(cut -f 1-7,9-11,16-23 $outpre.indel.Somatic.hc) > $outpre.indel.process &&\\
cut -f 1-7,9-11,16-23   $outpre.indel.Somatic.hc.pass > $outpre.indel.process &&\\
cut -f 1-7,9-11,16-23   $outpre.snv.Somatic.hc.filter.pass > $outpre.snv.process &&\\
echo "Gzip mpileup starts at `date`" &&\\
gzip $datadir/$normal/$normal.mpileup && gzip $datadir/$cancer/$cancer.mpileup &&\\
echo "Mission completes at `date`"

EOF
	close OUT;
}

if ( $vardict || $all){
	my $soft = "vardict";
	my $dir = "$outdir/$soft";
	my $outpre = "$prefix.$soft";
	if (-d $dir){
		die "ERROR:The dietectory $dir has exist!";
	}else{
		`mkdir $dir`;
	}
	open OUT, ">$dir/$outpre.slurm" or die $!;
	print OUT<<EOF;
#!/bin/bash
#SBATCH -J $outpre
#SBATCH -p cpu
#SBATCH --mail-type=end
#SBATCH --mail-user=ybliu1988\@163.com
#SBATCH -n 8
#SBATCH --ntasks-per-node=8
#SBATCH -o $outpre.slurm.o.%j
#SBATCH -e $outpre.slurm.e.%j
echo "Mission starts at `date`" &&\\
uname -n &&\\
JAVA_HOME=/usr
## WGS by 150bp overlapping 5kb segments
AF_THR="0.03" # minimum allele frequency
$VARDICT_HOME/VarDict -th 16 -I 150 -r 3 -G $GENOME \\
	-f \$AF_THR -N $prefix  \\
	-b "$cancer_bam|$normal_bam" \\
	-c 1 -S 2 -E 3 -g 4 $VARDICT_HOME/hg38_5k_150bpOL_seg.txt \\
	| $VARDICT_HOME/testsomatic.R \\
	| $VARDICT_HOME/var2vcf_paired.pl -N "$cancer|$normal" -f \$AF_THR > $outpre.vcf &&\\
grep "StrongSomatic" $outpre.vcf  | grep "PASS" > $outpre.somatic.vcf &&\\
awk '{if(\$7=="PASS"){print \$0}}' $outpre.somatic.vcf | grep -E "TYPE=Insertion|TYPE=Deletion"  > $outpre.indel.vcf &&\\
awk '{if(\$7=="PASS"){print \$0}}' $outpre.somatic.vcf | grep -E "TYPE=SNV"  > $outpre.snv.vcf &&\\
perl $BIN/vardict_vcf_process.pl $outpre.snv.vcf > $outpre.snv.process &&\\
perl $BIN/vardict_vcf_process.pl $outpre.indel.vcf > $outpre.indel.process &&\\

echo "Mission completes at `date`"
EOF
	close OUT;
}

if ( $somaticsniper || $all){
	my $soft = "somaticsniper";
	my $dir = "$outdir/$soft";
	my $outpre = "$prefix.$soft";
	if (-d $dir){
		die "ERROR:The dietectory $dir has exist!";
	}else{
		`mkdir $dir`;
	}
	open OUT, ">$dir/$outpre.slurm" or die $!;
	print OUT<<EOF;
#!/bin/bash
#SBATCH -J $outpre
#SBATCH -p cpu
#SBATCH --mail-type=end
#SBATCH --mail-user=ybliu1988\@163.com
#SBATCH -n 2
#SBATCH --ntasks-per-node=2
#SBATCH -o $outpre.slurm.o.%j
#SBATCH -e $outpre.slurm.e.%j
echo "Mission starts at `date`" &&\\
uname -n &&\\
bam-somaticsniper -Q 40 -q 20  -F vcf -G -L -f $GENOME $cancer_bam $normal_bam $outpre.vcf &&\\
perl /lustre/home/acct-medkkw/medlyb/soft/somatic-sniper-1.0.5.0/src/scripts/snpfilter.pl -snp-file $outpre.vcf -indel-file  $datadir/$cancer/$cancer.10col.mpileup &&\\
perl /lustre/home/acct-medkkw/medlyb/soft/somatic-sniper-1.0.5.0/src/scripts/prepare_for_readcount.pl -snp-file $outpre.vcf.SNPfilter &&\\
bam-readcount -b 20 -f $GENOME -l $outpre.vcf.SNPfilter.pos $cancer_bam > $outpre.vcf.SNPfilter.rc &&\\
perl /lustre/home/acct-medkkw/medlyb/soft/somatic-sniper-1.0.5.0/src/scripts/fpfilter.pl -snp-file $outpre.vcf.SNPfilter -readcount-file $outpre.vcf.SNPfilter.rc &&\
perl /lustre/home/acct-medkkw/medlyb/soft/somatic-sniper-1.0.5.0/src/scripts/highconfidence.pl -snp-file $outpre.vcf.SNPfilter.fp_pass &&\\
perl $BIN/somaticsniper_vcf_process.pl $outpre.vcf.SNPfilter.fp_pass.hc > $outpre.snv.process &&\\

echo "Mission completes at `date`"
EOF
	close OUT;
}

if ( $muse || $all){
	my $soft = "muse";
	my $dir = "$outdir/$soft";
	my $outpre = "$prefix.$soft";
	if (-d $dir){
		die "ERROR:The dietectory $dir has exist!";
	}else{
		`mkdir $dir`;
	}
	open OUT, ">$dir/$outpre.slurm" or die $!;
	print OUT<<EOF;
#!/bin/bash
#SBATCH -J $outpre
#SBATCH -p cpu
#SBATCH --mail-type=end
#SBATCH --mail-user=ybliu1988\@163.com
#SBATCH -n 2
#SBATCH --ntasks-per-node=2
#SBATCH -o $outpre.slurm.o.%j
#SBATCH -e $outpre.slurm.e.%j
echo "Mission starts at `date`" &&\\
uname -n &&\\

MuSE call -O $outpre -l $BIN/hg38.chrom.list -f $GENOME $cancer_bam $normal_bam &&\\
MuSE sump -I $outpre.MuSE.txt -G -O $outpre.MuSE.vcf  -D /lustre/home/acct-medkkw/medlyb/database/annotation/gatk_ann/hg38/dbsnp_146.hg38.vcf.gz &&\\
perl $BIN/muse_vcf_process.pl $outpre.MuSE.vcf > $outpre.snv.process &&\\

echo "Mission completes at `date`"
EOF
	close OUT;
}

if ( $manta_strelka || $all){
	my $soft = "manta_strelka";
	my $dir = "$outdir/$soft";
	my $outpre = "$prefix.$soft";
	if (-d $dir){
		die "ERROR:The dietectory $dir has exist!";
	}else{
		`mkdir $dir`;
	}
	open OUT, ">$dir/$outpre.slurm" or die $!;
	print OUT<<EOF;
#!/bin/bash
#SBATCH -J $outpre
#SBATCH -p cpu
#SBATCH --mail-type=end
#SBATCH --mail-user=ybliu1988\@163.com
#SBATCH -n 8
#SBATCH --ntasks-per-node=8
#SBATCH -o $outpre.slurm.o.%j
#SBATCH -e $outpre.slurm.e.%j
echo "Mission starts at `date`" &&\\
uname -n &&\\

python /lustre/home/acct-medkkw/medlyb/soft/manta-1.4.0.centos6_x86_64/bin/configManta.py --normalBam=$normal_bam --tumorBam=$cancer_bam --referenceFasta=$GENOME --callRegions=/lustre/home/acct-medkkw/medlyb/script/wgs_pipeline/hg38_chr.bed.gz --runDir=MantaAnalysis_$prefix &&\\
python MantaAnalysis_$prefix/runWorkflow.py -m local -j 8 &&\\
python /lustre/home/acct-medkkw/medlyb/soft/strelka-2.9.9.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py --normalBam=$normal_bam --tumorBam=$cancer_bam --referenceFasta=$GENOME --callRegions=/lustre/home/acct-medkkw/medlyb/script/wgs_pipeline/hg38_chr.bed.gz --indelCandidates=MantaAnalysis_$prefix/results/variants/candidateSmallIndels.vcf.gz --outputCallableRegions --runDir=StrelkaAnalysis_$prefix &&\\
python StrelkaAnalysis_$prefix/runWorkflow.py -m local -j 8 &&\\

perl $BIN/strelka_vcf_process.pl <(zcat StrelkaAnalysis_$prefix/results/variants/somatic.snvs.vcf.gz | grep PASS) > $prefix.strelka.snv.process &&\\
perl $BIN/strelka_vcf_process.pl <(zcat StrelkaAnalysis_$prefix/results/variants/somatic.indels.vcf.gz | grep PASS) > $prefix.strelka.indel.process &&\\

echo "Mission completes at `date`"
EOF
	close OUT;
}

if ( $cnvkit || $all){
	my $soft = "cnvkit";
	my $dir = "$outdir/$soft";
	my $outpre = "$prefix.$soft";
	if (-d $dir){
		die "ERROR:The dietectory $dir has exist!";
	}else{
		`mkdir $dir`;
	}
	open OUT, ">$dir/$outpre.slurm" or die $!;
	print OUT<<EOF;
#!/bin/bash
#SBATCH -J $outpre
#SBATCH -p cpu
#SBATCH --mail-type=end
#SBATCH --mail-user=ybliu1988\@163.com
#SBATCH -n 8
#SBATCH --ntasks-per-node=8
#SBATCH -o $outpre.slurm.o.%j
#SBATCH -e $outpre.slurm.e.%j
echo "Mission starts at `date`" &&\\
uname -n &&\\
cnvkit.py batch $cancer_bam --normal $normal_bam --method wgs --targets $BIN/hg38_refGene.bed --fasta $GENOME --annotate /lustre/home/acct-medkkw/medlyb/database/UCSChg38/refFlat.txt --output-reference $outpre.cnn --output-dir $outpre.out/

echo "Mission completes at `date`"
EOF
	close OUT;
}

if ( $gvcf || $all){
	my $soft = "gvcf";
	my $dir = "$outdir/$soft";
	my $outpre = "$prefix.$soft";
	my $cn_pre = "$cancer.$soft";
	my $nm_pre = "$normal.$soft";
	if (-d $dir){
		die "ERROR:The dietectory $dir has exist!";
	}else{
		`mkdir $dir`;
	}
	open OUT, ">$dir/$cn_pre.slurm" or die $!;
	print OUT<<EOF;
#!/bin/bash
#SBATCH -J $cn_pre
#SBATCH -p cpu
#SBATCH --mail-type=end
#SBATCH --mail-user=ybliu1988\@163.com
#SBATCH -n 2
#SBATCH --ntasks-per-node=2
#SBATCH -o $cn_pre.slurm.o.%j
#SBATCH -e $cn_pre.slurm.e.%j
uname -n &&\\
echo "Mission $cn_pre starts at  `date`" &&\\
java -Djava.io.tmpdir=$TMP -jar $GATK -T HaplotypeCaller -R $GENOME -I $cancer_bam --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -o ${cancer}.g.vcf.gz &&\\
echo "Mission completes at `date`"
EOF
	close OUT;

	open OUT, ">$dir/$nm_pre.slurm" or die $!;
	print OUT<<EOF;
#!/bin/bash
#SBATCH -J $nm_pre
#SBATCH -p cpu
#SBATCH --mail-type=end
#SBATCH --mail-user=ybliu1988\@163.com
#SBATCH -n 4
#SBATCH --ntasks-per-node=4
#SBATCH -o $nm_pre.slurm.o.%j
#SBATCH -e $nm_pre.slurm.e.%j
uname -n &&\\
echo "Mission $nm_pre starts at `date`" &&\\
java -Djava.io.tmpdir=$TMP -jar $GATK -T HaplotypeCaller -R $GENOME -I $normal_bam --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -o ${normal}.g.vcf.gz &&\\
echo "Mission completes at `date`"
EOF
	close OUT;
}

if ( $breakdancer || $all){
	my $soft = "breakdancer";
	my $dir = "$outdir/$soft";
	my $outpre = "$prefix.$soft";
	my $cn_pre = "$cancer.$soft";
	my $nm_pre = "$normal.$soft";
	if (-d $dir){
		die "ERROR:The dietectory $dir has exist!";
	}else{
		`mkdir $dir`;
	}
	open OUT, ">$dir/$cn_pre.slurm" or die $!;
	print OUT<<EOF;
#!/bin/bash
#SBATCH -J $cn_pre
#SBATCH -p cpu
#SBATCH --mail-type=end
#SBATCH --mail-user=ybliu1988\@163.com
#SBATCH -n 4
#SBATCH --ntasks-per-node=4
#SBATCH -o $cn_pre.slurm.o.%j
#SBATCH -e $cn_pre.slurm.e.%j
echo "Mission starts at `date`" &&\\
uname -n &&\\

perl /lustre/home/acct-medkkw/medlyb/bin/bam2cfg.pl -v 20 $cancer_bam > ${cn_pre}.cfg &&\\
breakdancer-max -h  ${cn_pre}.cfg > ${cn_pre}.ctx &&\\
perl /lustre/home/acct-medkkw/medlyb/shareliu/DNA_HiSeqWGS_2015a/breakdancer/breakdancer_v1.0/bin/filt_sv.pl  -m 100 -x 1000000 -s 30 -d 5 -i ${cn_pre}.ctx -o  ${cn_pre}.filter.ctx &&\\
perl ../ctx2avinput.pl ${cn_pre}.filter.ctx &&\\
perl /lustre/home/acct-medkkw/medlyb/soft/annovar/table_annovar.pl  ${cn_pre}.filter.avinput /lustre/home/acct-medkkw/medlyb/soft/annovar/humandb/ -buildver hg38 -protocol refGene,gff3 --gff3dbfile hg38_genehancer.txt   -operation g,r --remove  --outfile ${cn_pre}.filter &&\\
awk 'BEGIN{print "Chr\\tPos\\tINFO"}{OFS="\\t";print \$1,\$2,\$6}' ${cn_pre}.filter.avinput  > ${cn_pre}.filter.avinput.tmp &&\\
paste ${cn_pre}.filter.avinput.tmp ${cn_pre}.filter.hg38_multianno.txt | cut -f 1,2,3,9,10,11,14  > ${cn_pre}.filter.anno.txt &&\\
rm ${cn_pre}.filter.avinput.tmp &&\\

echo "Mission completes at `date`"
EOF
	close OUT;

	open OUT, ">$dir/$nm_pre.slurm" or die $!;
	print OUT<<EOF;
#!/bin/bash
#SBATCH -J $nm_pre
#SBATCH -p cpu
#SBATCH --mail-type=end
#SBATCH --mail-user=ybliu1988\@163.com
#SBATCH -n 4
#SBATCH --ntasks-per-node=4
#SBATCH -o $nm_pre.slurm.o.%j
#SBATCH -e $nm_pre.slurm.e.%j
echo "Mission starts at `date`" &&\\
uname -n &&\\

perl /lustre/home/acct-medkkw/medlyb/bin/bam2cfg.pl -v 20 $normal_bam > ${nm_pre}.cfg &&\\
breakdancer-max -h  ${nm_pre}.cfg > ${nm_pre}.ctx &&\\
perl /lustre/home/acct-medkkw/medlyb/shareliu/DNA_HiSeqWGS_2015a/breakdancer/breakdancer_v1.0/bin/filt_sv.pl  -m 100 -x 1000000 -s 30 -d 5 -i ${nm_pre}.ctx -o  ${nm_pre}.filter.ctx &&\\
perl ../ctx2avinput.pl ${nm_pre}.filter.ctx &&\\
perl /lustre/home/acct-medkkw/medlyb/soft/annovar/table_annovar.pl  ${nm_pre}.filter.avinput /lustre/home/acct-medkkw/medlyb/soft/annovar/humandb/ -buildver hg38 -protocol refGene,gff3 --gff3dbfile hg38_genehancer.txt   -operation g,r --remove  --outfile ${nm_pre}.filter &&\\
awk 'BEGIN{print "Chr\\tPos\\tINFO"}{OFS="\\t";print \$1,\$2,\$6}' ${nm_pre}.filter.avinput  > ${nm_pre}.filter.avinput.tmp &&\\
paste ${nm_pre}.filter.avinput.tmp ${nm_pre}.filter.hg38_multianno.txt | cut -f 1,2,3,9,10,11,14  > ${nm_pre}.filter.anno.txt &&\\
rm ${nm_pre}.filter.avinput.tmp &&\\

echo "Mission completes at `date`"
EOF
	close OUT;
}

if ( $lumpy || $all){
	my $soft = "lumpy";
	my $dir = "$outdir/$soft";
	my $outpre = "$prefix.$soft";
	if (-d $dir){
		die "ERROR:The dietectory $dir has exist!";
	}else{
		`mkdir $dir`;
	}
	open OUT, ">$dir/$outpre.slurm" or die $!;
	print OUT<<EOF;
#!/bin/bash
#SBATCH -J $outpre
#SBATCH -p cpu
#SBATCH --mail-type=end
#SBATCH --mail-user=ybliu1988\@163.com
#SBATCH -n 8
#SBATCH --ntasks-per-node=8
#SBATCH -o $outpre.slurm.o.%j
#SBATCH -e $outpre.slurm.e.%j
echo "Mission starts at `date`" &&\\
uname -n &&\\
lumpyexpress -B ${cancer_bam},${normal_bam} -D $datadir/$cancer/${cancer}.discordants.bam,$datadir/$normal/${normal}.discordants.bam -S $datadir/$cancer/${cancer}.splitters.bam,$datadir/$normal/${normal}.splitters.bam -o ${outpre}.vcf &&\\
svtyper -B ${cancer_bam},${normal_bam} -i ${outpre}.vcf  > ${outpre}.gt.vcf

echo "Mission completes at `date`"
EOF
	close OUT;
}

