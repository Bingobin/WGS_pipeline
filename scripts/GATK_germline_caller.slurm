#!/bin/bash

#SBATCH -J GATK_germline
#SBATCH -p cpu
#SBATCH --mail-type=end
#SBATCH --mail-user=ybliu1988@163.com
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -o GATK_germline.slurm.o.%j
#SBATCH -e GATK_germline.slurm.e.%j


PREFIX=$1
GENOME=/lustre/home/acct-medkkw/medlyb/database/annotation/gatk_ann/hg38/bwaindex/Homo_sapiens_assembly38.fasta
ANNO=/lustre/home/acct-medkkw/medlyb/database/annotation/gatk_ann/hg38/hg38bundle
TMP=/lustre/home/acct-medkkw/medlyb/tmp
INTERVAL=/lustre/home/acct-medkkw/medlyb/script/wgs_pipeline/hg38_chr.interval_list
VCF_LIST=""
VCF_LIST="$VCF_LIST -V /lustre/home/acct-medkkw/medlyb/project/$PREFIX/WGS/02.somatic_D0/gvcf/${PREFIX}_D0.g.vcf.gz"
VCF_LIST="$VCF_LIST -V /lustre/home/acct-medkkw/medlyb/project/$PREFIX/WGS/02.somatic_D0/gvcf/${PREFIX}_MN.g.vcf.gz"

gatk CombineGVCFs $VCF_LIST -R $GENOME -O ${PREFIX}.g.vcf.gz --tmp-dir=$TMP -L $INTERVAL &&\
gatk --java-options "-Xmx5g" GenotypeGVCFs -R $GENOME -V ${PREFIX}.g.vcf.gz -O ${PREFIX}.vcf.gz --tmp-dir=$TMP --dbsnp $ANNO/dbsnp_146.hg38.vcf &&\

#################SNP
gatk SelectVariants -R $GENOME -V ${PREFIX}.vcf.gz --tmp-dir=$TMP -O ${PREFIX}.raw.snp.vcf.gz --exclude-non-variants --select-type-to-include SNP &&\
gatk VariantRecalibrator -R $GENOME -V ${PREFIX}.raw.snp.vcf.gz \
	 --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $ANNO/hapmap_3.3.hg38.vcf.gz  \
	 --resource:omni,known=false,training=true,truth=false,prior=12.0 $ANNO/1000G_omni2.5.hg38.vcf.gz \
	 --resource:1000G,known=false,training=true,truth=false,prior=10.0 $ANNO/1000G_phase1.snps.high_confidence.hg38.vcf \
	 --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $ANNO/../dbsnp_151.hg38.vcf.gz \
	 -an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum \
	 -mode SNP \
	 --output ${PREFIX}.recalibrate_SNP.recal \
	 --tranches-file ${PREFIX}.recalibrate_SNP.tranches \
	 --rscript-file ${PREFIX}.recalibrate_SNP_plots.R &&\
gatk ApplyVQSR -R $GENOME -V ${PREFIX}.raw.snp.vcf.gz --tranches-file ${PREFIX}.recalibrate_SNP.tranches --recal-file ${PREFIX}.recalibrate_SNP.recal -O ${PREFIX}.VQSR.snp.vcf.gz -mode SNP &&\
rm ${PREFIX}.raw.snp.vcf.gz ${PREFIX}.raw.snp.vcf.gz.tbi &&\
#gatk SelectVariants -R $GENOME -V ${PREFIX}.filter.snp.vcf.gz --exclude-filtered -O ${PREFIX}.filtered.snp.vcf.gz &&\

################INDEL
gatk SelectVariants -R $GENOME -V ${PREFIX}.vcf.gz --tmp-dir=$TMP -O ${PREFIX}.raw.indel.vcf.gz --exclude-non-variants --select-type-to-include INDEL &&\
gatk VariantRecalibrator -R $GENOME -V ${PREFIX}.raw.indel.vcf.gz \
		 --resource:mills,known=true,training=true,truth=true,prior=12.0 $ANNO/Mills_and_1000G_gold_standard.indels.hg38.vcf  \
		 --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $ANNO/../dbsnp_151.hg38.vcf.gz \
		 -an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum --max-gaussians 4 \
		 -mode INDEL \
		 -O ${PREFIX}.recalibrate_INDEL.recal \
		 --tranches-file ${PREFIX}.recalibrate_INDEL.tranches \
		 --rscript-file ${PREFIX}.recalibrate_INDEL_plots.R &&\
gatk ApplyVQSR -R $GENOME -V ${PREFIX}.raw.indel.vcf.gz --tranches-file ${PREFIX}.recalibrate_INDEL.tranches --recal-file ${PREFIX}.recalibrate_INDEL.recal -O ${PREFIX}.VQSR.indel.vcf.gz -mode INDEL &&\
rm ${PREFIX}.raw.indel.vcf.gz  ${PREFIX}.raw.indel.vcf.gz.tbi &&\
#gatk SelectVariants -R $GENOME -V ${PREFIX}.filter.indel.vcf.gz --exclude-filtered -O ${PREFIX}.filtered.indel.vcf.gz

#################sort and merge
gatk SortVcf -I ${PREFIX}.VQSR.indel.vcf.gz -O ${PREFIX}.VQSR.indel.sort.vcf.gz &&\
gatk SortVcf -I ${PREFIX}.VQSR.snp.vcf.gz -O ${PREFIX}.VQSR.snp.sort.vcf.gz &&\
gatk MergeVcfs -I ${PREFIX}.VQSR.snp.sort.vcf.gz -I ${PREFIX}.VQSR.indel.sort.vcf.gz -O ${PREFIX}.VQSR.sort.vcf.gz &&\
rm ${PREFIX}.VQSR.indel.vcf.gz ${PREFIX}.VQSR.indel.vcf.gz.tbi &&\
rm ${PREFIX}.VQSR.indel.sort.vcf.gz   ${PREFIX}.VQSR.indel.sort.vcf.gz.tbi &&\
rm ${PREFIX}.VQSR.snp.vcf.gz ${PREFIX}.VQSR.snp.vcf.gz.tbi &&\
rm ${PREFIX}.VQSR.snp.sort.vcf.gz ${PREFIX}.VQSR.snp.sort.vcf.gz.tbi &&\
#################ANNO
zcat ${PREFIX}.VQSR.sort.vcf.gz | /usr/bin/java -jar /lustre/home/acct-medkkw/medlyb/soft/snpEff_latest_core/snpEff/SnpSift.jar annotate /lustre/home/acct-medkkw/medlyb/database/annotation/gatk_ann/hg38/dbsnp_151.hg38.vcf.gz  > ${PREFIX}.VQSR.sort.snp151.vcf &&\
perl /lustre/home/acct-medkkw/medlyb/soft/annovar/table_annovar.pl ${PREFIX}.VQSR.sort.snp151.vcf /lustre/home/acct-medkkw/medlyb/soft/annovar/humandb/ -buildver hg38 -protocol refGene,rmsk,gff3,cosmic87,avsnp150 --gff3dbfile hg38_genehancer.txt -operation g,r,r,f,f -vcfinput --remove --outfile ${PREFIX}.VQSR.sort.snp151 &&\
gatk IndexFeatureFile -F ${PREFIX}.VQSR.sort.snp151.hg38_multianno.vcf
