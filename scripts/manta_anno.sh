zcat somaticSV.vcf.gz | grep -v "^#" |  awk '{OFS="\t";print $1,$2,$2,0,0}'  > somaticSV.avinput
perl /lustre/home/acct-medkkw/medlyb/soft/annovar/table_annovar.pl somaticSV.avinput /lustre/home/acct-medkkw/medlyb/soft/annovar/humandb/ -buildver hg38 -protocol refGene --gff3dbfile hg38_genehancer.txt   -operation g --remove  --outfile somaticSV
paste <(awk '{OFS="\t";print $7,$6}'   somaticSV.hg38_multianno.txt) <(zcat somaticSV.vcf.gz | grep -v "^##"  | sed 's/^#//') > somaticSV.sv.txt

#zcat diploidSV.vcf.gz | grep -v "^#" |  awk '{OFS="\t";print $1,$2,$2,0,0}'  > diploidSV.avinput
#perl /lustre/home/acct-medkkw/medlyb/soft/annovar/table_annovar.pl diploidSV.avinput /lustre/home/acct-medkkw/medlyb/soft/annovar/humandb/ -buildver hg38 -protocol refGene --gff3dbfile hg38_genehancer.txt -operation g --remove  --outfile diploidSV
#paste <(awk '{OFS="\t";print $7,$6}'   diploidSV.hg38_multianno.txt) <(zcat diploidSV.vcf.gz | grep -v "^##"  | sed 's/^#//') > diploidSV.sv.txt

#zcat candidateSV.vcf.gz | grep -v "^#" |  awk '{OFS="\t";print $1,$2,$2,0,0}'  > candidateSV.avinput
#perl /lustre/home/acct-medkkw/medlyb/soft/annovar/table_annovar.pl candidateSV.avinput /lustre/home/acct-medkkw/medlyb/soft/annovar/humandb/ -buildver hg38 -protocol refGene --gff3dbfile hg38_genehancer.txt -operation g --remove  --outfile candidateSV
#paste <(awk '{OFS="\t";print $7,$6}'   candidateSV.hg38_multianno.txt) <(zcat candidateSV.vcf.gz | grep -v "^##"  | sed 's/^#//') > candidateSV.sv.txt
