#PREFIX=$1
#indel.Somatic.hc snv.Somatic.hc.filter

cut -f 1-7,9-11,16-23   VarScan_CHH_DO-MN.indel.Somatic.hc.pass > CHH_D0-MN.varscan.indel.process
cut -f 1-7,9-11,16-23   VarScan_CHH_DO-MN.snv.Somatic.hc.filter.pass > CHH_D0-MN.varscan.snv.process
