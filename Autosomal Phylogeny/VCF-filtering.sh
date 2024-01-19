target=/GATK-filter/target
mkdir -p $target
CHR=$1
RAW_vcf=745.raw.vcf.gz
header=name.conversion.WithGvcfnames.header.745
prefix="chr$CHR"  
date  
#tabix -h ${RAW_vcf}  ${CHR}:1-50000 |  bgzip > $target/$prefix.vcf.gz
#tabix -p vcf $prefix.vcf.gz 
 
BCFTOOLS_Filter_SNP=${prefix}.BCFTOOLS_Filter_SNP.vcf.gz 
 tabix -h ${RAW_vcf}  ${CHR}  |   bcftools filter  -s SnpGap --SnpGap 3 --set-GTs .   \
   -s lowQUAL -e 'QD<2.0 || MQRankSum < -12.5 || FS > 60.0 || ReadPosRankSum < -8.0 || MQ < 40.0 || SOR > 3.0  || AF < 0.003
' - \
  | bcftools annotate -x INFO,^FORMAT/GT,FORMAT/PL -  | bcftools reheader -s  $header - \
  | bcftools view -f PASS -v snps  -m2 -M2 -i 'F_MISSING<0.2'  -O z  -o $target/$BCFTOOLS_Filter_SNP -
tabix -p vcf $target/$BCFTOOLS_Filter_SNP  
date  
 # | grep -v "^#" | cut -f1-10 | more 