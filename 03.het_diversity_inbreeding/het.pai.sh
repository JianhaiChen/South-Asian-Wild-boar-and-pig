plink --vcf  S5.743.autosome.bcftool.vcf.gz  --het --out 745.Her

vcftools --gzvcf S5.743.autosome.bcftool.vcf.gz --keep ISEA.sampleName       --window-pi 50000 --window-pi-step 50000 --out ISEA.sampleName 