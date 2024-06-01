


~/soft/bcftools1.10/bcftools mpileup --threads 10 -C 50 -f ref mkdup.bam|~/soft/bcftools1.10/bcftools call -c --threads 10  |~/soft/bcftools1.10/misc/vcfutils.pl vcf2fq -d 5 -D 100| gzip > ${PWD##*/}.bwaallref.fq.gz ;

~/soft/psmc/utils/fq2psmcfa -q20 ${PWD##*/}.bwaallref.fq.gz >${PWD##*/}.bwaallref.psmcfa;

~/soft/psmc/psmc  -N25 -t15 -r5 -p "4+25*2+4+6" -o  ${PWD##*/}.bwaallref.psmc ${PWD##*/}.bwaallref.psmcfa
