for i in {1..18} X; do python -m momi.read_vcf     Chr$i.beagle.745.vcf.gz popsimp1  $i.snpAlleleCountscebouttop3.gz --bed $i.newneutral.top3.region  --outgroup sceb  & done;
wait;
  python -m momi.extract_sfs sfs.gz 100 *.snpAlleleCountscebouttop3.gz&


###this is samples and group
cat popsimp1 
ogcebig8 sceb
ogcebig1 sceb
ogcebig3 sceb
ogcebig6 sceb
awlaos1 mseaw
awlaos2 mseaw
awlaos3 mseaw
awpak1 SAS
awsumg1 ISEAWB
awsumg2 ISEAWB
ogsverrucosus1 sus
ogsverrucosus2 sus
ogverg1 sus
ogverg2 sus
awbanw3 SAS
adbann2 SAS
adbanp1 SAS
adbant1 SAS
