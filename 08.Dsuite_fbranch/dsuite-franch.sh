###VCF files preparation: both autosome and x

###1) for vcf of “top 2” tree 
plink -vcf samplewild2.top2.mafp5genop20onlytop21k.vcf -thin-count 1000000  --recode vcf-iid  --allow-extra-chr   --out samplewild3.top2.mafp5genop20.1m


###2) for vcf of “top 3” tree. Due to a lot of variants underlying this topology, we picked randomly only 1 million.
plink -vcf samplewild2.top3.mafp5genop20.vcf.vcf -thin-count 1000000  --recode vcf-iid  --allow-extra-chr   --out samplewild3.top2.mafp5genop20.1m



#step 1. assign population for individuals.
#cat sample.grwild3

awiran1	EUW
ewgreg3	EUW
ewitag1	EUW
ewnetg4	EUW
awsumg1	ISEA
awsumg2	ISEA
awsumg3	ISEA
awpak1	Lsa
awsri1	Lsa
awsri3	Lsa
awlaos2	MSEA
awlaos3	MSEA
awschig8	MSEA
awview1	MSEA
awjpwbg1	NAW
awkoreg6	NAW
awnchig5	NAW
awneag1	NAW
awnep15F	Nepalwild
awnep16F	Nepalwild
awnep9F	Nepalwild
ogsver1	OGVER
ogsver2	OGVER
ogverg1	OGVER
ogverg2	OGVER
ogBabyrs	Outgroup
pygmyhog1	Pygmy
pygmyhog2	Pygmy
pygmyhog3	Pygmy
pygmyhog4	Pygmy
pygmyhog5	Pygmy
pygmyhog6	Pygmy
ogafr1	Warthog
ogafr2	Warthog
ogafr3	Warthog
ogafr4	Warthog
ogafr6	Warthog

##step 2. manually edit a tree file

###cat sample.grwild3.tre
(Outgroup,(Warthog,(Pygmy,(OGVER,(ISEA,((Lsa,Nepalwild),(EUW,(NAW,MSEA))))))));

step 3. f4-ratio statistics with Dtrios in Dsuite

~/00.soft/Dsuite/Build/Dsuite Dtrios -c -n gr3top3 -t sample.grwild3.tre samplewild2.top3.mafp5genop20.1m.vcf sample.grwild3

step 4. f-branch statistics 
~/00.soft/Dsuite/Build/Dsuite Fbranch sample.grwild3.tre gr3top3_tree.txt > geneflow_Fbranchx.wildgr3.top3.txt

step 5. plot
python3 ~/00.soft/Dsuite/utils/dtools.py geneflow_Fbranchx.wildgr3.top3.txt sample.grwild3.tre
