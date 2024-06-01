
### 
###list and id_hap.gr are files with sampleid and group information
### for example, the format of list file is: adbann5,adbanp1,awlaos1,awlaos2,awpak1,awsri3,awsumg1,awview6,ogcebig2,ogcebig3,ogcebig4,sus01,sus02
## id_hap.gr format likes this: 
ogcebig2_A	sus
ogcebig2_B	sus
ogcebig3_A	sus
ogcebig3_B	sus
ogcebig4_A	sus
ogcebig4_B	sus
adbann5_A saw
adbann5_B saw
adbanp1_A saw
adbanp1_B saw
awlaos1_A msea
awlaos1_B msea
awlaos2_A msea
awlaos2_B msea
awpak1_A saw
awpak1_B saw
awsri3_A saw
awsri3_B saw
awsumg1_A isea
awsumg1_B isea
awview6_A msea
awview6_B msea
sus01_A isea
sus01_B isea
sus02_A isea
sus02_B isea
###

for i in {1..18} X; do
    (
        python ./genomics_general-master/VCF_processing/parseVCF.py \
            -i Chr$i.beagle.745.vcf.gz \
            --samples `cat list` \
            --skipIndels | gzip > $i.geno.gz

        python ./genomics_general-master/phyml_sliding_windows.py \
            -T 2 \
            -g $i.geno.gz \
            --prefix $i.phyml_bionj.w10000 \
            -w 10000 \
            --windType sites \
            --model GTR \
            --optimise n

        python ../twisst-master/twisst.py \
            -t $i.phyml_bionj.w10000.trees.gz \
            -w $i.phyml_bionj.w10000.weight.gz \
            --outputTopos $i.topologies.trees \
            -g isea \
            -g saw \
            -g msea \
            -g sus \
            --groupsFile id_hap.gr \
            --method complete \
            --outgroup sus
    ) &
done;