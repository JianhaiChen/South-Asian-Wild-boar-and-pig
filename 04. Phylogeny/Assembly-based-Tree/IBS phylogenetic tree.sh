### The phylogenetic trees in Figure 2f and Supplementary Figure 7 were generated with the following steps

#step 1 
#remove repeat regions defined by GCF_000003025.6_Sscrofa11.1_rm.out.gz
plink --bfile all745chr --exclude range repeat.bed --make-bed --out 745sample.nonrep

#step 2
#filtering the LD
plink -bfile 745sample.nonrep -indep 50 5 2 -out indep.nonog
plink -bfile 745sample.nonrep -extract indep.nonog.prune.in -make-bed -out 745sample.nonrep.ld

#step 3
#filtering by maf
plink -bfile 745sample.nonrep.ld -maf 0.05 -make-bed -out 745sample.nonrep.ld.maf0.05


#step 4
#focusing on SNPs within regions of specific recombination range (ranking by decile)
plink -bfile 745sample.nonrep.ld.maf0.05 -extract range xx.recombinationfile --distance square 1-ibs --out xx.recombinationfile  -chr 1-18

#step 5
#IBS matrix and neighbor-joining tree
awk '{print $1}'  xx.recombinationfile.mdist.id|paste - xx.recombinationfile.mdist|sed 's/\t/      /' >inf;echo  " " " " `wc -l inf|awk '{print $1}'` |cat - inf >infile;rm sortid inf outfile outtree;echo -ne 'Y'|neighbor;mv outtree all745.xx.recombinationfile.tre 



### The phylogenetic trees in Figure 3a were generated with the above code except the change in Step 3 (-chr 19).