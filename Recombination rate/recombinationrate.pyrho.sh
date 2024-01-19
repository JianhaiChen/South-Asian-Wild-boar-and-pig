1.Computing population history with Smc++
## For calculating population history with the X chromosome, use the 'estimate' parameter, and multiply the population history obtained from the CSV file by 3/4; For calculating autosomal chromosome population history, use the 'cv' parameter.
##smc++
smc++ vcf2smc -c 1000000 /home/cyxiao/00newpig_2020_june/snpdata/743sample_chr/chrX.vcf.gz smc/chrX.smc.gz X pop:awbanw6,adbann2,adbann5,awpak1,adbant4,adnephw15,adnepcm1,awnep4F,adnephe11,adnephe6,awsri3,awbanw2
smc++ estimate -o pop 2.5e-8 smc/*smc.gz --polarization-error 0.5
smc++ plot pop.pdf --csv pop/model.final.json



2.get sample vcf files with vcftools
##get_indv_vcf
vcftools --gzvcf /home/cyxiao/00newpig_2020_june/snpdata/743sample_chr/chrX.vcf.gz --keep keep --recode --recode-INFO-all --stdout | bgzip > chrX.vcf.gz


##pyrho
3.pyrho
Using lookuptable from smc++ 
##make_table
pyrho make_table -n 24 -N 30 --mu 2.5e-8 --logfile . --outfile pop_lookuptable.hdf --approx --smcpp_file pop.csv --decimate_rel_tol 0.1 --numthreads 10

##hyperparam
pyrho hyperparam -n 24 --mu 2.5e-8 --blockpenalty 10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100 --windowsize 10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100 --logfile . --tablefile pop_lookuptable.hdf --num_sims 10 --smcpp_file pop.csv --outfile pop_hyperparam_results.txt --numthreads 10

##cal_rec
pyrho optimize --tablefile pop_lookuptable.hdf --vcffile chrX.vcf.gz --ploidy 2 --outfile chrX.rmap --blockpenalty 100 --windowsize 100 --logfile . --numthreads 3
