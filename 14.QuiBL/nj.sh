cp x* snplist
plink -bfile ../../sample4 -extract snplist --distance square 1-ibs --out dis
awk '{print $2}' dis.mdist.id | paste - dis.mdist | sed 's/\t/      /' >inf
wc -l inf | awk '{print "  "$1}' | cat - inf |sed 's/rucosus//;s/hog//;s/chi//' >infile
rm -f sortid inf outfile outtree
echo Y | neighbor
mv outtree ibswild.tre
