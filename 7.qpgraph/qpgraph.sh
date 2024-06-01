### qpGraph ##

### The online version of ADMIXTOOLS 2 was used in R platform. Only 50K SNPs were randomly picked for analyses.
plink -bfile xxx -thin-count 50000 
##In R, execute the following code:
library(admixtools)
library(tidyverse)
admixtools::run_shiny_admixtools()


##Using Warthog as an outgroup. Optimize the model until the WR |Z| < 3. Change the saved tsv file to pdf, in shell.
dot -Tpdf ${i}.dot -o ${i}.pdf 
