#!/bin/bash
#CSUB -J myjobDXY
#CSUB -q c01
#CSUB -o %J.out
#CSUB -e %J.error
#CSUB -n 8
#CSUB -R "span[hosts=1]"
source ~/.bashrc
conda activate py27 
python2 QuIBL/QuIBL.py QuIBL_pygmy.cfg
