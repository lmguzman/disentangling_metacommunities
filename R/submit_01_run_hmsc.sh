#!/bin/bash

#$ -S /bin/bash
#$ -N hmsc_new 

#$ -o hmsc.out
#$ -e hmsc.err
#$ -j y

#$ -l h_rt=48:00:00
#$ -l h_vmem=8G

#$ -pe smp 4
#$ -binding linear:4

#$ -m a -M guzman@zoology.ubc.ca

#$ -wd /home/guzmanl

module load R/3.5.1-2

echo "loaded module"

Rscript run_hmsc_new_package.R 
