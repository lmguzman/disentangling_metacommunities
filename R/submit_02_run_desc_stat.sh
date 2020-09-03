#!/bin/bash

#$ -S /bin/bash
#$ -N SS 

#$ -o SS.out
#$ -e SS.err
#$ -j y

#$ -l h_rt=48:00:00
#$ -l h_vmem=8G

#$ -pe smp 1

#$ -m a -M guzman@zoology.ubc.ca

#$ -wd /home/guzmanl

module load R/3.5.1-2

Rscript run_ss_array.R 
