#!/bin/bash

#$ -S /bin/bash
#$ -N vp 

#$ -o /work/$USER/$JOB_NAME-$JOB_ID-$TASK_ID.log
#$ -j y

#$ -l h_rt=48:00:00
#$ -l h_vmem=8G

#$ -pe smp 1
#$ -binding linear:4

#$ -m a -M guzman@zoology.ubc.ca

#$ -wd /home/guzmanl

module load R/3.5.1-2

Rscript run_vp_array.R 
