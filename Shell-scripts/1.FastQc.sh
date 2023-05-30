#!/bin/bash

# Anik Dutta

#SBATCH --job-name=Fastqc
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10000
#SBATCH --time=48:00:00
#SBATCH --output=Fastqc.out
#SBATCH --error=Fastqc.err
#SBATCH --partition=cluster
#SBATCH --export=NONE

# Load module fastqc to produce quality stats and figures for the raw fastq files.

module load fastqc

for file in *.fq.gz
do
fastqc $file -o /gxfs_work1/cau/suaph275/Tomato_Illumina_Hiseq/FastQC_results
done
