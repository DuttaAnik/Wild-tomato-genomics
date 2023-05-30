#!/bin/bash

# Anik Dutta
# This bash script is used to perform analysis of genetic variation data using ANGSD (Analysis of Next Generation Sequencing Data) and RealSFS, 
# generating multiple statistics, such as allele frequencies, estimates of the site frequency spectrum (SFS), and Watterson's theta.

# The analysis is applied to a series of datasets, each denoted by a different number (e.g., 1958, 1963, etc.)

# ANGSD command options:
# -bam: Input bam file list
# -P: Number of threads
# -dosaf: Calculate the Site Frequency Spectrum
# -doCounts: Count the number A,C,G,T. All sites, All samples
# -GL: Genotype likelihood estimation method
# -remove_bads: Remove bad reads. Uses internal flag.
# -uniqueOnly: Discards reads that don't map uniquely
# -only_proper_pairs: Only use reads where the mate could be mapped as well
# -C: Adjust mapQ for excessive mismatches
# -baq: Perform base alignment quality (BAQ) computation
# -anc and -ref: Ancestral and reference states
# -minMapQ: Minimum mapping quality
# -minQ: Minimum base quality
# -minInd: Discard site if effective sample size below value
# -setMinDepth and -setMaxDepth: Min and max depth for site to be included
# -out: Output file prefix

for i in 1958 1963 2747 2931 2932 3111 4107 4117 4330
do
angsd -bam "$i".txt -P 32 -dosaf 1 -doCounts 1 -GL 1 -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -C 50 -baq 1 -anc S_chilense_Hirise.fasta -ref S_chilense_Hirise.fasta -minMapQ 30 -minQ 30 -minInd 9 -setMinDepth 20 -setMaxDepth 600 -out $i
done

# realSFS: Estimation of the SFS using the method of moments

for i in 1958 1963 2747 2931 2932 3111 4107 4117 4330
do
echo $i
realSFS $i.saf.idx -P 32 > $i.sfs
done

# Watterson's theta estimation using ANGSD
# -doThetas: estimate thetas
# -GL: Genotype likelihood estimation method
# -pest: Input external sfs estimate

for i in 1958 1963 2747 2931 2932 3111 4107 4117 4330
do
echo $i
angsd -P 16 -bam $i.txt -ref S_chilense_Hirise.fasta -anc S_chilense_Hirise.fasta -out $i -doThetas 1 -doSaf 1 -GL 2 -minMapQ 30 -minQ 20 -minInd 8 -pest $i.sfs
done

# Performing a sliding-window analysis of theta with thetaStat
# thetaStat is part of the ANGSD package and can estimate several population genetics parameters

for i in 1958 1963 2747 2931 2932 3111 4107 4117 4330
do
	echo $i
	thetaStat do_stat $i.thetas.idx &> /dev/null
	# perform a sliding-window analysis
	thetaStat do_stat $i.thetas.idx -win 10000 -step 5000 -outnames $i.thetas &> /dev/null
done
