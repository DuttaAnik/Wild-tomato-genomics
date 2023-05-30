#!/bin/bash

# Anik Dutta

#SBATCH --job-name=Admixture_0.2
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=120G
#SBATCH --qos=long
#SBATCH --time=10-00:00:00
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=anik.dutta@phytomed.uni-kiel.de
#SBATCH --partition=cluster


# Admixture analysis using a Pruned dataset which contain unlinked SNPs (r2=0.2) since SNPs in high LD affects the analysis. First, the main VCF was pruned following the procedure here https://speciationgenomics.github.io/pca/

## Pruning VCF file for admixture analysis (from speciation genomics website)

plink --vcf Chilense_99.genotyped.SNP.Clean.MAF.MAXMiss.MinDP.recode.dustmasked.headers.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.2 --out prunning

# make bed file with pruned data
plink --vcf Chilense_99.genotyped.SNP.Clean.MAF.MAXMiss.MinDP.recode.dustmasked.headers.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --extract prunning.prune.in --make-bed --out prunned_data

# ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0

awk '{$1="0";print $0}' prunned_data.bim > prunned_data.bim.tmp
mv prunned_data.bim.tmp prunned_data.bim

# Run Admixture

for K in {1..10} 
do 
	for repeat in {1..10} 
	do 
		admixture --cv Pruned.0.2.bed $K -j24 -s ${repeat} | tee log${K}.${repeat}.out
	done
done
