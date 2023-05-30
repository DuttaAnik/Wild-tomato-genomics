#!/bin/bash

# Anik Dutta

#SBATCH --job-name=LD_decay
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G
#SBATCH --qos=long
#SBATCH --time=10-00:00:00
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=anik.dutta@phytomed.uni-kiel.de
#SBATCH --partition=cluster

# Script for estimating LD decay.

# -MaxDist specifies the maximum distance (in bases) between variants to consider when calculating LD. Here, it is set to 500 bases.
# -Het specifies the maximum heterozygosity rate allowed for a SNP to be included in the analysis. Here, it is set to 1, meaning all SNPs (regardless of their heterozygosity rate) will be considered.
# -Miss specifies the maximum missing rate allowed for a SNP to be included in the analysis. Here, it is set to 0.2, meaning SNPs with a missing rate less than or equal to 20% will be considered.

PopLDdecay -InVCF Chilense_99.genotyped.SNP.Clean.MAF.MAXMiss.MinDP.recode.dustmasked.headers.vcf -MaxDist 500 -Het 1 -Miss 0.2 -OutStat 1958_pop -SubPop 1958_pop.txt

PopLDdecay -InVCF Chilense_99.genotyped.SNP.Clean.MAF.MAXMiss.MinDP.recode.dustmasked.headers.vcf -MaxDist 1000 -Het 1 -Miss 0.2 -OutStat 1958_pop_1000 -SubPop 1958_pop.txt

PopLDdecay -InVCF Chilense_99.genotyped.SNP.Clean.MAF.MAXMiss.MinDP.recode.dustmasked.headers.vcf -MaxDist 500 -Het 1 -Miss 0.2 -OutStat 1963_pop -SubPop 1963_pop.txt

PopLDdecay -InVCF Chilense_99.genotyped.SNP.Clean.MAF.MAXMiss.MinDP.recode.dustmasked.headers.vcf -MaxDist 1000 -Het 1 -Miss 0.2 -OutStat 1963_pop_1000 -SubPop 1963_pop.txt

PopLDdecay -InVCF Chilense_99.genotyped.SNP.Clean.MAF.MAXMiss.MinDP.recode.dustmasked.headers.vcf -MaxDist 500 -Het 1 -Miss 0.2 -OutStat 2747_pop -SubPop 2747_pop.txt

PopLDdecay -InVCF Chilense_99.genotyped.SNP.Clean.MAF.MAXMiss.MinDP.recode.dustmasked.headers.vcf -MaxDist 1000 -Het 1 -Miss 0.2 -OutStat 2747_pop_1000 -SubPop 2747_pop.txt

PopLDdecay -InVCF Chilense_99.genotyped.SNP.Clean.MAF.MAXMiss.MinDP.recode.dustmasked.headers.vcf -MaxDist 500 -Het 1 -Miss 0.2 -OutStat 2931_pop -SubPop 2931_pop.txt

PopLDdecay -InVCF Chilense_99.genotyped.SNP.Clean.MAF.MAXMiss.MinDP.recode.dustmasked.headers.vcf -MaxDist 1000 -Het 1 -Miss 0.2 -OutStat 2931_pop_1000 -SubPop 2931_pop.txt

PopLDdecay -InVCF Chilense_99.genotyped.SNP.Clean.MAF.MAXMiss.MinDP.recode.dustmasked.headers.vcf -MaxDist 500 -Het 1 -Miss 0.2 -OutStat 2932_pop -SubPop 2932_pop.txt

PopLDdecay -InVCF Chilense_99.genotyped.SNP.Clean.MAF.MAXMiss.MinDP.recode.dustmasked.headers.vcf -MaxDist 1000 -Het 1 -Miss 0.2 -OutStat 2932_pop_1000 -SubPop 2932_pop.txt

PopLDdecay -InVCF Chilense_99.genotyped.SNP.Clean.MAF.MAXMiss.MinDP.recode.dustmasked.headers.vcf -MaxDist 500 -Het 1 -Miss 0.2 -OutStat 3111_pop -SubPop 3111_pop.txt

PopLDdecay -InVCF Chilense_99.genotyped.SNP.Clean.MAF.MAXMiss.MinDP.recode.dustmasked.headers.vcf -MaxDist 1000 -Het 1 -Miss 0.2 -OutStat 3111_pop_1000 -SubPop 3111_pop.txt

PopLDdecay -InVCF Chilense_99.genotyped.SNP.Clean.MAF.MAXMiss.MinDP.recode.dustmasked.headers.vcf -MaxDist 500 -Het 1 -Miss 0.2 -OutStat 4107_pop -SubPop 4107_pop.txt

PopLDdecay -InVCF Chilense_99.genotyped.SNP.Clean.MAF.MAXMiss.MinDP.recode.dustmasked.headers.vcf -MaxDist 1000 -Het 1 -Miss 0.2 -OutStat 4107_pop_1000 -SubPop 4107_pop.txt

PopLDdecay -InVCF Chilense_99.genotyped.SNP.Clean.MAF.MAXMiss.MinDP.recode.dustmasked.headers.vcf -MaxDist 500 -Het 1 -Miss 0.2 -OutStat 4117_pop -SubPop 4117_pop.txt

PopLDdecay -InVCF Chilense_99.genotyped.SNP.Clean.MAF.MAXMiss.MinDP.recode.dustmasked.headers.vcf -MaxDist 1000 -Het 1 -Miss 0.2 -OutStat 4117_pop_1000 -SubPop 4117_pop.txt

PopLDdecay -InVCF Chilense_99.genotyped.SNP.Clean.MAF.MAXMiss.MinDP.recode.dustmasked.headers.vcf -MaxDist 500 -Het 1 -Miss 0.2 -OutStat 4330_pop -SubPop 4330_pop.txt

PopLDdecay -InVCF Chilense_99.genotyped.SNP.Clean.MAF.MAXMiss.MinDP.recode.dustmasked.headers.vcf -MaxDist 1000 -Het 1 -Miss 0.2 -OutStat 4330_pop_1000 -SubPop 4330_pop.txt

# Making LD plot

perl /gxfs_home/cau/suaph275/PopLDdecay/bin/Plot_MultiPop.pl -inList Multi_pop.list -output all_pop_500 -keepR

ls `pwd`/*_pop.stat.gz > pop.list












