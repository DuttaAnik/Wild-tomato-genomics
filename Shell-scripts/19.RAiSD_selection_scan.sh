#!/bin/bash

# Anik Dutta
# RAiSD selection scan method estimates the μ statistic, a composite evaluation test that scores genomic regions by quantifying changes in the SFS, the levels of LD, and the amount of genetic diversity along a chromosome. 
# I am running the RAiSD using the full vcf now which contains 30080652 SNPs (merged.vcf with 12 scaffolds). And I also discarded the missing SNPs.
# -y is ploidy level (diploid 2) 
# -M 0 (discards SNP missing data)
# -c 1 (Provides the slack for the SFS edges to be used for the calculation of mu_SFS. The default value is 1 (singletons and S-1 snp class, where S is the sample size))
# -w 50 (Provides the window size (integer value). The default value is 50 (empirically determined).
# -A 0.999 (significance threshold to define regions under selection, here 99.9%). With this option, a manhattan plot is generated automatically.

RAiSD -n 1958_C -f -y 2 -M 0 -R -D -c 1 -w 50 -I /gxfs_work1/cau/suaph275/Clean_VCF/GWAS/merged.vcf -A 0.999 -S /gxfs_work1/cau/suaph275/Clean_VCF/ABBA_pop_stats/sel_scan/1958_pop.txt

RAiSD -n 1963_C -f -y 2 -M 0 -R -D -c 1 -w 50 -I /gxfs_work1/cau/suaph275/Clean_VCF/GWAS/merged.vcf -A 0.999 -S /gxfs_work1/cau/suaph275/Clean_VCF/ABBA_pop_stats/sel_scan/1963_pop.txt

RAiSD -n 2931_C -f -y 2 -M 0 -R -D -c 1 -w 50 -I /gxfs_work1/cau/suaph275/Clean_VCF/GWAS/merged.vcf -A 0.999 -S /gxfs_work1/cau/suaph275/Clean_VCF/ABBA_pop_stats/sel_scan/2931_pop.txt

RAiSD -n 2932_SC -f -y 2 -M 0 -R -D -c 1 -w 50 -I /gxfs_work1/cau/suaph275/Clean_VCF/GWAS/merged.vcf -A 0.999 -S /gxfs_work1/cau/suaph275/Clean_VCF/ABBA_pop_stats/sel_scan/2932_pop.txt

RAiSD -n 4107_SC -f -y 2 -M 0 -R -D -c 1 -w 50 -I /gxfs_work1/cau/suaph275/Clean_VCF/GWAS/merged.vcf -A 0.999 -S /gxfs_work1/cau/suaph275/Clean_VCF/ABBA_pop_stats/sel_scan/4107_pop.txt

RAiSD -n Two_C -f -y 2 -M 0 -R -D -c 1 -w 50 -I /gxfs_work1/cau/suaph275/Clean_VCF/GWAS/merged.vcf -A 0.999 -S /gxfs_work1/cau/suaph275/Clean_VCF/ABBA_pop_stats/sel_scan/two_central.txt

RAiSD -n Two_SH -f -y 2 -M 0 -R -D -c 1 -w 50 -I /gxfs_work1/cau/suaph275/Clean_VCF/GWAS/merged.vcf -A 0.999 -S /gxfs_work1/cau/suaph275/Clean_VCF/ABBA_pop_stats/sel_scan/two_highland.txt


# Below is an example if you wish to run the RAiSD using individual scaffold/chromosome separately.

# Chromosome_1/Scaffold_1

RAiSD -n 1958_SC_1 -f -y 2 -M 0 -R -c 1 -w 50 -I /gxfs_work1/cau/suaph275/Clean_VCF/ABBA_pop_stats/sel_scan/scaffolds_imputed/Scaffold_1__2_contigs__length_78070536.vcf -B 78070536 2549116 -S /gxfs_work1/cau/suaph275/Clean_VCF/ABBA_pop_stats/sel_scan/1958_pop.txt

RAiSD -n 1963_SC_1 -f -y 2 -M 0 -R -c 1 -w 50 -I /gxfs_work1/cau/suaph275/Clean_VCF/ABBA_pop_stats/sel_scan/scaffolds_imputed/Scaffold_1__2_contigs__length_78070536.vcf -B 78070536 2549116 -S /gxfs_work1/cau/suaph275/Clean_VCF/ABBA_pop_stats/sel_scan/1963_pop.txt

RAiSD -n 2931_SC_1 -f -y 2 -M 0 -R -c 1 -w 50 -I /gxfs_work1/cau/suaph275/Clean_VCF/ABBA_pop_stats/sel_scan/scaffolds_imputed/Scaffold_1__2_contigs__length_78070536.vcf -B 78070536 2549116 -S /gxfs_work1/cau/suaph275/Clean_VCF/ABBA_pop_stats/sel_scan/2931_pop.txt

RAiSD -n 2932_SC_1 -f -y 2 -M 0 -R -c 1 -w 50 -I /gxfs_work1/cau/suaph275/Clean_VCF/ABBA_pop_stats/sel_scan/scaffolds_imputed/Scaffold_1__2_contigs__length_78070536.vcf -B 78070536 2549116 -S /gxfs_work1/cau/suaph275/Clean_VCF/ABBA_pop_stats/sel_scan/2932_pop.txt

RAiSD -n 4107_SC_1 -f -y 2 -M 0 -R -c 1 -w 50 -I /gxfs_work1/cau/suaph275/Clean_VCF/ABBA_pop_stats/sel_scan/scaffolds_imputed/Scaffold_1__2_contigs__length_78070536.vcf -B 78070536 2549116 -S /gxfs_work1/cau/suaph275/Clean_VCF/ABBA_pop_stats/sel_scan/4107_pop.txt

RAiSD -n Two_C_SC_1 -f -y 2 -M 0 -R -c 1 -w 50 -I /gxfs_work1/cau/suaph275/Clean_VCF/ABBA_pop_stats/sel_scan/scaffolds_imputed/Scaffold_1__2_contigs__length_78070536.vcf -B 78070536 2549116 -S /gxfs_work1/cau/suaph275/Clean_VCF/ABBA_pop_stats/sel_scan/two_central.txt

RAiSD -n Two_SH_SC_1 -f -y 2 -M 0 -R -c 1 -w 50 -I /gxfs_work1/cau/suaph275/Clean_VCF/ABBA_pop_stats/sel_scan/scaffolds_imputed/Scaffold_1__2_contigs__length_78070536.vcf -B 78070536 2549116 -S /gxfs_work1/cau/suaph275/Clean_VCF/ABBA_pop_stats/sel_scan/two_highland.txt
