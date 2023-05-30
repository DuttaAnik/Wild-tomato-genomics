# Population genomic analyses in Solanum chilense

In this project, we aimed to identify genes that are involved in wild plant adaptation to different geographical niches such as desert, mountainous and coastal ranges. We applied a population genomics approach to quantify the level of genetic diversity, demographic patterns and genomic regions under positive selection in nine distinct S. chilense populations.

## What is in the two folders?

The first folder **Accompanying_text_files** contains all the text files used in different analyses and make figures. The other two folders contain all the bash and R scripts used for analysing and visualizing Illumina Whole genome sequence data from 99 wild tomato plants.

### Analyses 

The following analyses were performed with the sequence data

* **FastQC** - *Quality check of raw sequence data*
* **Single Nucleotide Polymorphism calling** - *SNP calling using GATK and PacBio reference genome scaffolded with Hi-C*
* **GATK quality check** - *Quality check of called SNPs bash one liners*
* **SNP annotation** - *SNP annotation with SnpEff*
* **Breadth and Coverage of sequence data** - *Estimating average coverage and breadth of each BAM file after aligning to the reference genome*
* **Get various statistics** - *Extract different statistics like heterozygosity, number of singletons, length of each scaffold etc.*
* **Admixture** - *Admixture analysis using linkage pruned SNP data to identify genetic clusters*
* **Site frequency spectrum and theta** - *SAF and theta estimation using ANGSD*
* **LD decay** - *Linkage-disequilibrium decay analysis in each population*
* **Phylogenetics** - *Maximum likelihood tree, neighbour joining tree, Splitstree analysis*
* **Population genetics** - *Estimating nucleotide diversity (pi), TajimaD, pairwise Fst and Dxy (absolute measure of divergence)*
* **Treemix** - *Treemix analysis to evaluate gene flow and migration*
* **Demographic analysis** - *Demographic analysis to see patterns of effective population changes using PSMC and MSMC2*
* **Selection scan** - *Selection scan analysis using SweeD and RAiSD*
* **Finding genes with Bedtools** - *Extract gene information based on a gene annotation file and significant genomic regions under selection*
* **GWAS** - *A demo script to perform GWAS using flowering time in S. chilense*

## Tools used

* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [GATK](https://gatk.broadinstitute.org/hc/en-us) 
* [bwa](https://github.com/lh3/bwa) 
* [picard](https://github.com/broadinstitute/picard)
* [samtools](https://github.com/samtools/samtools) 
* [vcftools](https://vcftools.github.io/index.html)
* [bcftools](http://samtools.github.io/bcftools/)
* [plink](https://www.cog-genomics.org/plink/) 
* [Admixture](https://dalexander.github.io/admixture/download.html)
* [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD)
* [fasttree](http://www.microbesonline.org/fasttree/) 
* [PGDSpider](http://www.cmpg.unibe.ch/software/PGDSpider/)
* [PopLDdecay](https://github.com/BGI-shenzhen/PopLDdecay) 
* [Pixy](https://pixy.readthedocs.io/en/latest/)
* [SnpEff](https://pcingola.github.io/SnpEff/) 
* [Treemix](https://bitbucket.org/nygcresearch/treemix/wiki/Home)
* [PSMC](https://github.com/lh3/psmc) 
* [MSMC2](https://github.com/stschiff/msmc2)
* [RAiSD](https://github.com/alachins/raisd) 
* [SweeD](https://app.assembla.com/spaces/sweed/git/source)
* [Bedtools](https://bedtools.readthedocs.io/en/latest/) 
* [vcf2gwas](https://github.com/frankvogt/vcf2gwas)
* [R](https://posit.co/download/rstudio-desktop/)

## Author

* **Anik Dutta**
