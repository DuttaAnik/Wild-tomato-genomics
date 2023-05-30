#!/bin/bash

# Anik Dutta

fasttree -gtr -gamma -log log.file -nt Pruned.0.2.fasta > Fast_tree

# Runs the FastTree to create a phylogenetic tree
# -gtr : Use the General Time Reversible model of nucleotide substitution
# -gamma : Use the Gamma20 model of rate variation (which accommodates variable evolutionary rates across sites)
# -nt : Treat the input sequences as nucleotide sequences
# Pruned.0.2.fasta : The input fasta file that contains the pruned sequence data


# Create neighbor joining tree using the pruned data (adopted from https://github.com/jingwanglab/Populus-Genomic-Consequences-of-Hybridization/blob/main/1-population%20structure%20and%20phylogenetic%20analyses/outgroup_merge_vcf.populus227.ld_pruned.sh)
# Convert the plink bfiles that are pruned at r2>0.2

plink --allow-extra-chr --bfile Pruned.0.2 --recode vcf-iid --out Pruned.0.2

# Use plink 'ibs-1' for pairwise distance

plink --vcf Pruned.0.2.vcf --distance '1-ibs' --allow-extra-chr --allow-no-sex --out Neighbor.SNP

# Prepare file for neighbor-joining tree in MEGA X

cut -f 1 Neighbor.SNP.mdist.id | awk 'BEGIN{a="#"}{printf("%s",a);for(i=1;i<=NF;i++){printf($i);printf("")}printf("%s","\n")}' > Neighbor.SNP.mdist.id.mega.txt
echo -e "#mega" > Neighbor.SNP.mdist.meg

# here the Title event was not found. So I used the '\' in front of the '!' to escape the error (https://stackoverflow.com/questions/69475229/echo-displays-an-event-not-found-error)

echo -e "\!Title: Neighbor.SNP.mdist;" >> Neighbor.SNP.mdist.meg
echo -e "" >> Neighbor.SNP.mdist.meg
cat Neighbor.SNP.mdist.id.mega.txt >> Neighbor.SNP.mdist.meg
echo -e "" >> Neighbor.SNP.mdist.meg
cat Neighbor.SNP.mdist >> Neighbor.SNP.mdist.meg 

## Then use the Neighbor.SNP.mdist.meg in MEGAX to run NJ tree. Use the pairwise file opening in MEGA.
## Then export the tree in newick format to edit and visualize in Figtree.


# Convert a vcf file to phylip format for running phylogenetic analysis (optional). Source: https://github.com/edgardomortiz/vcf2phylip

module load python

python /gxfs_home/cau/suaph275/VCF2phylip/vcf2phylip.py -i imputed.vcf

# Run IQ-tree to create a phylogenetic tree (optional)
# -m Model Finder Plus
# -nt Number of CPUs set to AUTO
# -bb is for boostrap
# SH-alrt test
# -bnni to reduce the risk of overestimating branch supports with UFBoot due to severe model violations

iqtree -s Pruned.0.2.min4.phy -m MFP -nt AUTO -bb 1000 -alrt 1000 -bnni -redo
