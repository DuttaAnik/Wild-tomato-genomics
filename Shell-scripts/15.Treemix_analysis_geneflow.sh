#!/bin/bash

# Anik Dutta

# Treemix analysis to test for signatures of gene flow/migration events among 9 S. chilense populations. 

# Source: https://bitbucket.org/nygcresearch/treemix/wiki/Home

# Pop.clust contains three columns tab separated. ind ind pop. Make the input file using LD pruned (r2=0.2) SNP datafile 

plink --allow-extra-chr --bfile Pruned.0.2 --freq --missing --within Pop.clust --allow-no-sex

# Compress the output file

gzip plink.frq.strat

# Activate conda to run python version 2

conda activate py2

# Convert the compressed file to treemix input file format

python plink2treemix.py plink.frq.strat.gz treemix.frq.gz

# Activate the conda environment to run treemix

conda activate treemix

# Run treemix with 10 migration events {0-10}

for i in {0..10}
do 
treemix -i treemix.frq.gz -root P_3111 -k 12000 -m $i -bootstrap -k 500 -noss -se -o Chilense."$i" &
done

# You can run treemix for each migration event separately.

#######################
### zero migration events
#######################

echo "Test zero migration event"
#build tree assuming populations are independent
treemix -i treemix.frq.gz -o Chilense0
#place outgroup
treemix -i treemix.frq.gz -root P_3111 -o Chilense0
#use sliding window 1,000 SNPs to build the tree to account for LD
treemix -i treemix.frq.gz -root P_3111 -k 1000 -o Chilense0
#build the tree with allowed migration events
treemix -i treemix.frq.gz -root P_3111 -k 1000 -m 0 -o Chilense0
#bootstrap the tree
treemix -i treemix.frq.gz -root P_3111 -bootstrap  -m 0 -k 500  -o Chilense0

#######################
### one migration events
#######################
echo "Test one migration event"

#build tree assuming populations are independent
treemix -i treemix.frq.gz -o Chilense1

#place outgroup
treemix -i treemix.frq.gz -root P_3111 -o Chilense1

#use sliding window 1,000 SNPs to build the tree to account for LD
treemix -i treemix.frq.gz -root P_3111 -k 1000 -o Chilense1

#build the tree with allowed migration events
treemix -i treemix.frq.gz -root P_3111 -k 1000 -m 1 -o Chilense1

#bootstrap the tree
treemix -i treemix.frq.gz -root P_3111 -bootstrap  -m 1 -k 500  -o Chilense1


#######################
### two migration events
#######################
echo "Test two migration event"

#build tree assuming populations are independent
treemix -i treemix.frq.gz -o Chilense2

#place outgroup
treemix -i treemix.frq.gz -root P_3111 -o Chilense2

#use sliding window 1,000 SNPs to build the tree to account for LD
treemix -i treemix.frq.gz -root P_3111 -k 1000 -o Chilense2

#build the tree with allowed migration events
treemix -i treemix.frq.gz -root P_3111 -k 1000 -m 2 -o Chilense2

#bootstrap the tree
treemix -i treemix.frq.gz -root P_3111 -bootstrap  -m 2 -k 500  -o Chilense2

#######################
### three migration events
#######################
echo "Test three migration event"

#build tree assuming populations are independent
treemix -i treemix.frq.gz -o Chilense3

#place outgroup
treemix -i treemix.frq.gz -root P_3111 -o Chilense3

#use sliding window 1,000 SNPs to build the tree to account for LD
treemix -i treemix.frq.gz -root P_3111 -k 1000 -o Chilense3

#build the tree with allowed migration events
treemix -i treemix.frq.gz -root P_3111 -k 1000 -m 3 -o Chilense3

#bootstrap the tree
treemix -i treemix.frq.gz -root P_3111 -bootstrap  -m 3 -k 500  -o Chilense3

#######################
### four migration events
#######################
echo "Test four migration event"

#build tree assuming populations are independent
treemix -i treemix.frq.gz -o Chilense4

#place outgroup
treemix -i treemix.frq.gz -root P_3111 -o Chilense4

#use sliding window 1,000 SNPs to build the tree to account for LD
treemix -i treemix.frq.gz -root P_3111 -k 1000 -o Chilense4

#build the tree with allowed migration events
treemix -i treemix.frq.gz -root P_3111 -k 1000 -m 4 -o Chilense4

#bootstrap the tree
treemix -i treemix.frq.gz -root P_3111 -bootstrap  -m 4 -k 500  -o Chilense4

#######################
### five migration events
#######################
echo "Test five migration event"

#build tree assuming populations are independent
treemix -i treemix.frq.gz -o Chilense5

#place outgroup
treemix -i treemix.frq.gz -root P_3111 -o Chilense5

#use sliding window 1,000 SNPs to build the tree to account for LD
treemix -i treemix.frq.gz -root P_3111 -k 1000 -o Chilense5

#build the tree with allowed migration events
treemix -i treemix.frq.gz -root P_3111 -k 1000 -m 5 -o Chilense5

#bootstrap the tree
treemix -i treemix.frq.gz -root P_3111 -bootstrap  -m 5 -k 500  -o Chilense5

#######################
### six migration events
#######################
echo "Test six migration event"

#build tree assuming populations are independent
treemix -i treemix.frq.gz -o Chilense6

#place outgroup
treemix -i treemix.frq.gz -root P_3111 -o Chilense6

#use sliding window 1,000 SNPs to build the tree to account for LD
treemix -i treemix.frq.gz -root P_3111 -k 1000 -o Chilense6

#build the tree with allowed migration events
treemix -i treemix.frq.gz -root P_3111 -k 1000 -m 6 -o Chilense6

#bootstrap the tree
treemix -i treemix.frq.gz -root P_3111 -bootstrap  -m 6 -k 500  -o Chilense6

#######################
### seven migration events
#######################
echo "Test seven migration event"

#build tree assuming populations are independent
treemix -i treemix.frq.gz -o Chilense7

#place outgroup
treemix -i treemix.frq.gz -root P_3111 -o Chilense7

#use sliding window 1,000 SNPs to build the tree to account for LD
treemix -i treemix.frq.gz -root P_3111 -k 1000 -o Chilense7

#build the tree with allowed migration events
treemix -i treemix.frq.gz -root P_3111 -k 1000 -m 7 -o Chilense7

#bootstrap the tree
treemix -i treemix.frq.gz -root P_3111 -bootstrap  -m 7 -k 500  -o Chilense7

#######################
### eight migration events
#######################
echo "Test eight migration event"

#build tree assuming populations are independent
treemix -i treemix.frq.gz -o Chilense8

#place outgroup
treemix -i treemix.frq.gz -root P_3111 -o Chilense8

#use sliding window 1,000 SNPs to build the tree to account for LD
treemix -i treemix.frq.gz -root P_3111 -k 1000 -o Chilense8

#build the tree with allowed migration events
treemix -i treemix.frq.gz -root P_3111 -k 1000 -m 8 -o Chilense8

#bootstrap the tree
treemix -i treemix.frq.gz -root P_3111 -bootstrap  -m 8 -k 500  -o Chilense8

#######################
### nine migration events
#######################
echo "Test nine migration event"

#build tree assuming populations are independent
treemix -i treemix.frq.gz -o Chilense9

#place outgroup
treemix -i treemix.frq.gz -root P_3111 -o Chilense9

#use sliding window 1,000 SNPs to build the tree to account for LD
treemix -i treemix.frq.gz -root P_3111 -k 1000 -o Chilense9

#build the tree with allowed migration events
treemix -i treemix.frq.gz -root P_3111 -k 1000 -m 9 -o Chilense9

#bootstrap the tree
treemix -i treemix.frq.gz -root P_3111 -bootstrap  -m 9 -k 500  -o Chilense9

#######################
### ten migration events
#######################
echo "Test ten migration event"

#build tree assuming populations are independent
treemix -i treemix.frq.gz -o Chilense10

#place outgroup
treemix -i treemix.frq.gz -root P_3111 -o Chilense10

#use sliding window 1,000 SNPs to build the tree to account for LD
treemix -i treemix.frq.gz -root P_3111 -k 1000 -o Chilense10

#build the tree with allowed migration events
treemix -i treemix.frq.gz -root P_3111 -k 1000 -m 10 -o Chilense10

#bootstrap the tree
treemix -i treemix.frq.gz -root P_3111 -bootstrap  -m 10 -k 500  -o Chilense10


# Re-run treemix to find the optimum number of migration that supposedly took place to find out with optm package.

for m in {1..10}; do
    for i in {1..5}; do
         Generate random seed
        s=$RANDOM
        echo "Random seed = ${s}"
        treemix -i treemix.frq.gz -m ${m} -global -bootstrap -k 1000 -root P_1963 -se -noss -seed ${s} -o Chilense.${i}.${m}
    done
done

# Four population test to quantify the amount of geneflow. This test does not require an outgroup.

fourpop -i treemix.frq.gz -k 1000 > f4_out.txt
