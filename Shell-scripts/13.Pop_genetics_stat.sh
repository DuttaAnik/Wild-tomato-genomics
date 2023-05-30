#!/bin/bash

# Anik Dutta

# Making files for population genetics stats for each population. The final file has the first column as individual and the second column as the population

declare -a allNames=("2931" "2747" "2932" "1963" "1958" "4330" "3111" "4107" "4117")

for names in "${allNames[@]}"
do
echo $names
grep "$names" Population_ID.txt > $names.txt
done


# Estimating TajimaD and Pi with a window size of 10kb. First, upload those population text files in the vcf file directory and then run the analysis

# window_size=10000

for i in 1958 1963 2747 2931 2932 3111 4107 4117 4330
do
full_pop=""pop."$i"
vcftools --vcf Chilense_99.genotyped.SNP.Clean.MAXMiss.MinDP.recode.vcf --keep "$i".txt --maf 0.001 --window-pi $window_size --out "$full_pop".10kb.pi &
vcftools --vcf Chilense_99.genotyped.SNP.Clean.MAXMiss.MinDP.recode.vcf --keep "$i".txt --maf 0.001 --TajimaD $window_size --out "$full_pop".10kb.TajD &
vcftools --vcf Chilense_99.genotyped.SNP.Clean.MAXMiss.MinDP.recode.vcf --keep "$i".txt --SNPdensity 1000 --out "$full_pop".SNPdensity &
done

## Estimating Observed and expected heterozygosity from the all filtered applied VCF file (output.vcf)

vcftools --vcf output.vcf --het --out hetero

# then in the output file, do this
Observed.het = (N_SITES - O(HOM)) / N_SITES
Expected.het = (N_SITES - E(HOM)) / N_SITES

# Working with singleton outputs

sed 1d 2747.singleton.singletons | cut -f 3,5 | grep -v '^D' | sort | uniq -c > 2747.singletone.count
sed 1d 2931.singleton.singletons | cut -f 3,5 | grep -v '^D' | sort | uniq -c > 2931.singletone.count
sed 1d 2932.singleton.singletons | cut -f 3,5 | grep -v '^D' | sort | uniq -c > 2932.singletone.count
sed 1d 4107.singleton.singletons | cut -f 3,5 | grep -v '^D' | sort | uniq -c > 4107.singletone.count
sed 1d 4117.singleton.singletons | cut -f 3,5 | grep -v '^D' | sort | uniq -c > 4117.singletone.count
sed 1d 4330.singleton.singletons | cut -f 3,5 | grep -v '^D' | sort | uniq -c > 4330.singletone.count

# To convert the VCF file into Simon Martin .geno format for pi, Dxy and Fst estimation, first replace the missing values with the code

cat Chilense_99.genotyped.SNP.Clean.MAXMiss.MinDP.recode.vcf | perl -pe "s/\s\.:/\t.\/.:/g" > Coverted_Miss.vcf

# There was a problem with the parseVCF.py as it was not running because numpy was not found. So, I installed it with the code
 pip3 install -U numpy
# Somehow it did not work. Then, I had to create a conda environment with python 3. But bgzip did not work on py3
conda create --name py3 python=3
conda activate py3
pip3 install -U numpy

cat Chilense_99.genotyped.SNP.Invariant.Clean.vcf.recode.vcf | perl -pe "s/\s\.:/\t.\/.:/g" > Chilense_missingreplaced.vcf

python /gxfs_home/cau/suaph275/genomics_general/VCF_processing/parseVCF.py -i Chilense_missingreplaced.vcf > Chilense.geno

# Actiavte py3 before pixy. --bypass_invariant_check 'yes' \

tabix test_filtered.vcf.gz

pixy --stats pi fst dxy \
--vcf test_filtered.vcf.gz \
--populations Pop.info.txt \
--window_size 10000 \
--n_cores 16 \
--output_prefix Chilense_Pixy_output.out


## ABBA from hrazif (did not work)

# python /gxfs_home/cau/suaph275/genomics_general/freq.py -g Chilense.geno -o Chilense.derFreq.tsv.gz -f phased  -p P_2931 -p P_2747 -p P_1958 -p P_3111 -p P_4117 -p P_4330 -p P_2932 -p P_4107  -p P_1963 --popsFile Pop.info.txt --target derived

# Pi, Dxy, Fst (optional)

# python /gxfs_home/cau/suaph275/genomics_general/popgenWindows.py -g Chilense.geno -o Chilense.Fst.Dxy.pi.csv.gz -f phased -w 100000 -m 5 -s 100000 -T 32 -p P_2931 -p P_2747 -p P_2932 -p P_1963 -p P_1958 -p P_4330 -p P_3111 -p P_4107 -p P_4117 --popsFile Pop.info.txt --writeFailedWindows

#http://evomics.org/learning/population-and-speciation-genomics/2018-population-and-speciation-genomics/abba-baba-statistics/
#https://github.com/simonhmartin/tutorials/tree/master/ABBA_BABA_whole_genome