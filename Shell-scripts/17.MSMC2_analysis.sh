#!/bin/bash

# Demographic analysis using MSMC2. Source: https://github.com/stschiff/msmc2
# This is a kind of redundant analysis since Kai et al. (2022, New Phytologist) already did this thorough analysis using the same populations.

# This is an updated script for generating vcf file according the depth of each scaffold estimated from samtools (11/10/2022) 

## 2931_pop

for Scaffold in `cat data_chr.txt` 
do 
echo "Handling scaffold $Scaffold" 
MEANCOV=`samtools depth -r $Scaffold /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C004.bam | awk '{sum += $3} END {if (NR==0) print NR; else print sum / NR}' | tr ',' '.'` 
echo ${Scaffold} $MEANCOV >> C004.txt 
echo "Mean coverage for this individual, scaffold ${Scaffold}: $MEANCOV"
samtools mpileup -B -q 20 -Q 20 -C 50 -g -r $Scaffold -f /gxfs_work1/cau/suaph275/Demography/S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C004.bam | bcftools call -c -V indels | python /gxfs_home/cau/suaph275/msmc-tools/bamCaller.py $MEANCOV $Scaffold.C004.mask.bed.gz | bgzip -c > $Scaffold.C004.vcf.gz
done

for Scaffold in `cat data_chr.txt` 
do 
echo "Handling scaffold $Scaffold" 
MEANCOV=`samtools depth -r $Scaffold /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C006.bam | awk '{sum += $3} END {if (NR==0) print NR; else print sum / NR}' | tr ',' '.'` 
echo ${Scaffold} $MEANCOV >> C006.txt 
echo "Mean coverage for this individual, scaffold ${Scaffold}: $MEANCOV"
samtools mpileup -B -q 20 -Q 20 -C 50 -g -r $Scaffold -f /gxfs_work1/cau/suaph275/Demography/S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C006.bam | bcftools call -c -V indels | python /gxfs_home/cau/suaph275/msmc-tools/bamCaller.py $MEANCOV $Scaffold.C006.mask.bed.gz | bgzip -c > $Scaffold.C006.vcf.gz
done

for Scaffold in `cat data_chr.txt` 
do 
echo "Handling scaffold $Scaffold" 
MEANCOV=`samtools depth -r $Scaffold /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C008.bam | awk '{sum += $3} END {if (NR==0) print NR; else print sum / NR}' | tr ',' '.'` 
echo ${Scaffold} $MEANCOV >> C008.txt 
echo "Mean coverage for this individual, scaffold ${Scaffold}: $MEANCOV"
samtools mpileup -B -q 20 -Q 20 -C 50 -g -r $Scaffold -f /gxfs_work1/cau/suaph275/Demography/S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C008.bam | bcftools call -c -V indels | python /gxfs_home/cau/suaph275/msmc-tools/bamCaller.py $MEANCOV $Scaffold.C008.mask.bed.gz | bgzip -c > $Scaffold.C008.vcf.gz
done

######################

## 2747_pop

for Scaffold in `cat data_chr.txt` 
do 
echo "Handling scaffold $Scaffold" 
MEANCOV=`samtools depth -r $Scaffold /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C019.bam | awk '{sum += $3} END {if (NR==0) print NR; else print sum / NR}' | tr ',' '.'` 
echo ${Scaffold} $MEANCOV >> C019.txt 
echo "Mean coverage for this individual, scaffold ${Scaffold}: $MEANCOV"
samtools mpileup -B -q 20 -Q 20 -C 50 -g -r $Scaffold -f /gxfs_work1/cau/suaph275/Demography/S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C019.bam | bcftools call -c -V indels | python /gxfs_home/cau/suaph275/msmc-tools/bamCaller.py $MEANCOV $Scaffold.C019.mask.bed.gz | bgzip -c > $Scaffold.C019.vcf.gz
done

for Scaffold in `cat data_chr.txt` 
do 
echo "Handling scaffold $Scaffold" 
MEANCOV=`samtools depth -r $Scaffold /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C021.bam | awk '{sum += $3} END {if (NR==0) print NR; else print sum / NR}' | tr ',' '.'` 
echo ${Scaffold} $MEANCOV >> C021.txt 
echo "Mean coverage for this individual, scaffold ${Scaffold}: $MEANCOV"
samtools mpileup -B -q 20 -Q 20 -C 50 -g -r $Scaffold -f /gxfs_work1/cau/suaph275/Demography/S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C021.bam | bcftools call -c -V indels | python /gxfs_home/cau/suaph275/msmc-tools/bamCaller.py $MEANCOV $Scaffold.C021.mask.bed.gz | bgzip -c > $Scaffold.C021.vcf.gz
done

for Scaffold in `cat data_chr.txt` 
do 
echo "Handling scaffold $Scaffold" 
MEANCOV=`samtools depth -r $Scaffold /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C025.bam | awk '{sum += $3} END {if (NR==0) print NR; else print sum / NR}' | tr ',' '.'` 
echo ${Scaffold} $MEANCOV >> C025.txt 
echo "Mean coverage for this individual, scaffold ${Scaffold}: $MEANCOV"
samtools mpileup -B -q 20 -Q 20 -C 50 -g -r $Scaffold -f /gxfs_work1/cau/suaph275/Demography/S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C025.bam | bcftools call -c -V indels | python /gxfs_home/cau/suaph275/msmc-tools/bamCaller.py $MEANCOV $Scaffold.C025.mask.bed.gz | bgzip -c > $Scaffold.C025.vcf.gz
done

######################

## 2932_pop

for Scaffold in `cat data_chr.txt` 
do 
echo "Handling scaffold $Scaffold" 
MEANCOV=`samtools depth -r $Scaffold /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C033.bam | awk '{sum += $3} END {if (NR==0) print NR; else print sum / NR}' | tr ',' '.'` 
echo ${Scaffold} $MEANCOV >> C033.txt 
echo "Mean coverage for this individual, scaffold ${Scaffold}: $MEANCOV"
samtools mpileup -B -q 20 -Q 20 -C 50 -g -r $Scaffold -f /gxfs_work1/cau/suaph275/Demography/S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C033.bam | bcftools call -c -V indels | python /gxfs_home/cau/suaph275/msmc-tools/bamCaller.py $MEANCOV $Scaffold.C033.mask.bed.gz | bgzip -c > $Scaffold.C033.vcf.gz
done

for Scaffold in `cat data_chr.txt` 
do 
echo "Handling scaffold $Scaffold" 
MEANCOV=`samtools depth -r $Scaffold /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C035.bam | awk '{sum += $3} END {if (NR==0) print NR; else print sum / NR}' | tr ',' '.'` 
echo ${Scaffold} $MEANCOV >> C035.txt 
echo "Mean coverage for this individual, scaffold ${Scaffold}: $MEANCOV"
samtools mpileup -B -q 20 -Q 20 -C 50 -g -r $Scaffold -f /gxfs_work1/cau/suaph275/Demography/S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C035.bam | bcftools call -c -V indels | python /gxfs_home/cau/suaph275/msmc-tools/bamCaller.py $MEANCOV $Scaffold.C035.mask.bed.gz | bgzip -c > $Scaffold.C035.vcf.gz
done

for Scaffold in `cat data_chr.txt` 
do 
echo "Handling scaffold $Scaffold" 
MEANCOV=`samtools depth -r $Scaffold /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C037.bam | awk '{sum += $3} END {if (NR==0) print NR; else print sum / NR}' | tr ',' '.'` 
echo ${Scaffold} $MEANCOV >> C037.txt 
echo "Mean coverage for this individual, scaffold ${Scaffold}: $MEANCOV"
samtools mpileup -B -q 20 -Q 20 -C 50 -g -r $Scaffold -f /gxfs_work1/cau/suaph275/Demography/S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C037.bam | bcftools call -c -V indels | python /gxfs_home/cau/suaph275/msmc-tools/bamCaller.py $MEANCOV $Scaffold.C037.mask.bed.gz | bgzip -c > $Scaffold.C037.vcf.gz
done
######################

## 1963_pop
for Scaffold in `cat data_chr.txt` 
do 
echo "Handling scaffold $Scaffold" 
MEANCOV=`samtools depth -r $Scaffold /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C047.bam | awk '{sum += $3} END {if (NR==0) print NR; else print sum / NR}' | tr ',' '.'` 
echo ${Scaffold} $MEANCOV >> C047.txt 
echo "Mean coverage for this individual, scaffold ${Scaffold}: $MEANCOV"
samtools mpileup -B -q 20 -Q 20 -C 50 -g -r $Scaffold -f /gxfs_work1/cau/suaph275/Demography/S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C047.bam | bcftools call -c -V indels | python /gxfs_home/cau/suaph275/msmc-tools/bamCaller.py $MEANCOV $Scaffold.C047.mask.bed.gz | bgzip -c > $Scaffold.C047.vcf.gz
done

for Scaffold in `cat data_chr.txt` 
do 
echo "Handling scaffold $Scaffold" 
MEANCOV=`samtools depth -r $Scaffold /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C049.bam | awk '{sum += $3} END {if (NR==0) print NR; else print sum / NR}' | tr ',' '.'` 
echo ${Scaffold} $MEANCOV >> C049.txt 
echo "Mean coverage for this individual, scaffold ${Scaffold}: $MEANCOV"
samtools mpileup -B -q 20 -Q 20 -C 50 -g -r $Scaffold -f /gxfs_work1/cau/suaph275/Demography/S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C049.bam | bcftools call -c -V indels | python /gxfs_home/cau/suaph275/msmc-tools/bamCaller.py $MEANCOV $Scaffold.C049.mask.bed.gz | bgzip -c > $Scaffold.C049.vcf.gz
done

for Scaffold in `cat data_chr.txt` 
do 
echo "Handling scaffold $Scaffold" 
MEANCOV=`samtools depth -r $Scaffold /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C053.bam | awk '{sum += $3} END {if (NR==0) print NR; else print sum / NR}' | tr ',' '.'` 
echo ${Scaffold} $MEANCOV >> C053.txt 
echo "Mean coverage for this individual, scaffold ${Scaffold}: $MEANCOV"
samtools mpileup -B -q 20 -Q 20 -C 50 -g -r $Scaffold -f /gxfs_work1/cau/suaph275/Demography/S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C053.bam | bcftools call -c -V indels | python /gxfs_home/cau/suaph275/msmc-tools/bamCaller.py $MEANCOV $Scaffold.C053.mask.bed.gz | bgzip -c > $Scaffold.C053.vcf.gz
done
######################

## 1958_pop

for Scaffold in `cat data_chr.txt` 
do 
echo "Handling scaffold $Scaffold" 
MEANCOV=`samtools depth -r $Scaffold /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C061.bam | awk '{sum += $3} END {if (NR==0) print NR; else print sum / NR}' | tr ',' '.'` 
echo ${Scaffold} $MEANCOV >> C061.txt 
echo "Mean coverage for this individual, scaffold ${Scaffold}: $MEANCOV"
samtools mpileup -B -q 20 -Q 20 -C 50 -g -r $Scaffold -f /gxfs_work1/cau/suaph275/Demography/S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C061.bam | bcftools call -c -V indels | python /gxfs_home/cau/suaph275/msmc-tools/bamCaller.py $MEANCOV $Scaffold.C061.mask.bed.gz | bgzip -c > $Scaffold.C061.vcf.gz
done

for Scaffold in `cat data_chr.txt` 
do 
echo "Handling scaffold $Scaffold" 
MEANCOV=`samtools depth -r $Scaffold /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C063.bam | awk '{sum += $3} END {if (NR==0) print NR; else print sum / NR}' | tr ',' '.'` 
echo ${Scaffold} $MEANCOV >> C063.txt 
echo "Mean coverage for this individual, scaffold ${Scaffold}: $MEANCOV"
samtools mpileup -B -q 20 -Q 20 -C 50 -g -r $Scaffold -f /gxfs_work1/cau/suaph275/Demography/S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C063.bam | bcftools call -c -V indels | python /gxfs_home/cau/suaph275/msmc-tools/bamCaller.py $MEANCOV $Scaffold.C063.mask.bed.gz | bgzip -c > $Scaffold.C063.vcf.gz
done

for Scaffold in `cat data_chr.txt` 
do 
echo "Handling scaffold $Scaffold" 
MEANCOV=`samtools depth -r $Scaffold /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C067.bam | awk '{sum += $3} END {if (NR==0) print NR; else print sum / NR}' | tr ',' '.'` 
echo ${Scaffold} $MEANCOV >> C067.txt 
echo "Mean coverage for this individual, scaffold ${Scaffold}: $MEANCOV"
samtools mpileup -B -q 20 -Q 20 -C 50 -g -r $Scaffold -f /gxfs_work1/cau/suaph275/Demography/S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C067.bam | bcftools call -c -V indels | python /gxfs_home/cau/suaph275/msmc-tools/bamCaller.py $MEANCOV $Scaffold.C067.mask.bed.gz | bgzip -c > $Scaffold.C067.vcf.gz
done
######################

## 4330_pop

for Scaffold in `cat data_chr.txt` 
do 
echo "Handling scaffold $Scaffold" 
MEANCOV=`samtools depth -r $Scaffold /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C077.bam | awk '{sum += $3} END {if (NR==0) print NR; else print sum / NR}' | tr ',' '.'` 
echo ${Scaffold} $MEANCOV >> C077.txt 
echo "Mean coverage for this individual, scaffold ${Scaffold}: $MEANCOV"
samtools mpileup -B -q 20 -Q 20 -C 50 -g -r $Scaffold -f /gxfs_work1/cau/suaph275/Demography/S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C077.bam | bcftools call -c -V indels | python /gxfs_home/cau/suaph275/msmc-tools/bamCaller.py $MEANCOV $Scaffold.C077.mask.bed.gz | bgzip -c > $Scaffold.C077.vcf.gz
done

for Scaffold in `cat data_chr.txt` 
do 
echo "Handling scaffold $Scaffold" 
MEANCOV=`samtools depth -r $Scaffold /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C079.bam | awk '{sum += $3} END {if (NR==0) print NR; else print sum / NR}' | tr ',' '.'` 
echo ${Scaffold} $MEANCOV >> C079.txt 
echo "Mean coverage for this individual, scaffold ${Scaffold}: $MEANCOV"
samtools mpileup -B -q 20 -Q 20 -C 50 -g -r $Scaffold -f /gxfs_work1/cau/suaph275/Demography/S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C079.bam | bcftools call -c -V indels | python /gxfs_home/cau/suaph275/msmc-tools/bamCaller.py $MEANCOV $Scaffold.C079.mask.bed.gz | bgzip -c > $Scaffold.C079.vcf.gz
done

for Scaffold in `cat data_chr.txt` 
do 
echo "Handling scaffold $Scaffold" 
MEANCOV=`samtools depth -r $Scaffold /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C082.bam | awk '{sum += $3} END {if (NR==0) print NR; else print sum / NR}' | tr ',' '.'` 
echo ${Scaffold} $MEANCOV >> C082.txt 
echo "Mean coverage for this individual, scaffold ${Scaffold}: $MEANCOV"
samtools mpileup -B -q 20 -Q 20 -C 50 -g -r $Scaffold -f /gxfs_work1/cau/suaph275/Demography/S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C082.bam | bcftools call -c -V indels | python /gxfs_home/cau/suaph275/msmc-tools/bamCaller.py $MEANCOV $Scaffold.C082.mask.bed.gz | bgzip -c > $Scaffold.C082.vcf.gz
done
######################
## 3111_pop

for Scaffold in `cat data_chr.txt` 
do 
echo "Handling scaffold $Scaffold" 
MEANCOV=`samtools depth -r $Scaffold /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C090.bam | awk '{sum += $3} END {if (NR==0) print NR; else print sum / NR}' | tr ',' '.'` 
echo ${Scaffold} $MEANCOV >> C090.txt 
echo "Mean coverage for this individual, scaffold ${Scaffold}: $MEANCOV"
samtools mpileup -B -q 20 -Q 20 -C 50 -g -r $Scaffold -f /gxfs_work1/cau/suaph275/Demography/S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C090.bam | bcftools call -c -V indels | python /gxfs_home/cau/suaph275/msmc-tools/bamCaller.py $MEANCOV $Scaffold.C090.mask.bed.gz | bgzip -c > $Scaffold.C090.vcf.gz
done

for Scaffold in `cat data_chr.txt` 
do 
echo "Handling scaffold $Scaffold" 
MEANCOV=`samtools depth -r $Scaffold /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C092.bam | awk '{sum += $3} END {if (NR==0) print NR; else print sum / NR}' | tr ',' '.'` 
echo ${Scaffold} $MEANCOV >> C092.txt 
echo "Mean coverage for this individual, scaffold ${Scaffold}: $MEANCOV"
samtools mpileup -B -q 20 -Q 20 -C 50 -g -r $Scaffold -f /gxfs_work1/cau/suaph275/Demography/S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C092.bam | bcftools call -c -V indels | python /gxfs_home/cau/suaph275/msmc-tools/bamCaller.py $MEANCOV $Scaffold.C092.mask.bed.gz | bgzip -c > $Scaffold.C092.vcf.gz
done

for Scaffold in `cat data_chr.txt` 
do 
echo "Handling scaffold $Scaffold" 
MEANCOV=`samtools depth -r $Scaffold /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C094.bam | awk '{sum += $3} END {if (NR==0) print NR; else print sum / NR}' | tr ',' '.'` 
echo ${Scaffold} $MEANCOV >> C094.txt 
echo "Mean coverage for this individual, scaffold ${Scaffold}: $MEANCOV"
samtools mpileup -B -q 20 -Q 20 -C 50 -g -r $Scaffold -f /gxfs_work1/cau/suaph275/Demography/S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C094.bam | bcftools call -c -V indels | python /gxfs_home/cau/suaph275/msmc-tools/bamCaller.py $MEANCOV $Scaffold.C094.mask.bed.gz | bgzip -c > $Scaffold.C094.vcf.gz
done
######################
## 4107_pop

for Scaffold in `cat data_chr.txt` 
do 
echo "Handling scaffold $Scaffold" 
MEANCOV=`samtools depth -r $Scaffold /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C104.bam | awk '{sum += $3} END {if (NR==0) print NR; else print sum / NR}' | tr ',' '.'` 
echo ${Scaffold} $MEANCOV >> C104.txt 
echo "Mean coverage for this individual, scaffold ${Scaffold}: $MEANCOV"
samtools mpileup -B -q 20 -Q 20 -C 50 -g -r $Scaffold -f /gxfs_work1/cau/suaph275/Demography/S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C104.bam | bcftools call -c -V indels | python /gxfs_home/cau/suaph275/msmc-tools/bamCaller.py $MEANCOV $Scaffold.C104.mask.bed.gz | bgzip -c > $Scaffold.C104.vcf.gz
done

for Scaffold in `cat data_chr.txt` 
do 
echo "Handling scaffold $Scaffold" 
MEANCOV=`samtools depth -r $Scaffold /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C106.bam | awk '{sum += $3} END {if (NR==0) print NR; else print sum / NR}' | tr ',' '.'` 
echo ${Scaffold} $MEANCOV >> C106.txt 
echo "Mean coverage for this individual, scaffold ${Scaffold}: $MEANCOV"
samtools mpileup -B -q 20 -Q 20 -C 50 -g -r $Scaffold -f /gxfs_work1/cau/suaph275/Demography/S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C106.bam | bcftools call -c -V indels | python /gxfs_home/cau/suaph275/msmc-tools/bamCaller.py $MEANCOV $Scaffold.C106.mask.bed.gz | bgzip -c > $Scaffold.C106.vcf.gz
done

for Scaffold in `cat data_chr.txt` 
do 
echo "Handling scaffold $Scaffold" 
MEANCOV=`samtools depth -r $Scaffold /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C108.bam | awk '{sum += $3} END {if (NR==0) print NR; else print sum / NR}' | tr ',' '.'` 
echo ${Scaffold} $MEANCOV >> C108.txt 
echo "Mean coverage for this individual, scaffold ${Scaffold}: $MEANCOV"
samtools mpileup -B -q 20 -Q 20 -C 50 -g -r $Scaffold -f /gxfs_work1/cau/suaph275/Demography/S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C108.bam | bcftools call -c -V indels | python /gxfs_home/cau/suaph275/msmc-tools/bamCaller.py $MEANCOV $Scaffold.C108.mask.bed.gz | bgzip -c > $Scaffold.C108.vcf.gz
done

######################
## 4117_pop

for Scaffold in `cat data_chr.txt` 
do 
echo "Handling scaffold $Scaffold" 
MEANCOV=`samtools depth -r $Scaffold /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C118.bam | awk '{sum += $3} END {if (NR==0) print NR; else print sum / NR}' | tr ',' '.'` 
echo ${Scaffold} $MEANCOV >> C118.txt 
echo "Mean coverage for this individual, scaffold ${Scaffold}: $MEANCOV"
samtools mpileup -B -q 20 -Q 20 -C 50 -g -r $Scaffold -f /gxfs_work1/cau/suaph275/Demography/S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C118.bam | bcftools call -c -V indels | python /gxfs_home/cau/suaph275/msmc-tools/bamCaller.py $MEANCOV $Scaffold.C118.mask.bed.gz | bgzip -c > $Scaffold.C118.vcf.gz
done

for Scaffold in `cat data_chr.txt` 
do 
echo "Handling scaffold $Scaffold" 
MEANCOV=`samtools depth -r $Scaffold /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C121.bam | awk '{sum += $3} END {if (NR==0) print NR; else print sum / NR}' | tr ',' '.'` 
echo ${Scaffold} $MEANCOV >> C121.txt 
echo "Mean coverage for this individual, scaffold ${Scaffold}: $MEANCOV"
samtools mpileup -B -q 20 -Q 20 -C 50 -g -r $Scaffold -f /gxfs_work1/cau/suaph275/Demography/S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C121.bam | bcftools call -c -V indels | python /gxfs_home/cau/suaph275/msmc-tools/bamCaller.py $MEANCOV $Scaffold.C121.mask.bed.gz | bgzip -c > $Scaffold.C121.vcf.gz
done

for Scaffold in `cat data_chr.txt` 
do 
echo "Handling scaffold $Scaffold" 
MEANCOV=`samtools depth -r $Scaffold /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C123.bam | awk '{sum += $3} END {if (NR==0) print NR; else print sum / NR}' | tr ',' '.'` 
echo ${Scaffold} $MEANCOV >> C123.txt 
echo "Mean coverage for this individual, scaffold ${Scaffold}: $MEANCOV"
samtools mpileup -B -q 20 -Q 20 -C 50 -g -r $Scaffold -f /gxfs_work1/cau/suaph275/Demography/S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C123.bam | bcftools call -c -V indels | python /gxfs_home/cau/suaph275/msmc-tools/bamCaller.py $MEANCOV $Scaffold.C123.mask.bed.gz | bgzip -c > $Scaffold.C123.vcf.gz
done

# Generate input files for running MSMC2

dir=/gxfs_work1/cau/suaph275/MSMC_analysis

cat /gxfs_work1/cau/suaph275/MSMC_analysis/data_chr.txt | while read Scaffold 
do 
python /gxfs_home/cau/suaph275/msmc-tools/generate_multihetsep.py --chr ${Scaffold} --mask $dir/1958_pop/${Scaffold}.C061.mask.bed.gz --mask $dir/1958_pop/${Scaffold}.C063.mask.bed.gz --mask $dir/1958_pop/${Scaffold}.C067.mask.bed.gz --mask $dir/1963_pop/${Scaffold}.C047.mask.bed.gz --mask $dir/1963_pop/${Scaffold}.C049.mask.bed.gz --mask $dir/1963_pop/${Scaffold}.C053.mask.bed.gz --mask $dir/2747_pop/${Scaffold}.C019.mask.bed.gz --mask $dir/2747_pop/${Scaffold}.C021.mask.bed.gz --mask $dir/2747_pop/${Scaffold}.C025.mask.bed.gz --mask $dir/2931_pop/${Scaffold}.C004.mask.bed.gz --mask $dir/2931_pop/${Scaffold}.C006.mask.bed.gz --mask $dir/2931_pop/${Scaffold}.C008.mask.bed.gz --mask $dir/2932_pop/${Scaffold}.C033.mask.bed.gz --mask $dir/2932_pop/${Scaffold}.C035.mask.bed.gz --mask $dir/2932_pop/${Scaffold}.C037.mask.bed.gz --mask $dir/3111_pop/${Scaffold}.C090.mask.bed.gz --mask $dir/3111_pop/${Scaffold}.C092.mask.bed.gz --mask $dir/3111_pop/${Scaffold}.C094.mask.bed.gz --mask $dir/4107_pop/${Scaffold}.C104.mask.bed.gz --mask $dir/4107_pop/${Scaffold}.C106.mask.bed.gz --mask $dir/4107_pop/${Scaffold}.C108.mask.bed.gz --mask $dir/4117_pop/${Scaffold}.C118.mask.bed.gz --mask $dir/4117_pop/${Scaffold}.C121.mask.bed.gz --mask $dir/4117_pop/${Scaffold}.C123.mask.bed.gz --mask $dir/4330_pop/${Scaffold}.C077.mask.bed.gz --mask $dir/4330_pop/${Scaffold}.C079.mask.bed.gz --mask $dir/4330_pop/${Scaffold}.C082.mask.bed.gz $dir/1958_pop/${Scaffold}.C061.phased.vcf.gz $dir/1958_pop/${Scaffold}.C063.phased.vcf.gz $dir/1958_pop/${Scaffold}.C067.phased.vcf.gz $dir/1963_pop/${Scaffold}.C047.phased.vcf.gz $dir/1963_pop/${Scaffold}.C049.phased.vcf.gz $dir/1963_pop/${Scaffold}.C053.phased.vcf.gz $dir/2747_pop/${Scaffold}.C019.phased.vcf.gz $dir/2747_pop/${Scaffold}.C021.phased.vcf.gz $dir/2747_pop/${Scaffold}.C025.phased.vcf.gz $dir/2931_pop/${Scaffold}.C004.phased.vcf.gz $dir/2931_pop/${Scaffold}.C006.phased.vcf.gz $dir/2931_pop/${Scaffold}.C008.phased.vcf.gz $dir/2932_pop/${Scaffold}.C033.phased.vcf.gz $dir/2932_pop/${Scaffold}.C035.phased.vcf.gz $dir/2932_pop/${Scaffold}.C037.phased.vcf.gz $dir/3111_pop/${Scaffold}.C090.phased.vcf.gz $dir/3111_pop/${Scaffold}.C092.phased.vcf.gz $dir/3111_pop/${Scaffold}.C094.phased.vcf.gz $dir/4107_pop/${Scaffold}.C104.phased.vcf.gz $dir/4107_pop/${Scaffold}.C106.phased.vcf.gz $dir/4107_pop/${Scaffold}.C108.phased.vcf.gz $dir/4117_pop/${Scaffold}.C118.phased.vcf.gz $dir/4117_pop/${Scaffold}.C121.phased.vcf.gz $dir/4117_pop/${Scaffold}.C123.phased.vcf.gz $dir/4330_pop/${Scaffold}.C077.phased.vcf.gz $dir/4330_pop/${Scaffold}.C079.phased.vcf.gz $dir/4330_pop/${Scaffold}.C082.phased.vcf.gz > $dir/Combined_input/MSMC_input_combined.${Scaffold}.txt 
done

# Indices: [1958: 0,1,2,3,4,5] [1963: 6,7,8,9,10,11] [2747: 12,13,14,15,16,17] [2931: 18,19,20,21,22,23] [2932: 24,25,26,27,28,29] [3111: 30,31,32,33,34,35] [4107: 36,37,38,39,40,41] [4117: 42,43,44,45,46,47] [4330: 48,49,50,51,52,53]

# Run MSMC2 for all individual within one population (an example below)

msmc2_linux64bit -t 12 -o 4330_3 MSMC_input.Scaffold_*.txt


