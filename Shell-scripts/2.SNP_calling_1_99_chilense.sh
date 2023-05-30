####################################################################
# script used to trim, make bam, sam and g.vcf files for 99 Solanum chilense plants. 
# Anik Dutta (02/04/2022)    
####################################################################

# Before starting loop do the following:
#       1 - Genome reference fasta indexed (bt2) [code: bowtie2-build S_chilense_Hirise.fasta S_chilense_Hirise]
#       2 - Genome reference fasta index (.fai) [code: samtools faidx S_chilense_Hirise.fasta]
#       3 - Genome reference fasta dictionary (dict) [code: picard CreateSequenceDictionary -R S_chilense_Hirise.fasta]

### Add Readgroup for running GATK. This was done additionally because I ran into problem after aligning the raw reads to the reference genome using bwa (https://gatk.broadinstitute.org/hc/en-us/articles/360035532352-Errors-about-read-group-RG-information)####

# # 2) load all necessary software
#module load trimmomatic
module load bowtie2
module load bwa
module load samtools
module load picard
module load gatk
module load python

# 3) Set path to softwares and dependencies
#TRIM="/gxfs_home/sw/spack/spack0.16.0/usr/opt/spack/linux-rhel8-x86_64/gcc-10.2.0/trimmomatic-0.39-4mrdo2jaxvlekixczytlpeyezugj4su6/bin/"
#IlluminaAdaptor="/gxfs_home/sw/spack/spack0.16.0/usr/opt/spack/linux-rhel8-x86_64/gcc-10.2.0/trimmomatic-0.39-4mrdo2jaxvlekixczytlpeyezugj4su6/share/adapters/TruSeq3-PE.fa"

## Use the Pipe sign instead of tab or comma. Always start with "0" as the first column because in bash, 0 means the first column. Always put "placeholder" after the last column if you are trying to read the sample names from there. Also, after the last line, put an "enter" otherwise bash does not read the last line. Use cat to see that.

# 4) Set pipeline
IFS='|'
while read line
do
n=24
#echo $line
linearray=( $line )
givenRawR1=${linearray[0]} # first column of the input file
givenRawR2=${linearray[1]} # second column of the input file
givenSampleR1=${linearray[2]} # third column of the input file
givenSampleR2=${linearray[3]} # fourth column of the input file
Sample=${linearray[4]} # fifth column of the input file

RawReadsFolder="/gxfs_work1/cau/suaph275/Chilense_Hiseq"
rawR1=$( find $RawReadsFolder -name "*${givenRawR1}" )
rawR2=$( find $RawReadsFolder -name "*${givenRawR2}" )

#echo $rawR1

#echo $Sample $rawR1 $rawR2

#TrimReadsFolder="/gxfs_work1/cau/suaph275/Chilense_Hiseq/Trim_fastq"
#trimR1="$TrimReadsFolder/fp_$givenSampleR1"
#trimR2="$TrimReadsFolder/rp_$givenSampleR2"
#unpR1="$TrimReadsFolder/fu_$givenSampleR1"
#unpR2="$TrimReadsFolder/ru_$givenSampleR2"

#echo $trimR1 $trimR2

# Assign the input and output folders to short variables for making it easier. 

REF="/gxfs_work1/cau/suaph275/Chilense_Hiseq/S_chilense_Hirise"
REFERENCE="/gxfs_work1/cau/suaph275/Chilense_Hiseq/S_chilense_Hirise.fasta"
SamFolder="/gxfs_work1/cau/suaph275/Chilense_Hiseq/Sam_files"
BamFolder="/gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_files"
BamRGFolder="/gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files"
PicardFolder="/gxfs_work1/cau/suaph275/Chilense_Hiseq/Picard_files"
gVCFFolder="/gxfs_work1/cau/suaph275/Chilense_Hiseq/gvcf_files"
project="Chilense_99"

echo $Sample

echo -e '#!/bin/bash' > $Sample.sh

#### As my raw data were already cleaned, I did not need to use Trimmomatic to remove adapter contamination.

#echo "java -jar $TRIM/trimmomatic-0.39.jar PE -threads $n $rawR1 $rawR2 $trimR1 $unpR1 $trimR2 $unpR2 IILLUMINACLIP:$IlluminaAdaptor:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:5:10 MINLEN:50" >> $Sample.sh

# I previously used bowtie2 for small fungal genomes but I see most of the plant genomics analysis are done using bwa. So, I opted for BWA.

#echo "bowtie2 -p $n --very-sensitive-local --rg-id $Sample --rg SM:$Sample -x $REF -1 $rawR1 -2 $rawR2 -S $SamFolder/$Sample.sam" >> $Sample.sh

# Run bwa to map the raw reads to the reference genome

echo "bwa mem -t 24 -M $REFERENCE $rawR1 $rawR2 -o $SamFolder/$Sample.sam" >> $Sample.sh

# Creats Sam files

echo "samtools view -bS $SamFolder/$Sample.sam | samtools sort -o $BamFolder/$Sample.bam" >> $Sample.sh

# Index sam files

echo "samtools index $BamFolder/$Sample.bam" >> $Sample.sh

# Delete all the sam files because they used up a lot of space.

echo "rm $SamFolder/$Sample.sam" >> $Sample.sh

# Run Picard tools to add a header to each of the bam file for downstream analysis

echo "picard AddOrReplaceReadGroups \
    I=$BamFolder/$Sample.bam \
    O=$BamRGFolder/$Sample.bam \
    SORT_ORDER=coordinate \
    RGID=$Sample \
    RGLB=Chilense_Hiseq \
    RGPL=illumina \
    RGSM=$Sample \
    RGPU= $project \
    CREATE_INDEX=True" >> $Sample.sh

echo "rm $BamFolder/$Sample.bam" >> $Sample.sh    

# Remove duplicated reads using Picard Tools

echo "picard MarkDuplicates INPUT=$BamRGFolder/$Sample.bam OUTPUT=$PicardFolder/$Sample.DuplMark.bam METRICS_FILE=metrics_$isolate.txt CREATE_INDEX=true" >> $Sample.sh

# Run GATK for calling haplotypes

echo "gatk --java-options \"-Xmx4G\" HaplotypeCaller -R "$REF".fasta -ploidy 2 --emit-ref-confidence GVCF -I $PicardFolder/$Sample.DuplMark.bam -O $gVCFFolder/$Sample.g.vcf" >> $Sample.sh

echo "rm $PicardFolder/$Sample.DuplMark.bam" >> $Sample.sh

# Submit the job. Here all the above steps are run separately for each sample.

echo "sbatch --job-name=GATK_SSc --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=$n --qos=long --time=10-00:00:00 --mem=120G --error=job.%J.err --output=job.%J.out --mail-type=FAIL --mail-user=anik.dutta@phytomed.uni-kiel.de --partition=cluster $Sample.sh" >> sbatch_GATK

echo "sleep 0.2" >> sbatch_GATK

done < /gxfs_work1/cau/suaph275/Chilense_Hiseq/Chilense_accessions_2.txt


# submission list
#sh sbatch_GATK
