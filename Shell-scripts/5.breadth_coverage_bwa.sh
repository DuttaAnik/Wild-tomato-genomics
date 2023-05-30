#!/bin/bash

# Anik Dutta

# Write the bash script to the breadth_coverage.sh.
# Read each line from the bam_list.txt file that contains the name of a BAM file. 
# Append a command to the breadth_coverage.sh script that uses samtools to create a pileup for the BAM file.

echo -e '#!/bin/bash' > breadth_coverage.sh
while read line
do
  echo $line 
  echo "samtools mpileup /gxfs_work1/cau/suaph275/Bowtie2_results/Bam_files/$line | awk -v X="${MIN_COVERAGE_DEPTH}" '$4>=X' | wc -l > breadth_$line.txt" >> breadth_coverage.sh
done < bam_list.txt

# This script calculates the breadth of coverage of sequence data in .bam files.
# The breadth of coverage is a measure of the proportion of a reference sequence 
# that is covered by at least a certain number of reads (in this case, 5).

# The -v option is used to assign a value to a variable (X=5). '$4>=X' is a condition
# that checks if the 4th column value (which is the read depth) is greater or equal to 5.

# 'wc -l' counts the number of lines that pass the condition, effectively 
# giving us the number of bases with a coverage of 5 or more. 

samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C004.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C004.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C005.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C005.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C006.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C006.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C007.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C007.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C008.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C008.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C009.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C009.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C012.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C012.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C013.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C013.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C014.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C014.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C015.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C015.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C018.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C018.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C019.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C019.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C020.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C020.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C021.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C021.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C023.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C023.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C025.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C025.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C026.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C026.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C027.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C027.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C028.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C028.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C029.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C029.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C030.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C030.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C032.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C032.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C033.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C033.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C034.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C034.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C035.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C035.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C036.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C036.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C037.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C037.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C040.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C040.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C041.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C041.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C042.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C042.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C044.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C044.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C045.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C045.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C046.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C046.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C047.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C047.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C048.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C048.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C049.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C049.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C052.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C052.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C053.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C053.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C054.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C054.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C056.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C056.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C057.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C057.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C058.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C058.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C059.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C059.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C060.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C060.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C061.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C061.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C062.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C062.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C063.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C063.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C064.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C064.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C067.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C067.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C069.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C069.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C070.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C070.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C071.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C071.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C072.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C072.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C073.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C073.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C074.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C074.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C077.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C077.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C078.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C078.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C079.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C079.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C081.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C081.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C082.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C082.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C083.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C083.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C084.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C084.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C085.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C085.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C086.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C086.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C087.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C087.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C089.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C089.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C090.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C090.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C091.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C091.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C092.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C092.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C093.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C093.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C094.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C094.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C095.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C095.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C096.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C096.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C098.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C098.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C099.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C099.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C100.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C100.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C101.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C101.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C104.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C104.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C105.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C105.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C106.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C106.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C107.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C107.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C108.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C108.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C109.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C109.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C110.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C110.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C111.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C111.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C112.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C112.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C113.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C113.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C114.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C114.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C115.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C115.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C118.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C118.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C119.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C119.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C121.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C121.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C122.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C122.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C123.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C123.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C124.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C124.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C126.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C126.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C127.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C127.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C128.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C128.bam.txt
samtools mpileup /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C129.bam | awk -v X=5 '$4>=X' | wc -l > breadth_C129.bam.txt
