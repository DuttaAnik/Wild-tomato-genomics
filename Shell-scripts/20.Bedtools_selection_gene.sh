#!/bin/bash

# Working with bedtools to extract genes under selection from significant selective sweep regions. Can be used for SweeD and RAiSD results (Anik Dutta)

# Converting a gff file into bed format without sorting for chromosomes. gff2bed is coming from bedops package.

gff2bed < chilense_HiRise.helixer_land_plant_v0.3a_filtered.gff3 --do-not-sort > Chilense.bed

# Subsetting the converted GFF to BED format file for a set of chromosomes. First, make a txt file with only one column stating the chromosome names and then:

grep -f scaff_names.txt Chilense.bed > Chilense_12_Scaffold.bed

# Then I manually omitted the rest of scaffold text from the name in BBEdit.

# Here, I deleted all the exons mrnas from the converted annotation file to keep only the genes. $8 denotes the column with "gene"

awk '{
OFS="\t"
if ($8 == "gene")
print }' Chilense_12_Scaffold.bed > new_gene.bed

# Convert the .csv file containing significant regions under selection above 99.5% threshold from R.

awk -F, 'NR > 1 {gsub(/"/, "", $2); print $2"\t"$6"\t"$7}' Significant_SweeD_99.5%.csv > Significant.bed.txt

# bedtools sort -i Significant.bed.txt > sorted.significant.bed (do it later)

# Merged significant regions under selection

bedtools merge -i Significant.bed.txt > merged.bed

# Add flanking region in bp to each start and stop site based on LD. Here I added 1000bp on each side considering 1958_C LD i.e. 2.1 kb (https://bedtools.readthedocs.io/en/latest/content/tools/slop.html)

awk -v OFS='\t' '{print $1,$2-1050,$3+1050}' merged.bed > merged_flanked.bed # In case of Two_Sh populations, adding 15.3 kb on each side resulted in negative bp in some cases. so I put a 0 instead of negative bp.

# sort it and the annotation file
bedtools sort -i merged_flanked.bed > sorted.merged_flanked.bed

bedtools sort -i Chilense_12_Scaffold.bed > sorted.Chilense_12_Scaffold.bed

# Use intersect to find overlapping regions

#bedtools intersect -a sorted.merged.bed -b sorted.new_gene.bed -wa > results.txt # use the -u option and -wb as well for final analysis

#bedtools intersect -a sorted.new_gene.bed -b sorted.merged.bed -wa > results_1.txt (which one to choose from these two codes?)

#bedtools intersect -a sorted.merged.bed -b sorted.new_gene.bed -wb > results_2.txt (this and above give the same number of genes)

#bedtools intersect -a sorted.merged.bed -b sorted.new_gene.bed -u -wa | wc -l (count the number of genes)

## Final code to find the genes.

bedtools intersect -a sorted.new_gene.bed -b sorted.merged_flanked.bed -u -wa -wb > results.bed

# I used the following code to select only 12 scaffolds from the pacbio reference chilense genome. I have 12 scaffolds name in a text file. (https://www.biostars.org/p/319099/)

seqkit grep -n -f scaff_names.txt S_chilense_Hirise.fasta > S_chilense_Hirise_12_C.fasta # I then manually changed the Scaffolds name to match the annotation scaffold name e.g. Scaffold_1

# Transform a fasta genome file into bed format. This was done to add flanking region on both sides of the interval. But the above awk code worked.

pip install pyfaidx  
faidx --transform bed S_chilense_Hirise_12_C.fasta > S_chilense_Hirise_12_C.bed

# Extract fasta sequence from genome file using a txt file containing the interval or coordinates of the genes from the intersect output above. (https://bedtools.readthedocs.io/en/latest/content/tools/getfasta.html)

awk -v OFS='\t' '{print $1,$2,$3}' results.bed > results_1.bed # this code deletes all the other column from the result file containing chr, start, stop, gene id etc. It only keeps the chr, start and stop in a tab delimited bed file (https://unix.stackexchange.com/questions/595839/how-to-remove-column-from-file-without-altering-format)

bedtools getfasta -fi /Users/anikdutta/Documents/S_chilense_genomes/New_PB_HiRise_genome/S_chilense_Hirise_12_C.fasta -bed results_1.bed -fo results.fa.out # then submit this fasta file to Mercator4 for validation and subsequent annotation.

# If there are duplicate fasta entries, then use the following code to find duplicate entries in the fasta file

# Save the file names to an array
file_array=($(awk '/^>/ {print $1}' file.fasta))

# Find duplicate entries
duplicates=($(echo "${file_array[@]}" | tr ' ' '\n' | sort | uniq -d))

# Print duplicates
for duplicate in "${duplicates[@]}"; do
    echo "Duplicate: $duplicate"
done

# While running on Mercator, some sequences were identified as invalid. So, I created the following function to find out which sequence were invalid in the input file compared to Mercator generated input fasta file.

function extract_sequence {
  awk '/^>/ {if (seqlen) {print seq}; print ;seq="";next; } { seq = seq $0 } END {print seq}' $1
}

# Find the unique sequence in the original input fasta file

comm -23 <(extract_sequence results.fasta | sort) <(extract_sequence MercatorValidFasta.fa | sort)
##

# Only select fasta header from a fasta file

grep "^>" MercatorValidFasta.fa > Mercator_fasta_header.txt

# The fasta header in the amino acid annotation file was like >chilenese_blah_blah_voila gene=chilense_blah_blah. So I had to change it to >chilense_blah_blah

sed 's/^>.*gene=/>gene=/' chilense_HiRise.helixer_land_plant_v0.3a_filtered.aa.fa > new.aa.fa

sed 's/^>.*Scaffold_/>Scaffold_/' new.aa.fa > new_1.aa.fa # then i manually removed the traditional Scaffolds name and matched with my scaffold name e.g. Scaffold_1

# I had to edit the name Gene ID in the results.bed file from the intersect command. Because AA fasta file contains the gene ID that I need to match with my results to get the aa fasta sequences.

awk -v OFS='\t' '{sub(/ID=chilense_Hirise_/,"",$10); print}' results.bed > results_2.bed # here I deleted the ID=chilense_Hirise_ part from the gene ID in the output file.

awk -v OFS='\t' '{print $10}' results_2.bed > gene_ID.txt

seqkit grep -n -f gene_ID.txt /Users/anikdutta/Documents/S_chilense_genomes/New_PB_HiRise_genome/new_1.aa.fa > AA.fasta



## Codes for automation (See the 18.SweeD_selection_scan.sh last part). Always change the *_highland part for different genetic groups. Check the naming of the output file from SweeD

for file in SweeD_Report.*_highland.txt; do           
  new_name=$(echo "$file" | sed 's/SweeD_Report\.Scaffold_\([0-9]*\)__.*/Scaffold_\1/')         
  mv "$file" "$new_name"
done

for i in Scaffold*
do            
awk '{print FILENAME"\t"$0}' $i > $i.bk
mv $i.bk $i
done                    

for file in Scaffold*
do
  tail -n +2 "$file" > "$file.tmp" && mv "$file.tmp" "$file"
done

for file in Scaffold* 
do
  tail -n +2 "$file" > "$file.tmp" && mv "$file.tmp" "$file"
done

cat Scaffold_1 Scaffold_2 Scaffold_3 Scaffold_4 Scaffold_5 Scaffold_6 Scaffold_7 Scaffold_8 Scaffold_9 Scaffold_10 Scaffold_11 Scaffold_12 > All_combined_2931.txt
awk -F, 'NR > 1 {gsub(/"/, "", $2); print $2"\t"$6"\t"$7}' Significant_SweeD_99.5%.csv > Significant.bed.txt
bedtools sort -i Significant.bed.txt > sorted.significant.bed

#bedtools intersect -a sorted.new_gene.bed -b sorted.significant.bed -u -wa -wb > results_1.bed 

bedtools merge -i sorted.significant.bed > merged.bed
awk -v OFS='\t' '{print $1,$2-15300,$3+15300}' merged.bed > merged_flanked.bed
bedtools intersect -a sorted.new_gene.bed -b merged_flanked.bed -u -wa -wb > results.bed 



