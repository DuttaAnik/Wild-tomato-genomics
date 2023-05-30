#!/bin/bash

# Anik Dutta

# Selection scan using SweeD. SweeD was run on each genetic cluster identified by PCA and Admixture. There were 7 genetic clusters where the two_central.txt contained the populations 2747 and 3111, and the two_highland.txt contained the populations 4117 and 4330.

# Here, I have used 10 kb window measured by varying grid size depending on the length of the scaffold. For example, length 78070536/10000 bp = 7807 (-folded) is used as I dont know which allele is ancestral and which one is derived.

# Estimate grid value for each chromosome. Here, the grid value is set for 10 kb

head -n 100 Scaffold_1__2_contigs__length_78070536.vcf | grep -f scaff_names.txt | sed 's/##contig=<ID=//g; s/,length=/\t/g; s/>//g' | gawk '{printf "%s\t%d\n", $1,$2/10000+1}' > list-chr-gird.txt

# Define the path to the directory containing the input files

input_dir="path/to/input/files" # did not use

# Define the path to the directory to store the output files

output_dir="path/to/output/files" # did not use

# Define the grid values for each chromosome

grid_values=(7807 6173 11444 6022 8850 8359 8097 7689 7065 6970 7515 7144)

# Define the names of the population files

population_files=(1958_pop.txt 1963_pop.txt 2931_pop.txt 2932_pop.txt 4107_pop.txt two_central.txt two_highland.txt)

# Define the names of the chromosomes

chromosomes=(Scaffold_1__2_contigs__length_78070536 Scaffold_2__1_contigs__length_61731525 Scaffold_3__3_contigs__length_114446166 Scaffold_4__2_contigs__length_60228994 Scaffold_5__2_contigs__length_88504823 Scaffold_6__2_contigs__length_83591330 Scaffold_7__2_contigs__length_80977755 Scaffold_8__2_contigs__length_76899493 Scaffold_9__1_contigs__length_70652527 Scaffold_10__2_contigs__length_69708562 Scaffold_11__3_contigs__length_75158268 Scaffold_12__3_contigs__length_71445296)

# Run SweeD for each scaffold with each pop with different grid values set at 10kb. grid value= chr_length/10kb. 
# Main code

for i in "${!chromosomes[@]}" 
do 
chromosome="${chromosomes[$i]}" 
grid="${grid_values[$i]}" 
for population_file in "${population_files[@]}" 
do SweeD -input "$chromosome.vcf" -folded -grid "$grid" -noSeparator -sampleList "$population_file" -name "$chromosome-$population_file"  
done 
done


# Loop through each chromosome and corresponding grid value. Run this if you wish to use input_dir and output_dir.
#for i in "${!chromosomes[@]}"; do

  # Get the current chromosome and grid value
  #chromosome="${chromosomes[$i]}"
  #grid="${grid_values[$i]}"

  # Loop through each population file
  #for population_file in "${population_files[@]}"; do

    # Run SweeD on the current chromosome and population file
    #sweed -name "$chromosome-$population_file" -input "$input_dir/$chromosome.vcf" -pop "$input_dir/$population_file" -output "$output_dir/$chromosome-$population_file" -grid "$grid"

  #done

#done

# Script for formating the output files from SweeD to a readable version to work in R to get the significant regions under selection.

# Codes for automation. Always change the *_highland part for different genetic groups. Check the naming of the output file from SweeD

# Use 'sed' command to manipulate the filename, extracting 'Scaffold_' and the following digits.
# Rename the file with the new name.

for file in SweeD_Report.*_highland.txt; do
  new_name=$(echo "$file" | sed 's/SweeD_Report\.Scaffold_\([0-9]*\)__.*/Scaffold_\1/')
  mv "$file" "$new_name"
done

# Loop over all files starting with 'Scaffold'
# Use 'awk' to prepend each line of the file with the filename
# Rename the resulting file to the original filename

for i in Scaffold*
do
  awk '{print FILENAME"\t"$0}' $i > $i.bk
  mv $i.bk $i
done

# Loop over all files starting with 'Scaffold'
# Use 'tail' to remove the first line of each file (header) and write the result to a temporary file.
# Rename the temporary file to the original filename.

for file in Scaffold*
do
  tail -n +2 "$file" > "$file.tmp"
  mv "$file.tmp" "$file"
done

# The below block seems redundant as it's doing the same operation as the block above - removing the first line of each file.
# Use 'tail' to remove the first line of each file (header) and write the result to a temporary file.
# Rename the temporary file to the original filename
# Loop over all files starting with 'Scaffold'

for file in Scaffold* 
do
  tail -n +2 "$file" > "$file.tmp"
  mv "$file.tmp" "$file"
done

# Concatenate all the scaffolds into one result file

cat Scaffold_1 Scaffold_2 Scaffold_3 Scaffold_4 Scaffold_5 Scaffold_6 Scaffold_7 Scaffold_8 Scaffold_9 Scaffold_10 Scaffold_11 Scaffold_12 > All_combined_2931.txt

# Use this 'All_combined___.txt' file for plotting and finiding significant regions under selection in R.
