#!/bin/bash

# Anik Dutta

# Running MSMC2 for Cross Coalescent analysis. If you’re interested in looking at the split between two populations or species, then you’ll want to do things slightly differently, and will need to actually make three MSMC runs. You should generate one combined input file with all of the individuals you want to include in analyses. Then, you’ll run MSMC three times, using different indices: once for each of the two populations separately, and once for the two populations combined. 

# Here, I take two individuals (4 haplotypes from each population and then run pairwise msmc2 in all possible combination of those haplotypes to obtain Cross Coalescent estimates)


# First, 1958 and 1963
msmc2_linux64bit -t 12 -I 0,1,2,3 -s -o 1958_pop MSMC_input_combined.Scaffold_*.txt 
msmc2_linux64bit -t 12 -I 6,7,8,9 -s -o 1963_pop MSMC_input_combined.Scaffold_*.txt 
msmc2_linux64bit -t 12 -I 0-6,0-7,0-8,0-9,1-6,1-7,1-8,1-9,2-6,2-7,2-8,2-9,3-6,3-7,3-8,3-9 -s -o 1958_1963_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 12,13,14,15 -s -o 2747_pop MSMC_input_combined.Scaffold_*.txt 
msmc2_linux64bit -t 12 -I 0-12,0-13,0-14,0-15,1-12,1-13,1-14,1-15,2-12,2-13,2-14,2-15,3-12,3-13,3-14,3-15 -s -o 1958_2747_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 18,19,20,21 -s -o 2931_pop MSMC_input_combined.Scaffold_*.txt 
msmc2_linux64bit -t 12 -I 0-18,0-19,0-20,0-21,1-18,1-19,1-20,1-21,2-18,2-19,2-20,2-21,3-18,3-19,3-20,3-21 -s -o 1958_2931_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 24,25,26,27 -s -o 2932_pop MSMC_input_combined.Scaffold_*.txt 
msmc2_linux64bit -t 12 -I 0-24,0-25,0-26,0-27,1-24,1-25,1-26,1-27,2-24,2-25,2-26,2-27,3-24,3-25,3-26,3-27 -s -o 1958_2932_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 30,31,32,33 -s -o 3111_pop MSMC_input_combined.Scaffold_*.txt 
msmc2_linux64bit -t 12 -I 0-30,0-31,0-32,0-33,1-30,1-31,1-32,1-33,2-30,2-31,2-32,2-33,3-30,3-31,3-32,3-33 -s -o 1958_3111_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 36,37,38,39 -s -o 4107_pop MSMC_input_combined.Scaffold_*.txt 
msmc2_linux64bit -t 12 -I 0-36,0-37,0-38,0-39,1-36,1-37,1-38,1-39,2-36,2-37,2-38,2-39,3-36,3-37,3-38,3-39 -s -o 1958_4107_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 42,43,44,45 -s -o 4117_pop MSMC_input_combined.Scaffold_*.txt 
msmc2_linux64bit -t 12 -I 0-42,0-43,0-44,0-45,1-42,1-43,1-44,1-45,2-42,2-43,2-44,2-45,3-42,3-43,3-44,3-45 -s -o 1958_4117_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 48,49,50,51 -s -o 4330_pop MSMC_input_combined.Scaffold_*.txt 
msmc2_linux64bit -t 12 -I 0-48,0-49,0-50,0-51,1-48,1-49,1-50,1-51,2-48,2-49,2-50,2-51,3-48,3-49,3-50,3-51 -s -o 1958_4330_pop MSMC_input_combined.Scaffold_*.txt 


# 1963 and all
msmc2_linux64bit -t 12 -I 6-0,6-1,6-2,6-3,7-0,7-1,7-2,7-3,8-0,8-1,8-2,8-3,9-0,9-1,9-2,9-3 -s -o 1963_1958_pop MSMC_input_combined.Scaffold_*.txt

msmc2_linux64bit -t 12 -I 6-12,6-13,6-14,6-15,7-12,7-13,7-14,7-15,8-12,8-13,8-14,8-15,9-12,9-13,9-14,9-15 -s -o 1963_2747_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 6-18,6-19,6-20,6-21,7-18,7-19,7-20,7-21,8-18,8-19,8-20,8-21,9-18,9-19,9-20,9-21 -s -o 1963_2931_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 6-24,6-25,6-26,6-27,7-24,7-25,7-26,7-27,8-24,8-25,8-26,8-27,9-24,9-25,9-26,9-27 -s -o 1963_2932_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 6-30,6-31,6-32,6-33,7-30,7-31,7-32,7-33,8-30,8-31,8-32,8-33,9-30,9-31,9-32,9-33 -s -o 1963_3111_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 6-36,6-37,6-38,6-39,7-36,7-37,7-38,7-39,8-36,8-37,8-38,8-39,9-36,9-37,9-38,9-39 -s -o 1963_4107_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 6-42,6-43,6-44,6-45,7-42,7-43,7-44,7-45,8-42,8-43,8-44,8-45,9-42,9-43,9-44,9-45 -s -o 1963_4117_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 6-48,6-49,6-50,6-51,7-48,7-49,7-50,7-51,8-48,8-49,8-50,8-51,9-48,9-49,9-50,9-51 -s -o 1963_4330_pop MSMC_input_combined.Scaffold_*.txt 


# 2747 and all
msmc2_linux64bit -t 12 -I 12-0,12-1,12-2,12-3,13-0,13-1,13-2,13-3,14-0,14-1,14-2,14-3,15-0,15-1,15-2,15-3 -s -o 2747_1958_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 12-6,12-7,12-8,12-9,13-6,13-7,13-8,13-9,14-6,14-7,14-8,14-9,15-6,15-7,15-8,15-9 -s -o 2747_1963_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 12-18,12-19,12-20,12-21,13-18,13-19,13-20,13-21,14-18,14-19,14-20,14-21,15-18,15-19,15-20,15-21 -s -o 2747_2931_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 12-24,12-25,12-26,12-27,13-24,13-25,13-26,13-27,14-24,14-25,14-26,14-27,15-24,15-25,15-26,15-27 -s -o 2747_2932_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 12-30,12-31,12-32,12-33,13-30,13-31,13-32,13-33,14-30,14-31,14-32,14-33,15-30,15-31,15-32,15-33 -s -o 2747_3111_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 12-36,12-37,12-38,12-39,13-36,13-37,13-38,13-39,14-36,14-37,14-38,14-39,15-36,15-37,15-38,15-39 -s -o 2747_4107_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 12-42,12-43,12-44,12-45,13-42,13-43,13-44,13-45,14-42,14-43,14-44,14-45,15-42,15-43,15-44,15-45 -s -o 2747_4117_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 12-48,12-49,12-50,12-51,13-48,13-49,13-50,13-51,14-48,14-49,14-50,14-51,15-48,15-49,15-50,15-51 -s -o 2747_4330_pop MSMC_input_combined.Scaffold_*.txt 


# 2931 and all
msmc2_linux64bit -t 12 -I 18-0,18-1,18-2,18-3,19-0,19-1,19-2,19-3,20-0,20-1,20-2,20-3,21-0,21-1,21-2,21-3 -s -o 2931_1958_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 18-6,18-7,18-8,18-9,19-6,19-7,19-8,19-9,20-6,20-7,20-8,20-9,21-6,21-7,21-8,21-9 -s -o 2931_1963_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 18-12,18-13,18-14,18-15,19-12,19-13,19-14,19-15,20-12,20-13,20-14,20-15,21-12,21-13,21-14,21-15 -s -o 2931_2747_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 18-24,18-25,18-26,18-27,19-24,19-25,19-26,19-27,20-24,20-25,20-26,20-27,21-24,21-25,21-26,21-27 -s -o 2931_2932_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 18-30,18-31,18-32,18-33,19-30,19-31,19-32,19-33,20-30,20-31,20-32,20-33,21-30,21-31,21-32,21-33 -s -o 2931_3111_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 18-36,18-37,18-38,18-39,19-36,19-37,19-38,19-39,20-36,20-37,20-38,20-39,21-36,21-37,21-38,21-39 -s -o 2931_4107_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 18-42,18-43,18-44,18-45,19-42,19-43,19-44,19-45,20-42,20-43,20-44,20-45,21-42,21-43,21-44,21-45 -s -o 2931_4117_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 18-48,18-49,18-50,18-51,19-48,19-49,19-50,19-51,20-48,20-49,20-50,20-51,21-48,21-49,21-50,21-51 -s -o 2931_4330_pop MSMC_input_combined.Scaffold_*.txt 


# 2932 and all
msmc2_linux64bit -t 12 -I 24-0,24-1,24-2,24-3,25-0,25-1,25-2,25-3,26-0,26-1,26-2,26-3,27-0,27-1,27-2,27-3 -s -o 2932_1958_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 24-6,24-7,24-8,24-9,25-6,25-7,25-8,25-9,26-6,26-7,26-8,26-9,27-6,27-7,27-8,27-9 -s -o 2932_1963_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 24-12,24-13,24-14,24-15,25-12,25-13,25-14,25-15,26-12,26-13,26-14,26-15,27-12,27-13,27-14,27-15 -s -o 2932_2747_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 24-18,24-19,24-20,24-21,25-18,25-19,25-20,25-21,26-18,26-19,26-20,26-21,27-18,27-19,27-20,27-21 -s -o 2932_2931_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 24-30,24-31,24-32,24-33,25-30,25-31,25-32,25-33,26-30,26-31,26-32,26-33,27-30,27-31,27-32,27-33 -s -o 2932_3111_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 24-36,24-37,24-38,24-39,25-36,25-37,25-38,25-39,26-36,26-37,26-38,26-39,27-36,27-37,27-38,27-39 -s -o 2932_4107_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 24-42,24-43,24-44,24-45,25-42,25-43,25-44,25-45,26-42,26-43,26-44,26-45,27-42,27-43,27-44,27-45 -s -o 2932_4117_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 24-48,24-49,24-50,24-51,25-48,25-49,25-50,25-51,26-48,26-49,26-50,26-51,27-48,27-49,27-50,27-51 -s -o 2932_4330_pop MSMC_input_combined.Scaffold_*.txt 


# 3111 and all
msmc2_linux64bit -t 12 -I 30-0,30-1,30-2,30-3,31-0,31-1,31-2,31-3,32-0,32-1,32-2,32-3,33-0,33-1,33-2,33-3 -s -o 3111_1958_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 30-6,30-7,30-8,30-9,31-6,31-7,31-8,31-9,32-6,32-7,32-8,32-9,33-6,33-7,33-8,33-9 -s -o 3111_1963_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 30-12,30-13,30-14,30-15,31-12,31-13,31-14,31-15,32-12,32-13,32-14,32-15,33-12,33-13,33-14,33-15 -s -o 3111_2747_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 30-18,30-19,30-20,30-21,31-18,31-19,31-20,31-21,32-18,32-19,32-20,32-21,33-18,33-19,33-20,33-21 -s -o 3111_2931_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 30-24,30-25,30-26,30-27,31-24,31-25,31-26,31-27,32-24,32-25,32-26,32-27,33-24,33-25,33-26,33-27 -s -o 3111_2932_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 30-36,30-37,30-38,30-39,31-36,31-37,31-38,31-39,32-36,32-37,32-38,32-39,33-36,33-37,33-38,33-39 -s -o 3111_4107_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 30-42,30-43,30-44,30-45,31-42,31-43,31-44,31-45,32-42,32-43,32-44,32-45,33-42,33-43,33-44,33-45 -s -o 3111_4117_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 30-48,30-49,30-50,30-51,31-48,31-49,31-50,31-51,32-48,32-49,32-50,32-51,33-48,33-49,33-50,33-51 -s -o 3111_4330_pop MSMC_input_combined.Scaffold_*.txt 

# 4107 and all
msmc2_linux64bit -t 12 -I 36-0,36-1,36-2,36-3,37-0,37-1,37-2,37-3,38-0,38-1,38-2,38-3,39-0,39-1,39-2,39-3 -s -o 4107_1958_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 36-6,36-7,36-8,36-9,37-6,37-7,37-8,37-9,38-6,38-7,38-8,38-9,39-6,39-7,39-8,39-9 -s -o 4107_1963_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 36-12,36-13,36-14,36-15,37-12,37-13,37-14,37-15,38-12,38-13,38-14,38-15,39-12,39-13,39-14,39-15 -s -o 4107_2747_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 36-18,36-19,36-20,36-21,37-18,37-19,37-20,37-21,38-18,38-19,38-20,38-21,39-18,39-19,39-20,39-21 -s -o 4107_2931_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 36-24,36-25,36-26,36-27,37-24,37-25,37-26,37-27,38-24,38-25,38-26,38-27,39-24,39-25,39-26,39-27 -s -o 4107_2932_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 36-30,36-31,36-32,36-33,37-30,37-31,37-32,37-33,38-30,38-31,38-32,38-33,39-30,39-31,39-32,39-33 -s -o 4107_3111_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 36-42,36-43,36-44,36-45,37-42,37-43,37-44,37-45,38-42,38-43,38-44,38-45,39-42,39-43,39-44,39-45 -s -o 4107_4117_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 36-48,36-49,36-50,36-51,37-48,37-49,37-50,37-51,38-48,38-49,38-50,38-51,39-48,39-49,39-50,39-51 -s -o 4107_4330_pop MSMC_input_combined.Scaffold_*.txt 

# 4117 and all
msmc2_linux64bit -t 12 -I 42-0,42-1,42-2,42-3,43-0,43-1,43-2,43-3,44-0,44-1,44-2,44-3,45-0,45-1,45-2,45-3 -s -o 4117_1958_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 42-6,42-7,42-8,42-9,43-6,43-7,43-8,43-9,44-6,44-7,44-8,44-9,45-6,45-7,45-8,45-9 -s -o 4117_1963_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 42-12,42-13,42-14,42-15,43-12,43-13,43-14,43-15,44-12,44-13,44-14,44-15,45-12,45-13,45-14,45-15 -s -o 4117_2747_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 42-18,42-19,42-20,42-21,43-18,43-19,43-20,43-21,44-18,44-19,44-20,44-21,45-18,45-19,45-20,45-21 -s -o 4117_2931_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 42-24,42-25,42-26,42-27,43-24,43-25,43-26,43-27,44-24,44-25,44-26,44-27,45-24,45-25,45-26,45-27 -s -o 4117_2932_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 42-30,42-31,42-32,42-33,43-30,43-31,43-32,43-33,44-30,44-31,44-32,44-33,45-30,45-31,45-32,45-33 -s -o 4117_3111_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 42-36,42-37,42-38,42-39,43-36,43-37,43-38,43-39,44-36,44-37,44-38,44-39,45-36,45-37,45-38,45-39 -s -o 4117_4107_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 42-48,42-49,42-50,42-51,43-48,43-49,43-50,43-51,44-48,44-49,44-50,44-51,45-48,45-49,45-50,45-51 -s -o 4117_4330_pop MSMC_input_combined.Scaffold_*.txt 


# 4330 and all
msmc2_linux64bit -t 12 -I 48-0,48-1,48-2,48-3,49-0,49-1,49-2,49-3,50-0,50-1,50-2,50-3,51-0,51-1,51-2,51-3 -s -o 4330_1958_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 48-6,48-7,48-8,48-9,49-6,49-7,49-8,49-9,50-6,50-7,50-8,50-9,51-6,51-7,51-8,51-9 -s -o 4330_1963_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 48-12,48-13,48-14,48-15,49-12,49-13,49-14,49-15,50-12,50-13,50-14,50-15,51-12,51-13,51-14,51-15 -s -o 4330_2747_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 48-18,48-19,48-20,48-21,49-18,49-19,49-20,49-21,50-18,50-19,50-20,50-21,51-18,51-19,51-20,51-21 -s -o 4330_2931_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 48-24,48-25,48-26,48-27,49-24,49-25,49-26,49-27,50-24,50-25,50-26,50-27,51-24,51-25,51-26,51-27 -s -o 4330_2932_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 48-30,48-31,48-32,48-33,49-30,49-31,49-32,49-33,50-30,50-31,50-32,50-33,51-30,51-31,51-32,51-33 -s -o 4330_3111_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 48-36,48-37,48-38,48-39,49-36,49-37,49-38,49-39,50-36,50-37,50-38,50-39,51-36,51-37,51-38,51-39 -s -o 4330_4107_pop MSMC_input_combined.Scaffold_*.txt 

msmc2_linux64bit -t 12 -I 48-42,48-43,48-44,48-45,49-42,49-43,49-44,49-45,50-42,50-43,50-44,50-45,51-42,51-43,51-44,51-45 -s -o 4330_4117_pop MSMC_input_combined.Scaffold_*.txt 


# You then need to combine all three of these runs into a combined msmc output file using combineCrossCoal.py

python /gxfs_home/cau/suaph275/msmc-tools/combineCrossCoal.py 4330_1958_pop.final.txt 4330_pop.final.txt 1958_pop.final.txt > 4330_1958_combined.txt
python /gxfs_home/cau/suaph275/msmc-tools/combineCrossCoal.py 4330_1963_pop.final.txt 4330_pop.final.txt 1963_pop.final.txt > 4330_1963_combined.txt
python /gxfs_home/cau/suaph275/msmc-tools/combineCrossCoal.py 4330_2747_pop.final.txt 4330_pop.final.txt 2747_pop.final.txt > 4330_2747_combined.txt
python /gxfs_home/cau/suaph275/msmc-tools/combineCrossCoal.py 4330_2931_pop.final.txt 4330_pop.final.txt 2931_pop.final.txt > 4330_2931_combined.txt
python /gxfs_home/cau/suaph275/msmc-tools/combineCrossCoal.py 4330_2932_pop.final.txt 4330_pop.final.txt 2932_pop.final.txt > 4330_2932_combined.txt
python /gxfs_home/cau/suaph275/msmc-tools/combineCrossCoal.py 4330_3111_pop.final.txt 4330_pop.final.txt 3111_pop.final.txt > 4330_3111_combined.txt
python /gxfs_home/cau/suaph275/msmc-tools/combineCrossCoal.py 4330_4107_pop.final.txt 4330_pop.final.txt 4107_pop.final.txt > 4330_4107_combined.txt
python /gxfs_home/cau/suaph275/msmc-tools/combineCrossCoal.py 4330_4117_pop.final.txt 4330_pop.final.txt 4117_pop.final.txt > 4330_4117_combined.txt

# This output can then be plotted in R, using the same code in "MSMC2_plot.R"


















