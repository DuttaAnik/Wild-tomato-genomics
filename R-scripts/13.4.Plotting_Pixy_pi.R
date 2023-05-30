# Anik Dutta
# Script to process output file from Pixy and plot nucleotide diversity (pi)

library(tidyverse)
library(tidyverse)
library(pheatmap)
library(corrplot)
library(ComplexHeatmap)
library(circlize)
library(lazyeval)

# Function to convert Pixy output files to long format for easy manipulation and analysis
pixy_to_long <- function(pixy_files){
  
  # Initializing empty list to store dataframes
  pixy_df <- list()
  
  # Loop over all Pixy files
  for(i in 1:length(pixy_files)){
    
    # Extract file type from file name
    stat_file_type <- gsub(".*_|.txt", "", pixy_files[i])
    
    # Different handling based on file type
    if(stat_file_type == "pi"){
      
      # Read the file and convert it to long format
      df <- read_delim(pixy_files[i], delim = "\t")
      df <- df %>%
        gather(-pop, -window_pos_1, -window_pos_2, -chromosome,
               key = "statistic", value = "value") %>%
        rename(pop1 = pop) %>%
        mutate(pop2 = NA)
      
      # Add the dataframe to the list
      pixy_df[[i]] <- df
      
    } else{
      
      # Read the file and convert it to long format
      df <- read_delim(pixy_files[i], delim = "\t")
      df <- df %>%
        gather(-pop1, -pop2, -window_pos_1, -window_pos_2, -chromosome,
               key = "statistic", value = "value")
      
      # Add the dataframe to the list
      pixy_df[[i]] <- df
      
    }
    
  }
  
  # Bind all dataframes in the list to a single dataframe and arrange it by several variables
  bind_rows(pixy_df) %>%
    arrange(pop1, pop2, chromosome, window_pos_1, statistic)
  
}

# Specify directory containing Pixy files
pixy_folder <- "Pixy"
# Get the full paths of the Pixy files
pixy_files <- list.files(pixy_folder, full.names = TRUE)

# Call the function to process the Pixy files and store the result in a dataframe
pixy_df <- pixy_to_long(pixy_files)
# Check column names
names(pixy_df)
# Check first few rows of dataframe
head(pixy_df)

# Read "pi" file and compute total diversity per population
allHaps_pi <- read.csv("Chilense_Pixy_output.out_pi.txt", sep = "\t")
group_by(allHaps_pi, pop) %>% dplyr::summarise(tot = sum(count_diffs, na.rm = T)/sum(count_comparisons, na.rm = T))

# Compute "new_pi" by dividing count_diffs by count_comparisons
allHaps_pi$new_pi<-(allHaps_pi$count_diffs, na.rm = T)/(allHaps_pi$count_comparisons, na.rm = T)
allHaps_pi$new_pi<- (allHaps_pi$count_diffs)/(allHaps_pi$count_comparisons)

# Read "dxy" file and compute total divergence per population pair
filt_Haps_DXY <- read.csv("Chilense_Pixy_output.out_dxy.txt", sep = "\t")
allDXY_sum <- group_by(filt_Haps_DXY, pop1, pop2) %>% summarise(tot = sum(count_diffs, na.rm = T)/sum(count_comparisons, na.rm = T))

# Read "fst" file and compute average FST per population pair
filt_Haps_Fst <- read.csv("Chilense_Pixy_output.out_fst.txt", sep = "\t")
names(filt_Haps_Fst)

# Set negative FST values to 0
filt_Haps_Fst$avg_wc_fst[filt_Haps_Fst$avg_wc_fst < 0] <- 0
allfst_sum <- group_by(filt_Haps_Fst, pop1, pop2) %>% summarise(tot = mean(avg_wc_fst, na.rm = T))

# Plot the "new_pi" values using ggplot2
library(ggplot2)
pdf(file="Pi_pixy.pdf")
ggplot(allHaps_pi,aes(x=reorder(factor(pop),new_pi,FUN=mean),y=new_pi)) +
  geom_violin(aes(fill = factor(pop)))+
  stat_summary(fun=mean, geom="point", size=2, color="black") + 
  scale_fill_manual(values=mypal2) + ylab(expression(pi)) + 
  theme(axis.title.y = element_text(size = rel(1.8)))+ theme(legend.position="none") +
  xlab("Population")  + theme(axis.title = element_text(size = 16),axis.text.x = element_text(size = 12, angle = 45, hjust = 1),axis.text.y = element_text(size = 12),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                              panel.background = element_blank(),panel.border = element_rect(color = "black",
                                                                                             fill = NA,
                                                                                             size = 1))
dev.off()

