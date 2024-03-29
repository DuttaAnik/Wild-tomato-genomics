# Anik Dutta

#This script is for calculating genome coverage for a batch of files post read-mapping
#Change to the appropriate working directory were all coverage files are stored

library(dplyr)
library(ggplot2)
library(ggridges)

#Place all files with .cov extension in vector
cov_files <- list.files(pattern="*.cov", full.names=TRUE)

#Can include 2012 like this
#cov_files <- list.files(pattern="^[19].*\\.cov$", full.names=TRUE)

##Read in files and rename columns
#Rename columns vector
covcols <- c("contig", "start", "end", "coverage")

#Make a function to read in files
loadCovFile <- function(x) {
  #read in a coverage file, extract the sample name from the file, and add it as a column
  
  df <- read.delim(x, header=FALSE, col.names=covcols)
  df$sample_name <- sub("(.+)\\.cov", "\\1", basename(x))
  df
}

#Then use lapply to apply this function to all files and load in as list
cov <- lapply(cov_files, loadCovFile)

#Then, make this list into one long dataframe
cov_df <- do.call(rbind, cov)

#Average coverage per sample
avg_cov <- cov_df %>% group_by(sample_name) %>%
  summarize(average_cov = mean(coverage)) %>%
  arrange(desc(average_cov))

write.table(avg_cov, file = "average_coverage.txt", sep = "\t", row.names = FALSE)

#Average coverage per sample/per contig
avg_cov_contig <- cov_df %>% group_by(.dots = c("sample_name", "contig")) %>%
  summarize(average_cov = mean(coverage)) %>%
  arrange(desc(average_cov))

write.table(avg_cov_contig, file = "average_coverage_per_contig.txt", sep = "\t", row.names = FALSE)

#now make plots
#Set theme to be classic
theme_set(theme_classic() + 
            theme(
              legend.text=element_text(size=11),
              axis.text=element_text(size=9),
              axis.title = element_text(size=15),
              axis.line = element_line(size=0.5),
              strip.text.y = element_blank()
            )
)

#Ridge plots
plot_cov <- ggplot(cov_df, aes(x=coverage, y=sample_name)) +
  geom_density_ridges(scale=1.75, fill="tomato2") + 
  scale_y_discrete(expand=c(0.001, 0)) +
  #facet_wrap(~sample_year) +
  xlim(0, 100) +
  theme(strip.text.x = element_blank()) +
  labs(x = "Average Depth Per Contig", y="")

# Save the plots
png("coverage_plot.png", units="in", width=15, height=10, res=300)
plot_cov
dev.off()

pdf("coverage_plot.pdf")
plot_cov
dev.off()