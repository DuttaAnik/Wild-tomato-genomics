# Anik Dutta

# Script for making manhattan plot from SweeD selection scan outputs. The same can be used for plotting RAiSD result and get the significant regions under selection as a .csv file. Just need to change the headers for RAiSD.

library(ggplot2)
library(dplyr)

# load data
mydata <- read.table("All_combined_2931.txt", header = T, stringsAsFactors = FALSE, na.strings=c("", "NA"),fill = T)
names(mydata)
names(mydata)[1] <- "Contig"
names(mydata)[2] <- "Position"
names(mydata)[3] <- "Likelihood"
names(mydata)[4] <- "Alpha"
names(mydata)[5] <- "StartPos"
names(mydata)[6] <- "EndPos"

#mydata_noNA <- mydata[rowSums(is.na(mydata)) == 0,] # remove rows with NA

# calculate quantile for setting significant threshold

quantile(mydata$Likelihood, c(.99,.995,.999)) # 99%, 99.5%, 99.9%
class(mydata$Likelihood)
print(unique(mydata$Likelihood))

# subset for plot test
# mydata_chr23 <- subset(mydata, Chr == "1" | Chr == "2" | Chr == "3" | Chr == "4" | Chr == "5" | Chr == "6")

qq_999th <- unname(quantile(mydata$Likelihood, .999))
qq_999th
highlight_markers <- subset(mydata, Likelihood >= qq_999th) # regions above threshold 

# Optional (since I'll have a lot of contigs, subset and display only contigs with signigicant peak)
selected_contigs = unique(highlight_markers$Contig)
selected_contigs[2]
mydata_noNA_clean_subset <- filter(mydata, Contig %in% selected_contigs)
unique(mydata_noNA_clean_subset$Contig) == selected_contigs # all should be true

# Plot manhattan

ggplot(mydata, aes(x=Position/1000000, y=Likelihood)) + geom_point(color='darkgray') + 
  geom_point(data=highlight_markers, aes(x=Position/1000000,y=Likelihood), color='red',size=1) + 
  labs(fill = "US BT") + labs(x ="Position (Mbp)", y = "Likelihood") +
  geom_hline(yintercept=qq_999th, linetype="dashed", color = "black", size=0.8) +
    
  theme(text = element_text(size=10),axis.text.x = element_text(angle =90, hjust=1),panel.background = element_rect(fill = NA),
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.line = element_line(size = 1, colour = "black")) +
   facet_grid(. ~Contig, scale='free_x')
# facet_grid(~factor(Contig, levels = c("Scaff_1", "Scaff_2", "Scaff_3", "Scaff_4", "Scaff_5", "Scaff_6", "Scaff_7","Scaff_8", "Scaff_9", "Scaff_10", "Scaff_11", "Scaff_12")), scale='free_x')+

# Save the plot

ggsave(file = paste0("2931.pdf"), width = 14, height = 6, dpi = 300, units = "in", device='pdf', useDingbats=F)

# Save the significant regions in a .csv file to work with bedtools

write.csv(highlight_markers, file = "Significant_SweeD_99.9%.csv")

