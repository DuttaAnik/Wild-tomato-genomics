# Anik Dutta

# Plotting PCA from plink output. I need two outputs from the plink pca and then followed this link for the analysis
# https://speciationgenomics.github.io/pca/

# Running Principal Component Analysis (PCA) in terminal/HPC using the software PLINK on a VCF file, outputting the results in BED format and saving the PCA results.
'plink --vcf Pruned_0.2.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# \
  --extract Pruned_0.2.prune.in \
--make-bed --pca --out PCA_0.2'

# Load the required libraries for data manipulation and plotting.
library(tidyverse)
library(RColorBrewer)
library(ggplot2)

# Read the output eigenvec file generated by PLINK.
pca <- read_table("PCA_0.2.eigenvec", col_names = FALSE)

# Read the eigenvalues.
eigenval<- scan("PCA_0.2.eigenval")

# Remove the first column from the PCA dataframe.
pca<- pca[,-1]

# Renaming the columns of the pca dataframe to 'ind' for the first column and 'PC1', 'PC2' etc. for the rest of the columns.
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# Reading the population codes from the text file 'Population_ID.txt'.
pop_code <- read.table("Population_ID.txt", header=T)

# Merging the PCA data with the population codes based on individual id ('ind').
tab=merge(pca,pop_code, by="ind")

# Converting the 'Population' column to character type.
tab$Population=as.character(tab$Population)

# Define a color palette for the populations.
my_pal <- rcartocolor::carto_pal(n = 9, name = "Vivid")

# Compute the proportion of variance explained by each principal component.
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

# Generate a scree plot to visualize the proportion of variance explained by each principal component.
pdf("PCA_Screeplot_0.2.pdf", height = 4, width = 6, useDingbats = F)
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + theme(axis.title = element_text(size = 16),axis.text.x = element_text(size = 12, hjust = 1),axis.text.y = element_text(size = 12),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),panel.border = element_rect(color = "black",
                                                                         fill = NA,
                                                                         size = 1))+xlab("Principal component") + ylab("Percent variance explained")
dev.off()

# Compute the cumulative proportion of variance explained.
cumsum(pve$pve)

# Generate a plot for the first two principal components.
pdf("PCA_0.2_new.pdf", width = 4, height = 4)
b <- ggplot(tab, aes(PC1, PC2, col = Population)) + geom_point(size = 3)
b <- b + scale_colour_manual(values = c("#E58606","#5D69B1","#52BCA3","#99C945","#CC61B0","#24796C","#DAA51B","#2F8AC4","#A5AA99"))
b <- b + coord_equal() + theme_light()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))#+ggtitle("99 Chilense plants, 29921339 SNPs" )
dev.off()

# Generate a plot for the first and third principal components.
pdf("PCA1_3_0.2_new.pdf", width = 4, height = 4)
b <- ggplot(tab, aes(PC1, PC3, col = Population)) + geom_point(size = 3)
b <- b + scale_colour_manual(values = c("#7F3C8D", "#11A579","#3969AC", "#F2B701", "#E73F74", "#80BA5A", "#E68310", "#008695","#A5AA99"))
b <- b + coord_equal() + theme_light()
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)"))#+ggtitle("99 Chilense plants, 29921339 SNPs" )
dev.off()

# The following code is another version to make a multipanel plot of different PCs.
# Here, two different population files are read, labeled accordingly and merged with the PCA data to plot.

# Reload the required libraries for data manipulation and plotting.
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(cowplot)

# Read the output eigenvec file generated by PLINK.
pca <- read_table("../2_output/PCA.eigenvec", col_names = FALSE)

# Read the eigenvalues.
eigenval<- scan("../2_output/PCA.eigenval")

# Remove the first column from the PCA dataframe.
pca<- pca[,-1]

# Renaming the columns of the pca dataframe to 'ind' for the first column and 'PC1', 'PC2' etc. for the rest of the columns.
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# Read the population codes from two text files and label them accordingly.
pop_code1 <- read.table("US_pop_2015+_asexual.txt", header=F)
pop_code1$V2 <- "US_pop_2015_asexual"
pop_code2 <- read.table("US_pop_2015+_buckthorn.txt", header=F)
pop_code2$V2 <- "US_pop_2015_buckthorn"

# Combine the two population code dataframes.
pop_code <- rbind(pop_code1,pop_code2)

# Renaming the first column of the population code dataframe to 'ind'.
names(pop_code)[1] <- "ind"

# Merge the PCA data with the population codes.
tab=merge(pca,pop_code, by="ind")

# Converting the 'V2' column to character type and renaming it to 'Groups'.
tab$V2=as.character(tab$V2)
names(tab)[22] <- "Groups"

# Define a color palette for the groups.
pop1_Off_Color_Color_Palette2=c("US_pop_2015_asexual" = "#1E88E5","US_pop_2015_buckthorn" = "#004D40")

# Compute the proportion of variance explained by each principal component.
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

# Generate a bar plot to visualize the proportion of variance explained by each principal component.
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity") + xlab("Principal component") + ylab("Percent variance explained") + theme(axis.title = element_text(size = 16),axis.text.x = element_text(size = 12, hjust = 1), axis.text.y = element_text(size = 12),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),panel.border = element_rect(color = "black",fill = NA, size = 1))
a

# Generate a plot for the first two principal components.
b <- ggplot(tab, aes(PC1, PC2, col = Groups)) + geom_point(size = 3) + scale_color_manual(values=pop1_Off_Color_Color_Palette2) + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + theme_light() + theme(legend.position="none")
b

# Generate a plot for the first and third principal components.
c <- ggplot(tab, aes(PC1, PC3, col = Groups)) + geom_point(size = 3) + scale_color_manual(values=pop1_Off_Color_Color_Palette2) + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)")) + theme_light()
c

# Combine all plots in a single multi-panel plot.
p <- plot_grid(a, b, c, nrow = 1, labels = "AUTO")
p

# Save the multi-panel plot as a PDF file.
ggsave(file = paste0("pca203_PCA.pdf"), width = 16, height = 5, units = "in", device='pdf', useDingbats=F)
