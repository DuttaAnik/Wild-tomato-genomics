#!/bin/bash

# Anik Duta

# Use after finding genes using bedtools with selection scan results.

# How to find overlapping genes from text files from multiples groups in R. This is for making the upset plot to show how many genes under selection are unique in each population and shared among populations. 

# load the data frames for each group
group1 <- read.table("group1.txt", header = TRUE)
group2 <- read.table("group2.txt", header = TRUE)
group3 <- read.table("group3.txt", header = TRUE)

# combine the data frames into a single data frame
all_groups <- rbind(group1, group2, group3)

# extract the unique gene names
group1_genes <- unique(group1$gene_name)
group2_genes <- unique(group2$gene_name)
group3_genes <- unique(group3$gene_name)

# find the overlapping genes
overlapping_genes <- intersect(group1_genes, group2_genes, group3_genes)

# store the overlapping genes
write.table(overlapping_genes, "overlapping_genes.txt", row.names = FALSE, col.names = FALSE)

## Plot it
# count the number of overlapping genes in each group
group1_overlap <- sum(group1_genes %in% overlapping_genes)
group2_overlap <- sum(group2_genes %in% overlapping_genes)
group3_overlap <- sum(group3_genes %in% overlapping_genes)

# create a data frame with the number of overlapping genes
overlap_counts <- data.frame(group = c("group1", "group2", "group3"),
                             overlap = c(group1_overlap, group2_overlap, group3_overlap))

# plot the results using a bar plot
library(ggplot2)
ggplot(overlap_counts, aes(x = group, y = overlap)) +
  geom_col() +
  xlab("Group") +
  ylab("Number of Overlapping Genes") +
  ggtitle("Number of Overlapping Genes in Each Group")

# Use upset plot

# install the UpSetR library if necessary
if (!require(UpSetR)) {
  install.packages("UpSetR")
}
library(UpSetR)

# create a list of the gene sets for each group
gene_sets <- list(group1 = group1_genes, group2 = group2_genes, group3 = group3_genes)

# create the upset plot
upset(gene_sets, sets = c("group1", "group2", "group3"), main.bar.size = 0.4)
