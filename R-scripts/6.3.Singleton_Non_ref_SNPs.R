# Anik Dutta

# Plotting Singleton and number of nonref variants (https://z3tt.github.io/beyond-bar-and-box-plots/)

install.packages("devtools") 
devtools::install_github("jbgb13/peRReo") 
library(ggplot2)
library(peRReo)

library(tidyverse)     ## data wrangling + ggplot2
library(colorspace)    ## adjust colors
library(rcartocolor)   ## Carto palettes
library(ggforce)       ## sina plots
library(ggdist)        ## halfeye plots
library(ggridges)      ## ridgeline plots
library(ggbeeswarm)    ## beeswarm plots
library(gghalves)      ## off-set jitter
library(systemfonts) 

# This data file is made by running to following code using bcftools. Then I extracted, the number of singletons and non-reference SNP manually to make a .csv file.

'for i in 1958_pop 1963_pop 2747_pop 2931_pop 2932_pop 3111_pop 4107_pop 4117_pop 4330_pop
do
bcftools stats -S "$i".txt /gxfs_work1/cau/suaph275/Clean_VCF/ABBA_pop_stats/Chilense_99.genotyped.SNP.Clean.vcf.recode.vcf > "$i".stats
done'

# Read the data

a<-read.csv("Singleton.csv", header = T)
names(a)

my_pal <- rcartocolor::carto_pal(n = 9, name = "Bold")

pdf("Singleton.pdf",height = 6, width = 8,useDingbats = F)
g <- ggplot(a, aes(x = reorder(Population, Singleton, FUN = mean), y = Singleton, color = Population, fill = Population)) +
  
  scale_color_manual(values = my_pal, guide = "none") +ylim(0,650000)+
  scale_fill_manual(values = my_pal, guide = "none")
#g + geom_boxplot(alpha = .5, size = 1.5, outlier.size = 5)
g + stat_summary(
    geom = "point", 
    fun = mean,
    shape =11, size = 2, color = "red", stroke = 1.5
  ) + 
  ggforce::geom_sina(
    maxwidth = .5, scale = "count", 
    size = 3, alpha = .5, seed = 0
  ) + 
  ggforce::geom_sina(
    maxwidth = .5, scale = "count", 
    size = 3, alpha = .5, seed = 0,
    shape = 1, color = "black", stroke = .8
  ) + theme(axis.title = element_text(size = 16),axis.text.x = element_text(size = 12, angle = 45, hjust = 1),axis.text.y = element_text(size = 12),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(),panel.border = element_rect(color = "black",
                                                        fill = NA,
                                                        size = 1))+xlab("Population") + ylab("Number of Singleton")
dev.off()

# Read the data of non-reference SNP
b<-read.csv("Non_ref_SNP.csv", header = T)
names(b)

pdf("Non_ref_SNP.pdf",height = 6, width = 8,useDingbats = F)
g <- ggplot(b, aes(x = reorder(Population, Non_ref_SNP, FUN = mean), y = Non_ref_SNP, color = Population, fill = Population)) +
  
  scale_color_manual(values = my_pal, guide = "none") +ylim(3300000,8000000)+
  scale_fill_manual(values = my_pal, guide = "none")
#g + geom_boxplot(alpha = .5, size = 1.5, outlier.size = 5)
g + stat_summary(
  geom = "point", 
  fun = mean,
  shape =11, size = 2, color = "red", stroke = 2.5
) + 
  ggforce::geom_sina(
    maxwidth = .5, scale = "count", 
    size = 3, alpha = .5, seed = 0
  ) + 
  ggforce::geom_sina(
    maxwidth = .5, scale = "count", 
    size = 3, alpha = .5, seed = 0,
    shape = 1, color = "black", stroke = .8
  ) + theme(axis.title = element_text(size = 16),axis.text.x = element_text(size = 12, angle = 45, hjust = 1),axis.text.y = element_text(size = 12),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(),panel.border = element_rect(color = "black",
                                                                           fill = NA,
                                                                           size = 1))+xlab("Population") + ylab("Number of variants")
dev.off()


