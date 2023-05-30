# Anik Dutta
# Making genetic relatedness matrix using LD pruned (r2=0.2) SNP dataset. 

library(gdsfmt)
library(SNPRelate)
library(dplyr)
library(wesanderson)
library(viridisLite)
colors<-rainbow(9)
my_pal
annoCol <- list(group=c(Pop_1958=colors[1], Pop_1963=colors[2], Pop_2747=colors[3], Pop_2931=colors[4],Pop_2932=colors[5], Pop_3111=colors[6], Pop_4107=colors[7], Pop_4117=colors[8],Pop_4330=colors[9]))

# colors 
color=wes_palette("Darjeeling2", 250, type = "continuous")

# load vcffile 
bed.fn <- ("Pruned.0.2.bed")
fam.fn <- ("Pruned.0.2.fam")
bim.fn <- ("Pruned.0.2.bim")

# convert
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn,cvt.chr="char", "Pruned_0.2.gds")
SNPRelate::snpgdsVCF2GDS("imputed.vcf"
                         , out.fn = "~test.gds"
                         , method = "copy.num.of.ref")

# load ids for genetic pools

# GROUP 1
Pop_1958 <- read.table("1958_pop.txt")
colnames(Pop_1958) <- "id"
row.names(Pop_1958) <- Pop_1958$id
Pop_1958 <- Pop_1958 %>% mutate(group = "Pop_1958") %>% select(group)

# GROUP 2
Pop_1963 <- read.table("1963_pop.txt")
colnames(Pop_1963) <- "id"
row.names(Pop_1963) <- Pop_1963$id
Pop_1963 <- Pop_1963 %>% mutate(group = "Pop_1963") %>% select(group)

# GROUP 3
Pop_2747 <- read.table("2747_pop.txt")
colnames(Pop_2747) <- "id"
row.names(Pop_2747) <- Pop_2747$id
Pop_2747<- Pop_2747 %>% mutate(group = "Pop_2747") %>% select(group)

# GROUP 4
Pop_2931 <- read.table("2931_pop.txt")
colnames(Pop_2931) <- "id"
row.names(Pop_2931) <- Pop_2931$id
Pop_2931 <- Pop_2931 %>% mutate(group = "Pop_2931") %>% select(group)

# GROUP 5
Pop_2932 <- read.table("2932_pop.txt")
colnames(Pop_2932) <- "id"
row.names(Pop_2932) <- Pop_2932$id
Pop_2932 <- Pop_2932 %>% mutate(group = "Pop_2932") %>% select(group)

# GROUP 6
Pop_3111 <- read.table("3111_pop.txt")
colnames(Pop_3111) <- "id"
row.names(Pop_3111) <- Pop_3111$id
Pop_3111 <- Pop_3111 %>% mutate(group = "Pop_3111") %>% select(group)

# GROUP 7
Pop_4107 <- read.table("4107_pop.txt")
colnames(Pop_4107) <- "id"
row.names(Pop_4107) <- Pop_4107$id
Pop_4107<- Pop_4107 %>% mutate(group = "Pop_4107") %>% select(group)

# GROUP 8
Pop_4117 <- read.table("4117_pop.txt")
colnames(Pop_4117) <- "id"
row.names(Pop_4117) <- Pop_4117$id
Pop_4117 <- Pop_4117 %>% mutate(group = "Pop_4117") %>% select(group)

# GROUP 9
Pop_4330 <- read.table("4330_pop.txt")
colnames(Pop_4330) <- "id"
row.names(Pop_4330) <- Pop_4330$id
Pop_4330 <- Pop_4330 %>% mutate(group = "Pop_4330") %>% select(group)

# BIND DATAFRAMES
anot <- rbind(Pop_1958,Pop_1963,Pop_2747,Pop_2931,Pop_2932,Pop_3111,Pop_4107,Pop_4117,Pop_4330)

# LOAD VCF IN GDS FORMAT
file <- SNPRelate::snpgdsOpen("Pruned.gds")

# CALCULATE KINDSHIP MATRIX
IBD <- SNPRelate::snpgdsIBS(file, autosome.only = F, remove.monosnp = T,sample.id = as.factor(row.names(anot)))
a <- read.table("IBS.mdist")
df <- as.data.frame(IBD$ibs)
colnames(df) <- rownames(anot)
row.names(df) <- rownames(anot)

# HEATMAP
library(pheatmap)
pdf("heatmap_0.5_IBS.pdf",width = 14, height = 14, useDingbats = F)

#color=wes_palette("Zissou1", 250, type = "continuous")
out<-pheatmap(df,annotation_row = anot, annotation_col =  anot, color = color, annotation_colors = annoCol, cellwidth = 4.5, cellheight =  4.5, fontsize_row=3.5,fontsize_col = 3.5,cex=1)

#hist(IBD$ibs, breaks = 30)

#png("heatmap.png",  width = 10, height = 10, res = 300, units = "cm")
dev.off()

