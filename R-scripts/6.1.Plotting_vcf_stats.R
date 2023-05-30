# Anik Dutta

# Script for plotting the preliminary statistics resulted from the VCF file

# Mean depth of SNPs
var_depth <- read_delim("mean_depth_site.ldepth.mean", delim = "\t",
                        col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
pdf("mean_depth_1.pdf", height = 6, width = 8, useDingbats = F)
a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
a + theme_light() + xlim(0, 100)
dev.off()
summary(var_depth$mean_depth)

# Heterzygosity
pdf("hetero.pdf", height = 6, width = 8, useDingbats = F)
ind_het <- read_delim("het.het", delim = "\t",
                      col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)
a <- ggplot(ind_het, aes(f)) + geom_histogram(fill = "tomato3", colour = "black", alpha = 0.3)
a + theme_light()+xlab("Heterozygosity")
dev.off()

# Missingness per individual
pdf("miss_ind.pdf", height = 6, width = 8, useDingbats = F)
ind_miss  <- read_delim("miss_indv.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "skyblue", colour = "black", alpha = 0.3)
a + theme_light()+xlab("proportion of missing data per individual")
dev.off()

# Depth per individual
pdf("depth_ind.pdf", height = 6, width = 8, useDingbats = F)
ind_depth <- read_delim("depth_ind.idepth", delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1)  
a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "orange", colour = "black", alpha = 0.3)
a + theme_light()+xlab("Depth per individual")
dev.off()


var_freq <- read_delim("allele_freq.frq", delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)

# Minor allele frequency
var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))
pdf("MAF.pdf", height = 6, width = 8, useDingbats = F)
a <- ggplot(var_freq, aes(maf)) + geom_density(fill = "turquoise2", colour = "black", alpha = 0.3)
a + theme_light()+xlab("Minor allele frequency")
dev.off()
summary(var_freq$maf)

# Missingness in SNP
var_miss <- read_delim("miss_site.lmiss", delim = "\t",
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)

pdf("Missingness_SNP.pdf", height = 6, width = 8, useDingbats = F)
a <- ggplot(var_miss, aes(fmiss)) + geom_density(fill = "darkmagenta", colour = "black", alpha = 0.3)
a + theme_light()+xlab("SNP missingness")
dev.off()
summary(var_miss$fmiss)

# Qualit
var_qual <- read_delim("site_qual.lqual", delim = "\t",
                       col_names = c("chr", "pos", "qual"), skip = 1)

pdf("Quality.pdf", height = 6, width = 8, useDingbats = F)
a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "lightsalmon1", colour = "black", alpha = 0.3)
a + theme_light()
dev.off()

pdf("Zoom_Quality.pdf", height = 6, width = 8, useDingbats = F)
a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "lightsalmon1", colour = "black", alpha = 0.3)
a + theme_light()+xlim(0, 10000)
dev.off()












