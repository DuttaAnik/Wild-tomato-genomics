# Anik Dutta

# Plotting Admixture graph using Pophelper

# Install the 'pophelper' package from GitHub.
remotes::install_github('royfrancis/pophelper')

# Load the 'pophelper' and 'ggplot2' libraries.
library(pophelper)
library(ggplot2)

# List all the files in the current directory.
afiles <- list.files(path=".", full.names=T)

# Read Q-files (a type of data file) in the current directory.
alist <- readQ(files=afiles)

# Read Q-files in 'basic' format.
readQ(files=afiles, filetype = 'basic')

# Read a text file containing population ID labels.
labels<-read.csv("Population_ID.txt",header=F, stringsAsFactors=F, sep = "\t")

# List all the files in the current directory again.
sfiles <- list.files(path=".", full.names=T)

# Read Q-files and sort them based on individual labels in the file.
slist <- sortQ(readQ(sfiles,indlabfromfile=T))

# Plot admixture graph using specific parameters, return the plot instead of exporting it.
p1 <- plotQ(slist[c(1,2,3,4)], imgoutput="join", returnplot=T, exportplot=F, basesize=3,
            grplab=labels, grplabsize=3.5, showyaxis=T, showticks=T, linesize=0.5, pointsize=3, panelratio=c(5,1),
            splabsize=8, grplabangle=45, grplabpos=0.6, grplabheight=1, grplabjust=1,
            selgrp="V2", ordergrp=T, indlabvjust=1, clustercol = my_pal)

# Arrange and display the plot on a grid.
grid.arrange(p1$plot[[1]])

# Read the population ID labels again.
labels<-read.csv("Population_ID.txt",header=F, stringsAsFactors=F, sep = "\t")

# Create a new variable 'twolabelset' for two label sets.
twolabelset <- labels
twolabelset$V1 <- NULL  # Remove the 'V1' column from the label set.

# Create a final figure in PDF format using specific parameters for the admixture plot.
pdf("Admixture.pdf",width = 10, height = 10, useDingbats = F)
p1<-plotQ(slist[c(1,2,3,4)], imgoutput="join", barsize = 1, panelspacer = 0.3, divtype = 1, returnplot=T, exportplot=F, showyaxis=T, basesize=10,
          grplab=labels, grplabangle=45, grplabheight = 0.2, grplabpos = 0.8, grplabsize=2, linesize=0.8, pointsize=3, clustercol = my_pal)

# Arrange and display the plot on a grid.
grid.arrange(p1$plot[[1]])
dev.off()  # Close the PDF device.

# Read a text file 'CV.txt' containing data for Cross Validation error.
b<-read.table("CV.txt", header = T)

# Convert the 'K' column to factor with levels defined by 'K'.
b$K <- factor(b$K , levels = b$K )

# Create a final figure in PDF format for Cross validation error to show which K value is the best.
pdf("Admixture_CV.pdf", height = 4, width = 6, useDingbats = F)
ggplot(data = b, aes(x = K, y = CV,group = 1)) + geom_point(size = 3) + geom_line(lty = 3) + 
  xlab("K") + ylab("Cross Validation error") + theme(axis.title = element_text(size = 16),axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                     panel.background = element_blank(),panel.border = element_rect(color = "black",
                                                                                                                    fill = NA,
                                                                                                                    size = 1))+xlab("K") + ylab("Cross validation error")
dev.off()  # Close the PDF device.