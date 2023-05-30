# Anik Dutta

# Script for plotting Fst and Dxy from the output of Pixy tool

a<-read.csv("Chilense.Fst.Dxy.pi.csv",header = T)
a[a < 0] <- 0
colMeans(x=a, na.rm = TRUE)

library(dplyr)
a$sites<-NULL
t(apply(a, 1, function(x) tapply(x, colnames(a), mean)))
aa<-as.data.frame (colMeans(a, na.rm = T))          

# Get lower triangle of the correlation matrix
get_lower_tri<-function(a){
  a[upper.tri(a)] <- NA
  return(a)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(a){
  a[lower.tri(a)]<- NA
  return(a)
}

a<-read.csv("Dxy_pairwise_Pixy.csv", row.names = 1)
a<-read.csv("Fst_pairwise_Pixy.csv", row.names = 1)
a<-as.matrix(a)
upper_tri <- get_lower_tri(a)
upper_tri
library(reshape2)
melted_cormat <- melt(upper_tri, na.rm = TRUE)

upper_tri<-round(upper_tri,2)

#Plot the figure
pdf("pairwise_fst.pdf",height = 7, width = 7, useDingbats = F)
corrplot(upper_tri,method="circle",col.lim = c(0.10,0.51),addCoef.col = 'tomato2',number.cex = 0.8,col = COL1('Blues'),
         tl.col="black", tl.srt=45, is.corr = FALSE,diag = FALSE,addgrid.col = NA,
         na.label = "square",na.label.col = "white")
dev.off()
