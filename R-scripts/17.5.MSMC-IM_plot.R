# Anik Dutta

# Script to plot the results from MSMC2-IM (Isolation-migration model). It is an incomplete analysis. Do not rely on this. Rerun if you want to publish.
# The output file from the code '17.3.MSMC2_IM.sh' looks like the following:

# left_time_boundary	im_N1	im_N2	m	M
# 0.0	40304.31166516613	9081.699592009832	7.720610339544217e-33	0.0
# 482.908	40304.31166516613	9081.699592009832	7.720610339544217e-33	0.0

# Read the output file

a<-read.table("2747_4330.haps.b1_1e-09.b2_1e-05.MSMC_IM.estimates.txt", header = T)

# This code works
library(stringr)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(tidyr)
library(RColorBrewer)
library(grid)
gen<-5 # generation time
a$time<-gen*a$left_time_boundary

#plot((a$time),a$M,type='s',log = "x",col="black",ylab='Cumulative migration',xlab='log10 years before present')

q <- quantile(a$left_time_boundary, probs = c(0.25, 0.5, 0.75))

plot_m <- ggplot(a,aes(x = left_time_boundary, y =  m))  + 
  geom_step(mapping = aes(x=left_time_boundary, y=m), linetype="solid", color="red") +
  scale_x_continuous( trans = 'log10') + scale_y_continuous(labels=function(x)x*100) +
  theme_bw() +
  theme( axis.text.x = element_text(size=12), axis.title = element_text(size = 14), 
        axis.text.y = element_text(size=12), strip.text.y = element_text(angle = 180,size=8,face="bold"), strip.placement = "outside",
        strip.background = element_blank(), legend.position="none")+ labs(y="migration rate (m)", x = "Years ago")


ggsave(file = paste0("2747_4330_09_05.pdf"), width = 8, height = 6, dpi = 300, units = "in", device='pdf', useDingbats=F)

# geom_vline(xintercept = q[2],linetype=4, color="red") +

