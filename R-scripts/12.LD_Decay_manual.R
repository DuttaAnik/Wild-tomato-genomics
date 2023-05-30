# Anik Dutta
# Script for plotting linkage-disequilibrium decay manually with the output files from PopLDdecay tool

# Load required libraries
library(tidyverse)
library(viridisLite)

# Define color palette
my_pal <- rcartocolor::carto_pal(n = 9, name = "Vivid")

# Read in files as data frames
read.table("all_pop_1000.Pop_1958")->c1
read.table("all_pop_1000.Pop_1963")->c2
read.table("all_pop_1000.Pop_2747")->c3
read.table("all_pop_1000.Pop_2931")->c4
read.table("all_pop_1000.Pop_2932")->c5
read.table("all_pop_1000.Pop_3111")->c6
read.table("all_pop_1000.Pop_4107")->c7
read.table("all_pop_1000.Pop_4117")->c8
read.table("all_pop_1000.Pop_4330")->c9

# Define a function that computes expected R-squared values for Linkage Disequilibrium (LD) decay
# based on the Hill and Weir (1988) equation for a finite sample size

get_LD <- function (Distance, Rsquared, n, rho.start, group) {
  hill.weir.eq=(Rsquared~(((10+(rho*Distance))/((2+(rho*Distance))*(11+(rho*Distance))))*(1+(((3+(rho*Distance))*(12+(12*(rho*Distance))+((rho*Distance)**2)))/(n*(2+(rho*Distance))*(11+(rho*Distance)))))))
  m=nls(formula=hill.weir.eq, start=list(rho=rho.start))
  results.m=summary(m)
  cor(Rsquared, predict(m))
  rho.estimate=results.m$parameters[1]
  Distance=sort(Distance)
  exp.rsquared=(((10+(rho.estimate*Distance))/((2+(rho.estimate*Distance))*(11+(rho.estimate*Distance))))*(1+(((3+(rho.estimate*Distance))*(12+(12*(rho.estimate*Distance))+((rho.estimate*Distance)**2)))/(n*(2+(rho.estimate*Distance))*(11+(rho.estimate*Distance))))))
  return(output = data.frame(R_2 = exp.rsquared, Distance = Distance, group = group, points = Rsquared))
}

# Apply the function to each dataset

c1_smooth <- get_LD(c1$V1, c1$V2, 11, 0.001, "1958")
c2_smooth <- get_LD(c2$V1, c2$V2, 11, 0.001, "1963")
c3_smooth <- get_LD(c3$V1, c3$V2, 11, 0.001, "2747")
c4_smooth <- get_LD(c4$V1, c4$V2, 11, 0.001, "2931")
c5_smooth <- get_LD(c5$V1, c5$V2, 11, 0.001, "2932")
c6_smooth <- get_LD(c6$V1, c6$V2, 11, 0.001, "3111")
c7_smooth <- get_LD(c7$V1, c7$V2, 12, 0.001, "4107")
c8_smooth <- get_LD(c8$V1, c8$V2, 10, 0.001, "4117")
c9_smooth <- get_LD(c9$V1, c9$V2, 11, 0.001, "4330")

# Combine all the data frames
df <- rbind(c1_smooth, c2_smooth,c3_smooth,c4_smooth,c5_smooth,c6_smooth,c7_smooth,c8_smooth,c9_smooth)

# Convert distance from base pairs to kilobase pairs
df %>% mutate(Distance = Distance/(10**3)) -> df

# Plot LD decay
library(ggplot2)
p <- ggplot(data = df)  +
  geom_line(aes(x=Distance, y=R_2, color = group),  alpha = 1, size = 1)+
  
  labs(x="Distance (Kilobase)",y=expression(LD~(r^{2})))+
  labs(title = "", col="Group")+
  scale_color_manual(values = my_pal) +
  theme_bw()+
  theme( text = element_text(family = "Times New Roman", size=22),
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.line = element_line(colour = "black"))
p
