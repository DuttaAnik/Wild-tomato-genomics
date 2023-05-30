# Anik Dutta
# Script for comparing whether 9 populations differ significantly from each other in terms of various population genetics statistics.

# Load necessary libraries
library(dplyr)
library(lattice)
library(survival)
library(Formula)
library(Hmisc)
library(gmodels)
library(Rmisc)
library(tidyverse)

#Loading Pixy data

allHaps_pi <- read.csv("Chilense_Pixy_output.out_pi.txt", sep = "\t")
names(allHaps_pi)

allHaps_pi$Pi<-(allHaps_pi$count_diffs)/(allHaps_pi$count_comparisons)

group_by(allHaps_pi, pop) %>% summarise(tot = sum(count_diffs, na.rm = T)/sum(count_comparisons, na.rm = T))

##Significance testing
kruskal.test(Pi~ as.factor(pop), data = allHaps_pi)
pairwise.wilcox.test(allHaps_pi$Pi, allHaps_pi$pop,
                     p.adjust.method = "bonferroni")

### When used the TajimaD R script, you can use this one to estimate bootstrapped confidence interval 
data %>% 
  select(TajimaD, Group) %>% 
  group_by(Group) %>% 
  dplyr::summarise(data = list(smean.cl.boot(cur_data(), conf.int = .95, B = 1000, na.rm = TRUE))) %>%
  tidyr::unnest_wider(data)

#group_map(~ smean.cl.boot(., conf.int = .95, B = 1000, na.rm = TRUE)) %>%
  #bind_rows()

# Confidence interval for Fst values

filt_Haps_Fst <- read.csv("Chilense_Pixy_output.out_fst.txt", sep = "\t")
names(filt_Haps_Fst)
filt_Haps_Fst$avg_wc_fst[filt_Haps_Fst$avg_wc_fst < 0] <- 0
allfst_sum <- group_by(filt_Haps_Fst, pop1, pop2) %>% summarise(tot = mean(avg_wc_fst, na.rm = T))

allFst_sum<-filt_Haps_Fst %>%
  select(avg_wc_fst, pop1, pop2) %>%
  group_by(pop1, pop2) %>%
  dplyr::summarise(data = list(smean.cl.boot(cur_data(), conf.int = .95, B = 1000, na.rm = TRUE))) %>%
  tidyr::unnest_wider(data)

# Confidence interval for Dxy values (absolute divergence measure)
filt_Haps_DXY <- read.csv("Chilense_Pixy_output.out_dxy.txt", sep = "\t")

allDXY_sum <- group_by(filt_Haps_DXY, pop1, pop2) %>% summarise(tot = sum(count_diffs, na.rm = T)/sum(count_comparisons, na.rm = T))

allDXY_sum<-filt_Haps_DXY %>%
  select(avg_dxy, pop1, pop2) %>%
  group_by(pop1, pop2) %>%
  dplyr::summarise(data = list(smean.cl.boot(cur_data(), conf.int = .95, B = 1000, na.rm = TRUE))) %>%
  tidyr::unnest_wider(data)
write.csv(allDXY_sum, file = "DXY.csv")

DXY_sd <- aggregate(filt_Haps_DXY$avg_dxy,list(filt_Haps_DXY$pop1,filt_Haps_DXY$pop2),sd,na.rm=T)
