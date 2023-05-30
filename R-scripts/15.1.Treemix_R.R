# Anik Dutta
# Plotting treemix results like a tree. The functions used for plotting the treemix outputs are available in the treemix download package.
# It depends on how many migration events you want to plot. So, first check what is the optimum number of migration events using the optm package. Then plot that specific migration numbered output.

source("~/Documents/treemix-1.13/src/plotting_funcs.R")

pdf(file="Migration_5.pdf")
plot_tree("Chilense5") # treemix results with 5 migrations events Chilense5
#plot_tree("Chilense.1.3")
dev.off()

plot_resid("Chilense0",pop_order = "Pop_list.txt") # plot residues with 0 migration events

plot_resid("Chilense4",pop_order = "Pop_list.txt") # plot residues with 4 migration events

# Plotting likelihood of migration events
library(BioGeoBEARS)

get_f("Chilense0")
get_f("Chilense1")
get_f("Chilense2")
get_f("Chilense3")
get_f("Chilense4")
get_f("Chilense5")
get_f("Chilense6")
get_f("Chilense7")
get_f("Chilense8")
get_f("Chilense9")
get_f("Chilense10")
plot_resid("Chilense10", pop_order = "Pop_list.txt")
prefix="Chilense"

library(R.utils)

#lrttest( 376.67 ,368.451,9,8, returnwhat="all")


a<-read.csv("Likelihood_var_treemix.csv", header = T)
names(a)

library(ggplot2)

a$Migration <- factor(a$Migration, levels = a$Migration)

pdf("Likelihood_treemix.pdf", height = 4, width = 6, useDingbats = F)

ggplot(data = a, aes(x = Migration, y = Variance,group = 1)) + geom_point(size = 3) + geom_line(lty = 3) + 
  xlab("Migration event") + ylab("Variance explained") + theme(axis.title = element_text(size = 16),axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(),panel.border = element_rect(color = "black",
                                                                     fill = NA,
                                                                     size = 1))+xlab("Migration event") + ylab("Variance explained")


dev.off()

