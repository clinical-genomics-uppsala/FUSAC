# Made by: Hugo Swenson, 2019-01-25

library(tidyverse)
library(ggplot2)

# Loads the score results as a separate variable, uses bayes score due to better G2M separation
# Arranges these scores with respect to the ordIndex
fumic_stats <- read.csv(file='Fumic_Stats/fumic_stats.csv')
print(fumic_stats)

# Stores each phase score as a new parameter
var_stats <- fumic_stats$Var
ffpe_stats <- fumic_stats$FFPE
perc_stats <- fumic_stats$Perc
change_stats <- fumic_stats$BaseChange
x_ind <- fumic_stats$X

scatter_stats <- data.frame(x_ind, var_stats , ffpe_stats)
bar_stats <- data.frame(change_stats, bar_ind = rep(1, length(change_stats)))

# Simple scatter-plot of counts vs var.no
ggplot(scatter_stats, aes(x=x_ind, y=Counts, color=variable)) + 
  geom_point(aes(y=var_stats, col="var_stats")) + 
  geom_point(aes(y=ffpe_stats, col="ffpe_stats"))

# Simple vs. scatter-plot Var vs. FFPE for the base-change
ggplot(scatter_stats, aes(x=var_stats, y=ffpe_stats)) + 
  geom_point(aes(color=factor(change_stats)))

# Simple vs. scatter-plot Var vs. FFPE for the base-change
ggplot(scatter_stats, aes(x=var_stats, y=ffpe_stats)) + 
  geom_point(aes(color=factor(change_stats)))+
  xlim(0, 50) +
  ylim(0, 50)

# Simple vs. scatter-plot Var vs. FFPE for the base-change
ggplot(bar_stats, aes(x=change_stats, y=bar_ind,  fill=factor(change_stats))) + 
  geom_bar(stat="identity")