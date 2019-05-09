# Made by: Hugo Swenson, 2019-01-25

library(tidyverse)
library(ggplot2)

# Loads the score results as a separate variable, uses bayes score due to better G2M separation
# Arranges these scores with respect to the ordIndex
fumic_stats <- read.csv(file='Fumic_Stats/fumic_stats.csv')
print(fumic_stats)

# Stores each phase score as a new parameter
ref_stats <- fumic_stats$Ref
var_stats <- fumic_stats$Var
ffpe_stats <- fumic_stats$FFPE
perc_stats <- fumic_stats$Perc
change_stats <- fumic_stats$BaseChange
x_ind <- fumic_stats$X

scatter_stats <- data.frame(x_ind, var_stats , ffpe_stats)
bar_stats <- data.frame(change_stats, bar_ind = rep(1, length(change_stats)))
ref_perc_val = ((ffpe_stats/ref_stats)*100)
ref_perc_stats <- data.frame(x_ind, ref_perc_val)
var_perc_val = ((ffpe_stats/var_stats)*100)
var_perc_stats <- data.frame(x_ind, var_perc_val)

pdf('Fumic_Stats/count_vs_var_no.pdf')
# Simple scatter-plot of counts vs var.no
ggplot(scatter_stats, aes(x=x_ind, y=Counts, color=variable)) + 
  geom_point(aes(y=var_stats, col="var_stats")) + 
  geom_point(aes(y=ffpe_stats, col="ffpe_stats"))
dev.off()

pdf('Fumic_Stats/ffpe_vs_ref_perc.pdf')
# Line-plot for the percentual ratio of ffpe/ref  
ggplot(ref_perc_stats, aes(x=x_ind, y=ref_perc_val)) + 
  geom_line(color = '#2d98da') +
  geom_point(color = '#FEA47F')
dev.off()

pdf('Fumic_Stats/ffpe_vs_var_perc.pdf')
# Line-plot for the percentual ratio of ffpe/var  
ggplot(var_perc_stats, aes(x=x_ind, y=var_perc_val)) + 
  geom_line(color = '#fd9644') +
  geom_point(color = '#45aaf2')
dev.off()

pdf('Fumic_Stats/ffpe_vs_var_sca.pdf')
# Simple vs. scatter-plot Var vs. FFPE for the base-change
ggplot(scatter_stats, aes(x=var_stats, y=ffpe_stats)) + 
  geom_point(aes(color=factor(change_stats)))
dev.off()

pdf('Fumic_Stats/ffpe_vs_ref_sca.pdf')
# Simple vs. scatter-plot Ref vs. FFPE for the base-change
ggplot(scatter_stats, aes(x=ref_stats, y=ffpe_stats)) + 
  geom_point(aes(color=factor(change_stats)))
dev.off()

pdf('Fumic_Stats/ffpe_count_bar.pdf')
# Bar plot showing the counts for each FFPE-artefact
ggplot(bar_stats, aes(x=change_stats, y=bar_ind,  fill=factor(change_stats))) + 
  geom_bar(stat="identity")
dev.off()