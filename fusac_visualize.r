# Made by: Hugo Swenson, 2019-01-25

library(tidyverse)
library(ggplot2)

# Loads the score results as a separate variable, uses bayes score due to better G2M separation
# Arranges these scores with respect to the ordIndex
fusac_stats <- read.csv(file='FUSAC_Stats/fusac_stats.csv')
print(fusac_stats)

# Stores each phase score as a new parameter
ref_stats <- fusac_stats$Ref
var_stats <- fusac_stats$Var
ffpe_stats <- fusac_stats$FFPE
perc_stats <- fusac_stats$Perc
change_stats <- fusac_stats$BaseChange
x_ind <- fusac_stats$X

scatter_stats <- data.frame(x_ind, var_stats , ffpe_stats)
bar_stats <- data.frame(change_stats, bar_ind = rep(1, length(change_stats)))
perc_stats_amp <- perc_stats*100
line_perc_stats <- data.frame(x_ind, perc_stats)

ref_perc_val = ((ffpe_stats/ref_stats)*100)
ref_perc_stats <- data.frame(x_ind, ref_perc_val)
var_perc_val = ((ffpe_stats/var_stats)*100)
var_perc_stats <- data.frame(x_ind, var_perc_val)

pdf('FUSAC_Stats/ffpe_vs_var_sca.pdf')
# Simple vs. scatter-plot Var vs. FFPE for the base-change
ggplot(scatter_stats, aes(x=var_stats, y=ffpe_stats)) + 
  geom_point(aes(color=factor(change_stats)))
dev.off()

pdf('FUSAC_Stats/ffpe_vs_ref_sca.pdf')
# Simple vs. scatter-plot Ref vs. FFPE for the base-change
ggplot(scatter_stats, aes(x=ref_stats, y=ffpe_stats)) + 
  geom_point(aes(color=factor(change_stats)))
dev.off()

pdf('FUSAC_Stats/ffpe_count_bar.pdf')
# Bar plot showing the counts for each FFPE-artefact
ggplot(bar_stats, aes(x=change_stats, y=bar_ind,  fill=factor(change_stats))) + 
  geom_bar(stat="identity")
dev.off()

pdf('FUSAC_Stats/ffpe_vs_var_bar.pdf')
# Stacked Bar Plot with Colors and Legend
bar_counts <- table(fusac_stats$Var, fusac_stats$FFPE)
barplot(bar_counts,  x_lab="Variant-record no.", y_lab="Counts"  col=c("#fd9644","#45aaf2"),
  legend = rownames(counts))
dev.off()

pdf('FUSAC_Stats/perc_rat_ffpe.pdf')
# Line-plot for the percentual ratio of ffpe/total-counts 
ggplot(line_perc_stats, aes(x=x_ind, y=line_perc_stats)) + 
  geom_line(color = '#fd9644') +
  geom_point(color = '#45aaf2')
dev.off()
