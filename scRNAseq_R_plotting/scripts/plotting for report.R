library(ggplot2)
library(ggrepel)

source(file = "scripts/volcano plot function.R")
pat_ant = read.csv(file = 'C:/Work cache/py_projs/scRNAseq_AGR2/project data cache/testing integration with separation and the stem cells part 2/saved files/metaplastic_ant.csv', row.names = 1)

geom_text_repel(max.overlaps = 100)

plot_volcano(df = pat_ant, x_axis = "Log.Fold.Change", pval_col = "Adjusted.P.Value", top_n = 150, plot_title = "Metaplastic LGR5+ cells vs Gastric LGR5+ cells (Wilcoxon)")

