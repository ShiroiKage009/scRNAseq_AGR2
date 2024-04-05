# This script requires running the file that defines the volcano plot function
library(ggplot2)
library(ggrepel)

source(file = "scripts/volcano plot function.R")

pat_LGR5 = read.csv(file = "C:/Work cache/py_projs/scRNAseq_AGR2/project data cache/testing integration with separation and the stem cells part 2/saved files/pat_LGR5.csv", row.names = 1)
cont_LGR5 = read.csv(file = "C:/Work cache/py_projs/scRNAseq_AGR2/project data cache/testing integration with separation and the stem cells part 2/saved files/cont_LGR5.csv", row.names = 1)
enterocytes = read.csv(file = "C:/Work cache/py_projs/scRNAseq_AGR2/project data cache/testing integration with separation and the stem cells part 2/saved files/enterocytes_df.csv", row.names = 1)

pat_LGR5_procant = read.csv(file = "C:/Work cache/py_projs/scRNAseq_AGR2/project data cache/testing integration with separation and the stem cells part 2/saved files/pat_LGR5_procant.csv", row.names = 1)
cont_LGR5_procant = read.csv(file = "C:/Work cache/py_projs/scRNAseq_AGR2/project data cache/testing integration with separation and the stem cells part 2/saved files/cont_LGR5_procant.csv", row.names = 1)
enterocytes_procant = read.csv(file = "C:/Work cache/py_projs/scRNAseq_AGR2/project data cache/testing integration with separation and the stem cells part 2/saved files/enterocytes_procant_df.csv", row.names = 1)


geom_text_repel(max.overlaps = 100)

plot_volcano(df = pat_LGR5, 
             x_axis = "Log.Fold.Change", 
             pval_col = "Adjusted.P.Value", 
             top_n = 160, 
             plot_title = "Gated Patient LGR5 cluster vs rest (Wilcoxon)")

plot_volcano(df = cont_LGR5, 
             x_axis = "Log.Fold.Change", 
             pval_col = "Adjusted.P.Value", 
             top_n = 160, 
             plot_title = "Gated Control LGR5 cluster vs rest (Wilcoxon)")

plot_volcano(df = enterocytes, 
             x_axis = "Log.Fold.Change", 
             pval_col = "Adjusted.P.Value", 
             top_n = 160, 
             plot_title = "Gated Metaplastic Enterocytes vs rest (Wilcoxon)")

plot_volcano(df = pat_LGR5_procant, 
             x_axis = "Log.Fold.Change", 
             pval_col = "Adjusted.P.Value", 
             top_n = 160, 
             plot_title = "Patient LGR5 cluster vs rest (Wilcoxon)")

plot_volcano(df = cont_LGR5_procant, 
             x_axis = "Log.Fold.Change", 
             pval_col = "Adjusted.P.Value", 
             top_n = 160, 
             plot_title = "Control LGR5 cluster vs rest (Wilcoxon)")

plot_volcano(df = enterocytes_procant, 
             x_axis = "Log.Fold.Change", 
             pval_col = "Adjusted.P.Value", 
             top_n = 160, 
             plot_title = "Metaplastic Enterocytes vs rest (Wilcoxon)")
