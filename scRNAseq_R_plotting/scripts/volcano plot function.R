# Load the libraries ggrepel and ggplot2
# This function assumes that the gene names are row names
# x_axis is the string name of the column that contains the fold change data
# y_axis is automatically assigned as MinusLog10P calculated out of pval_col
# highlight_gene is a string that looks to highlight the gene specified by the symbol
# top_n specifies the number of most variable genes to be labeled determined by p value
# plot_title is a string that gives the plot its title

plot_volcano <- function(df, x_axis, pval_col, highlight_gene = NULL, top_n = 0, plot_title = "Volcano Plot") {
  
  y_axis = "MinusLog10P"
  # Add row names as a new column for easy access
  df$GeneName <- rownames(df)
  # Add a -log10 column
  df$MinusLog10P = -log10(df[[pval_col]])
  
  # Compute significance
  df$Significance <- ifelse(df[[pval_col]] < 0.05 & abs(df[[x_axis]]) > 1, ifelse(df[[x_axis]] < 0, "significant negative", "significant"), "not significant")
  
  # Start plotting
  p <- ggplot(df, aes_string(x = x_axis, y = y_axis)) +
    geom_point(aes_string(color = "Significance"), alpha = 0.5) +
    scale_color_manual(values = c("significant negative" = "#2C56FF", "significant" = "red", "not significant" = "grey")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkgrey") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "darkgrey") +
    annotate("text", x = Inf, y = -log10(0.05) + 0.5, label = "P-value Threshold", color = "black", hjust = 1) +
    labs(title = plot_title, x = "Log Fold Change", y = "-log10(Adjusted P-value)")
  
  # Label top n most significant genes
  if (top_n > 0) {
    top_genes <- df[order(df[[pval_col]]), ][1:top_n, ]
    p <- p + geom_text_repel(data = top_genes,
                             aes_string(x = x_axis, y = y_axis, label = "GeneName"),
                             size = 3, 
                             vjust = 1.5, 
                             hjust = 0.5)
  }
  
  # Highlight a specific gene if provided
  if (!is.null(highlight_gene) && highlight_gene %in% df$GeneName) {
    p <- p + 
      geom_point(data = subset(df, GeneName == highlight_gene), aes_string(x = x_axis, y = y_axis), color = "#18B728", size = 3, shape = 1, stroke = 2) +
      geom_label(data = subset(df, GeneName == highlight_gene), aes_string(x = x_axis, y = y_axis, label = sprintf("'%s'", highlight_gene)), hjust = 0.5, vjust = -1.5, label.size = NA, color = "#18B728")
  }
  
  # Return the plot
  return(p)
}

# Example of use:
# library(ggplot2)
# library(ggrepel)
# metaplastic_ant = read.csv(file = "C:/Work cache/py_projs/scRNAseq_AGR2/project data cache/testing integration with separation and the stem cells part 2/saved files/metaplastic_ant.csv", row.names = 1)
# plot_volcano(df = metaplastic_ant, x_axis = "Log.Fold.Change", pval_col = "Adjusted.P.Value", highlight_gene = 'AQP1', top_n = 40, plot_title = "Metaplastic LGR5+ cells vs Gastric LGR5+ cells (Wilcoxon)")
