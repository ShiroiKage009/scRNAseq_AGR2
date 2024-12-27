# Install clustree if you haven't already
# install.packages("clustree")
library(clustree)

# Read the CSV file
clustree_data <- read.csv("S:/data cache/code_in_out/agr2/script_4_out/clustree_prep_full_ant.csv")

# Use gsub to replace underscores in the numeric part of the column names
colnames(clustree_data) <- gsub("leiden_(\\d+)_(\\d+)", "leiden_\\1.\\2", colnames(clustree_data))


# Generate the clustree plot
clustree(clustree_data, prefix = "leiden_")

