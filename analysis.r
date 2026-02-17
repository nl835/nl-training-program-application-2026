# ---------------------------------------------------------

# Melbourne Bioinformatics Training Program

# This exercise to assess your familiarity with R and git. Please follow
# the instructions on the README page and link to your repo in your application.
# If you do not link to your repo, your application will be automatically denied.

# Leave all code you used in this R script with comments as appropriate.
# Let us know if you have any questions!


# You can use the resources available on our training website for help:
# Intro to R: https://mbite.org/intro-to-r
# Version Control with Git: https://mbite.org/intro-to-git/

# ----------------------------------------------------------

# Load libraries -------------------
# You may use base R or tidyverse for this exercise

library(tidyverse)

# Load data here ----------------------
# Load each file with a meaningful variable name.
setwd("~/Desktop/nicky/nl-training-program-application-2026/")
metadata <- read.csv("data/GSE60450_filtered_metadata.csv")
cpm_df <- read.csv("data/GSE60450_GeneLevel_Normalized(CPM.and.TMM)_data.csv")


# Inspect the data -------------------------

# What are the dimensions of each data set? (How many rows/columns in each?)
# Keep the code here for each file.

## Expression data
dim(cpm_df)

## Metadata
dim(metadata)

# Prepare/combine the data for plotting ------------------------
# How can you combine this data into one data.frame?
cpm_clean <- cpm_df %>% select(-X)
cpm_long <- cpm_clean %>% pivot_longer(cols = -gene_symbol, 
                                       names_to = "sampleID",
                                       values_to = "cpm")
metadata <- metadata %>% rename(sampleID = X)
merged_long <- cpm_long %>% left_join(metadata, by = "sampleID")
# utilise log CPM for plotting - otherwise impossible to see most datapoints
merged_long <- merged_long %>% mutate(log_cpm = log10(cpm + 1))

# Plot the data --------------------------
## Plot the expression by cell type
## Can use boxplot() or geom_boxplot() in ggplot2
library(ggplot2)

# for readability, rename axes labels
axis_labels <- c(
  "mammary gland, luminal cells, virgin"             = "luminal, virgin",
  "mammary gland, luminal cells, 18.5 day pregnancy" = "luminal, pregnancy",
  "mammary gland, luminal cells, 2 day lactation"    = "luminal, lactation",
  "mammary gland, basal cells, virgin"               = "basal, virgin",
  "mammary gland, basal cells, 18.5 day pregnancy"   = "basal, pregnancy",
  "mammary gland, basal cells, 2 day lactation"      = "basal, lactation"
)


# i've plotted both cpm & log cpm just in case
p <- ggplot(merged_long, aes(x = characteristics, y = cpm, fill = characteristics)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "Cell Type", y = "CPM", fill = "Characteristics", 
       title = "CPM of all genes in mammary gland luminal and basal cells") +
  scale_x_discrete(labels = axis_labels)
p

p_log <- ggplot(merged_long, aes(x = characteristics, y = log_cpm, fill = characteristics)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "Cell Type", y = "Log10(CPM)", fill = "Characteristics", 
       title = "Log10 CPM of all genes in mammary gland luminal and basal cells") +
  scale_x_discrete(labels = axis_labels)
p_log
## Save the plot
### Show code for saving the plot with ggsave() or a similar function
ggsave("results/celltype_CPM_boxplot.pdf", plot = p)
ggsave("results/celltype_logCPM_boxplot.pdf", plot = p_log)
