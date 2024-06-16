# get all packages
packages <- c(
  "shiny", "ggplot2", "pheatmap", "readr", "visdat", "reshape2", "dplyr", "tidyr",
  "ggpubr", "RColorBrewer", "imputeLCMD", "rrcovNA", "standby", "missForest", "DT",
  "shinyWidgets", "shinycssloaders", "tibble"
)
packages_bioc <- c("limma", "ROTS", "EnhancedVolcano", "pcaMethods")
for (i in packages) {
  if (!require(i, character.only = TRUE)) {
    install.packages(i, dependencies = TRUE)
    library(i, character.only = TRUE)
  }
}
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(packages_bioc)
for (i in packages_bioc) {
  library(i, character.only = TRUE)
}



path_data <- "~/Proteomics_shiny/20240422_130246_MCF001701_Report_pivotProteinGroups_noCrossRunNormalization_short.txt"
data <- read.delim2(path_data)
data <- data[, c(1, 3, 4, 7:18)]
path_sample <- "~/Proteomics_shiny/sample_info.txt"
sample_info <- read.delim2(path_sample)
data <- data %>% mutate_all(~ ifelse(is.nan(.), NA, .))
groups <- as.factor(sample_info$Condition)
groups_color <- groups
levels(groups_color) <- brewer.pal(length(levels(groups)), "Dark2")

colnames(data)[colnames(data) == "PG.ProteinGroups"] <- "uniprot"
colnames(data)[colnames(data) == "PG.Genes"] <- "gene"
colnames(data)[colnames(data) == "PG.ProteinDescriptions"] <- "full_name"

for (i in 1:nrow(sample_info)) {
  for (j in 1:ncol(data)) {
    if (colnames(data)[j] == sample_info$Sample_col_name[i]) {
      colnames(data)[j] <- sample_info$Sample_name[i]
    }
  }
}

# creating long dataframe for input to ggplot
data_long <- data %>% pivot_longer(names_to = "sample", values_to = "intensity", cols = sample_info$Sample_name)
group_col <- c()
for (j in (1:nrow(data_long))) {
  for (i in (1:nrow(sample_info))) {
    if (data_long$sample[j] == sample_info$Sample_name[i]) {
      group_col <- c(group_col, sample_info$Condition[i])
    }
  }
}
data_long$group <- group_col

# creating matrix for the analysis
mat <- as.matrix(data[, c(sample_info$Sample_name)])
rownames(mat)<-data$gene

# function to preprocess data based on selected condtions
preprocessing <- function(mat, log, norm, imp) {
  if (log == TRUE) {
    mat_log <- log2(mat)
  } else {
    mat_log <- mat
  }
  
  if (norm == "None") {
    mat_norm <- mat_log
  } else if (norm == "Quantile") {
    mat_norm <- limma::normalizeQuantiles(mat_log)
  } else if (norm == "Cloess") {
    mat_norm <- normalizeCyclicLoess(mat_log, method = "fast")
  }
  
  if (imp == "None") {
    mat_imp <- mat_norm
  } else if (imp == "MinDet") {
    mat_imp <- impute.MinDet(mat_norm)
  } else if (imp == "ImpSeq") {
    mat_imp <- impSeq(mat_norm)
  } else if (imp == "MissForest") {
    mat_imp_mF <- missForest(mat_norm, maxiter = 1)
    mat_imp <- mat_imp_mF$ximp
  } else if (imp == "MinProb") {
    mat_imp <- impute.MinProb(mat_norm)
  }
  return(mat_imp)
}