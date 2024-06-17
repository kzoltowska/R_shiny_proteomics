# get all packages
packages <- c(
  "shiny", "ggplot2", "pheatmap", "readr", "visdat", "reshape2", "dplyr", "tidyr",
  "ggpubr", "RColorBrewer", "imputeLCMD", "rrcovNA", "standby", "missForest", "DT",
  "shinyWidgets", "shinycssloaders", "tibble", "stringr", "name"
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


# function to preprocess data based on selected conditions
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
