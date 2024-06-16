#
# This is the server logic of a Shiny web application.

# path to the data
path_data <- "~/Proteomics_shiny/20240422_130246_MCF001701_Report_pivotProteinGroups_noCrossRunNormalization.txt"
data <- read.delim2(path_data)
data <- data[, c(1, 3, 4, 7:18)]
path_sample <- "~/Proteomics_shiny/sample_info.txt"
sample_info <- read_delim(path_sample, delim = "\t", escape_double = FALSE, trim_ws = TRUE)
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

server <- function(input, output) {
  # create entire data table
  output$entire_data <- DT::renderDT({
    data %>%
      mutate_if(is.numeric, round, digits = 2)
  })

  # create summary table
  output$summary_table <- DT::renderDT({
    df <- as.data.frame(summary(data))
    df <- df[, c(2, 3)]
    colnames(df) <- c("Variable", "info")
    df$Variable <- as.factor(df$Variable)
    ls <- list()
    levels <- levels(df$Variable)
    for (i in levels) {
      ls[[i]] <- subset(df, df[, 1] == i, select = c("info"))
    }
    df2 <- as.data.frame(ls)
    colnames(df2) <- levels
    df2[is.na(df2) == TRUE] <- ""
    df2
  })

  # create summary plot
  output$summary_plot <- renderPlot({
    vis_dat(data[, sample_info$Sample_name]) +
      theme(
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 15),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        plot.title = element_text(size = 16, face = "bold.italic")
      ) +
      ggtitle("Graphical summary of the data")
  })

  output$boxplot_raw <- renderPlot({
    ggplot(data_long) +
      geom_boxplot(aes(x = sample, y = intensity, fill = group)) +
      theme_pubr() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      scale_fill_manual(values = levels(groups_color))
  })

  output$boxplot_log <- renderPlot({
    data_long$intensity <- log2(data_long$intensity)
    ggplot(data_long) +
      geom_boxplot(aes(x = sample, y = intensity, fill = group)) +
      theme_pubr() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      scale_fill_manual(values = levels(groups_color))
  })

  mat_imp <- reactive({
    if (input$log == TRUE) {
      log <- TRUE
    } else {
      log <- FALSE
    }
    if (input$normalisation == "None") {
      norm <- "None"
    } else if (input$normalisation == "Quantile") {
      norm <- "Quantile"
    } else if (input$normalisation == "Cloess") {
      norm <- "Cloess"
    }
    if (input$imputation == "None") {
      imp <- "None"
    } else if (input$imputation == "MinDet") {
      imp <- "MinDet"
    } else if (input$imputation == "ImpSeq") {
      imp <- "ImpSeq"
    } else if (input$imputation == "MissForest") {
      imp <- "MissForest"
    } else if (input$imputation == "MinProb") {
      imp <- "MinProb"
    }
    if (is.null(input$group1) == TRUE | is.null(input$group2) == TRUE |
      length(input$group1) == 1 | length(input$group2) == 1) {
      mat_imp <- preprocessing(mat = mat, log = log, norm = norm, imp = imp)
    } else {
      mat1 <- mat[, c(input$group1)]
      mat_imp1 <- preprocessing(mat = mat1, log = log, norm = norm, imp = imp)
      mat2 <- mat[, c(input$group2)]
      mat_imp2 <- preprocessing(mat = mat2, log = log, norm = norm, imp = imp)
      mat_imp <- cbind(mat_imp1, mat_imp2)
    }
  })

  output$boxplot_data <- renderPlot({
    boxplot(mat_imp(), col = as.character(groups_color), las = 2)
  })

  output$heatmap_data <- renderPlot({
    pheatmap(mat_imp(), cluster_cols = TRUE, cluster_rows = FALSE, scale = "row", show_rownames = FALSE)
  })

  table<-reactive({
    cbind(gene=data$gene, uniprot = data$uniprot, round(mat_imp(), 2)) %>%
      as.data.frame() %>% remove_rownames()
  })
  
  output$matrix_log_norm_imp <- DT::renderDT({
  table()
  })

  output$downloadData <- downloadHandler(
    filename = function() {
      "pre_preprocessed.csv"
    },
    content = function(fname) {
      write.csv(table(), fname, row.names = FALSE)
    }
  )
}
