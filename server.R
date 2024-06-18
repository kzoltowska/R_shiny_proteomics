server <- function(input, output, session) {
  data_init <- reactive({
    req(input$file1)
    req(input$file2)
    infile <- input$file1
    read.delim2(infile$datapath, header = TRUE)
  })

  sample_info <- reactive({
    req(input$file2)
    infile2 <- input$file2
    read.delim2(infile2$datapath, header = TRUE)
  })

  observe({
    updateSelectInput(session,
      inputId = "group1", label = NULL,
      choices = sample_info()$Sample_name,
      selected = NULL
    )
  })
  observe({
    updateSelectInput(session,
      inputId = "group2", label = NULL,
      choices = sample_info()$Sample_name,
      selected = NULL
    )
  })
  observe({
    updateSelectInput(session,
      inputId = "group1_stat", label = NULL,
      choices = sample_info()$Sample_name,
      selected = NULL
    )
  })
  observe({
    updateSelectInput(session,
      inputId = "group2_stat", label = NULL,
      choices = sample_info()$Sample_name,
      selected = NULL
    )
  })

  groups_color <- reactive({
    groups <- as.factor(sample_info()$Condition)
    groups_color_tmp <- groups
    levels(groups_color_tmp) <- brewer.pal(length(levels(groups)), "Dark2")
    groups_color_tmp
  })

  data <- reactive({
    data_init2 <- data_init()[, c(1, 3, 4, 7:18)] %>% mutate_all(~ ifelse(is.nan(.), NA, .))
    colnames(data_init2)[colnames(data_init2) == "PG.ProteinGroups"] <- "uniprot"
    colnames(data_init2)[colnames(data_init2) == "PG.Genes"] <- "gene"
    colnames(data_init2)[colnames(data_init2) == "PG.ProteinDescriptions"] <- "full_name"
    for (i in 1:nrow(sample_info())) {
      for (j in 1:ncol(data_init2)) {
        if (colnames(data_init2)[j] == sample_info()$Sample_col_name[i]) {
          colnames(data_init2)[j] <- sample_info()$Sample_name[i]
        }
      }
    }
    if (input$Keep_one == TRUE) {
      position <- as.vector(str_split(data_init2$uniprot, pattern = ";"))
      gene_1 <- as.vector(str_split(data_init2$gene, pattern = ";"))
      full_1 <- as.vector(str_split(data_init2$full_name, pattern = ";"))
      pos <- c()
      for (i in c(1:length(position))) {
        pos <- c(pos, match(min(nchar(position[[i]])), as.vector(nchar(position[[i]]))))
      }
      uniprot <- c()
      for (i in c(1:length(pos))) {{ j <- pos[i]
        uniprot <- c(uniprot, position[[i]][j]) }}
      gene <- c()
      for (i in c(1:length(pos))) {{ j <- pos[i]
        if (length(gene_1[i]) >= j) {
          gene <- c(gene, gene_1[[i]][j])
        } else {
          gene <- c(gene, gene_1[[i]][1])
        } }}
      data_init2$uniprot <- uniprot
      data_init2$gene <- gene
      full_name <- c()
      for (i in c(1:length(pos))) {{ j <- pos[i]
        full_name <- c(full_name, full_1[[i]][j]) }}
      data_init2$full_name <- full_name
    } else {
      data_init2 <- data_init2
    }
    if (input$filter == "") {
      data_init2 <- data_init2
    } else {
      i <- as.character(input$filter)
      data_init2 <- data_init2[!grepl(i, data_init2$gene), ]
    }
    data_init2
  })



  data_long <- reactive({
    df_long <- data() %>% pivot_longer(
      names_to = "sample", values_to = "intensity",
      cols = sample_info()$Sample_name
    )
    group_col <- c()
    for (j in (1:nrow(df_long))) {
      for (i in (1:nrow(sample_info()))) {
        if (df_long$sample[j] == sample_info()$Sample_name[i]) {
          group_col <- c(group_col, sample_info()$Condition[i])
        }
      }
    }
    df_long$group <- group_col
    df_long
  })


  # create entire data table
  output$entire_data <- DT::renderDT({
    req(data())
    data() %>%
      mutate_if(is.numeric, round, digits = 2)
  })

  # create summary table
  output$summary_table <- DT::renderDT({
    req(data())
    df <- as.data.frame(summary(data()))
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
    req(data())
    vis_dat(data()[, sample_info()$Sample_name]) +
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

  output$barplot <- renderPlot({
    req(data())
    p <- barplot(colSums(!is.na(data() %>% dplyr::select(-c(gene, full_name, uniprot)))),
      col = as.character(groups_color()), las = 2,
      ylim = c(0, max(colSums(!is.na(data() %>% dplyr::select(-c(gene, full_name, uniprot))))) * 1.5)
    )
    text(
      x = p, y = 20 + colSums(!is.na(data() %>% dplyr::select(-c(gene, full_name, uniprot)))),
      labels = colSums(!is.na(data() %>% dplyr::select(-c(gene, full_name, uniprot))))
    )
  })

  mat_imp <- reactive({
    mat_tmp <- as.matrix(data()[, c(sample_info()$Sample_name)])
    rownames(mat_tmp) <- data()$gene
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
      mat_imp_all <- preprocessing(mat = mat_tmp, log = log, norm = norm, imp = imp)
    } else {
      mat1 <- mat_tmp[, c(input$group1)]
      mat_imp1 <- preprocessing(mat = mat1, log = log, norm = norm, imp = imp)
      mat2 <- mat_tmp[, c(input$group2)]
      mat_imp2 <- preprocessing(mat = mat2, log = log, norm = norm, imp = imp)
      mat_imp_all <- cbind(mat_imp1, mat_imp2)
    }
    mat_imp_all
  })

  output$boxplot_data <- renderPlot({
    req(input$file1)
    req(input$file2)
    boxplot(mat_imp(), col = as.character(groups_color()), las = 2)
  })

  output$heatmap_data <- renderPlot({
    req(input$file1)
    req(input$file2)
    pheatmap(mat_imp(),
      cluster_cols = TRUE, cluster_rows = FALSE,
      scale = "row", show_rownames = FALSE
    )
  })

  table <- reactive({
    req(input$file1)
    req(input$file2)
    cbind(gene = data()$gene, uniprot = data()$uniprot, round(mat_imp(), 2)) %>%
      as.data.frame() %>%
      remove_rownames()
  })

  output$matrix_log_norm_imp <- DT::renderDT({
    req(input$file1)
    req(input$file2)
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



  stat_tmp1 <- reactive({
    if (input$stat == "ROTS") {
      req(input$file1)
      req(input$file2)
      validate(need(
        max(rowSums(is.na(mat_imp()))) < 2,
        "Data contains more than two missing values per row"
      ))
      req(length(input$group1_stat) > 1 & length(input$group2_stat) > 1)

      ROTS_tmp1 <- ROTS(mat_imp()[, c(input$group1_stat, input$group2_stat)],
        groups = c(rep(1, length(input$group1_stat)), rep(0, length(input$group2_stat)))
      )
      ROTS_tmp1 <- as.data.frame(cbind(
        gene = rownames(ROTS_tmp1$data),
        logfc = ROTS_tmp1$logfc,
        pvalue = ROTS_tmp1$pvalue, FDR = ROTS_tmp1$FDR, ROTS_tmp1$data
      )) %>%
        remove_rownames() %>%
        mutate(across(2:last_col(), ~ as.numeric(.)))
      ROTS_tmp1
    } else if (input$stat == "limma") {
      req(length(input$group1_stat) > 1 & length(input$group2_stat) > 1)
      data_raw <- data() %>% rename_with(~ paste(., "_raw", sep = "_"))
      condition <- c(rep("group1", length(input$group1_stat)), rep("group2", length(input$group2_stat)))
      design <- model.matrix(~ 0 + condition)
      colnames(design) <- gsub("condition", "", colnames(design))
      fit1 <- lmFit(mat_imp()[, c(input$group1_stat, input$group2_stat)], design = design)
      cont <- makeContrasts(group2 - group1, levels = design)
      fit2 <- contrasts.fit(fit1, contrasts = cont)
      fit3 <- eBayes(fit2)
      limma.results <- topTable(fit3, adjust = "BH", n = Inf)
      colnames(limma.results) <- c("logfc", "AveExpr", "t", "pvalue", "FDR", "B")
      limma.results <- cbind(
        gene = rownames(limma.results),
        limma.results[, c("logfc", "pvalue", "FDR")]
      )
      rownames(limma.results) <- limma.results$gene
      limma.results <- merge(limma.results, mat_imp(), by = "row.names")
      limma.results <- subset(limma.results, select = -c(Row.names))
      limma.results
    }
  })

  stat <- reactive({
    merge(data()[, c("gene", "uniprot")], stat_tmp1(), by = "gene")
  })

  output$volcano <- renderPlot({
    req(stat())
    pCutoff <- input$slider_p
    FCcutoff <- input$slider_fc

    if (input$pval == "p value") {
      EnhancedVolcano(stat(),
        lab = stat()$gene,
        x = "logfc", y = "pvalue", pCutoff = pCutoff, FCcutoff = FCcutoff
      )
    } else if (input$pval == "adjusted p value / FDR") {
      EnhancedVolcano(stat(),
        lab = stat()$gene,
        x = "logfc", y = "FDR", pCutoff = pCutoff, FCcutoff = FCcutoff
      )
    }
  })
  
stat_subset<-reactive({
  stat() %>%
    filter((logfc<input$logfc_heat_tab_min | logfc>input$logfc_heat_tab_max) & FDR<=input$pval_heat_tab)
})  

output$heatmap_stat<-renderPlot({
  req(stat())
  heat_stat_mat<-as.matrix(stat_subset() %>% dplyr::select(c(input$group1_stat, input$group2_stat)) )
  rownames(heat_stat_mat)<-stat_subset()$gene
  validate(need(nrow(heat_stat_mat)>=2, "Nothing to display"))
  pheatmap(heat_stat_mat, show_rownames = input$Rownames_heat, cluster_cols=TRUE, 
           cluster_rows = TRUE, scale="row")
})



  output$stat_table <- DT::renderDT({
    req(stat_subset())
    stat_tmp1 <- stat_subset() %>%
      remove_rownames() %>%
      mutate_if(is.numeric, round, 2)
  })

  output$downloadData_stat <- downloadHandler(
    filename = function() {
      "stat.csv"
    },
    content = function(fname) {
      write.csv(stat_subset(), fname, row.names = FALSE)
    }
  )
}
