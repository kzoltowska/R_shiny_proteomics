#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Author: Katarzyna Zoltowska
# Aim: Testing proteomics data analysis workflows
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


ui <- fluidPage(
  # Application title
  titlePanel("Proteomics data analysis"),
  h2("Hobby project"),
  # create tabset
  tabsetPanel(
    # create first tab
    tabPanel(
      title = "Data overview", fluid = TRUE,
      sidebarLayout(
        sidebarPanel(
          width = 2,
          # create radiobuttons to select how to view the data
          radioButtons(
            inputId = "Description",
            label = "How would you like to view the data?",
            choices = c("Entire dataset", "Tabular summary", "Graphical summary", "Boxplot - raw intensities", "Boxplot - log2 intensities")
          )
        ),
        mainPanel(
          width = 10,
          br(),
          # show entire data
          conditionalPanel(
            condition = "input.Description == 'Entire dataset'",
            withSpinner(DTOutput("entire_data"))
          ),
          # Show a summary table
          conditionalPanel(
            condition = "input.Description == 'Tabular summary'",
            DTOutput("summary_table")
          ),
          # Show a graphical summary
          conditionalPanel(
            condition = "input.Description == 'Graphical summary'",
            plotOutput("summary_plot", height = "600px")
          ),
          # Show boxplot with raw values
          conditionalPanel(
            condition = "input.Description == 'Boxplot - raw intensities'",
            plotOutput("boxplot_raw", height = "600px")
          ),
          # Show boxplot with log values
          conditionalPanel(
            condition = "input.Description == 'Boxplot - log2 intensities'",
            plotOutput("boxplot_log", height = "600px")
          )
        )
      )
    ),
    # create tab for log, normalisation and imputation
    tabPanel(
      title = "Data preprocessing",
      fluidRow(
        br(),
        column(
          2,
          checkboxInput(
            inputId = "log",
            label = HTML("<b>Convert to log2</b>")
          )
        ),
        column(
          3,
          pickerInput(
            inputId = "normalisation",
            label = "Choose normalisation type",
            choices = c("None", "Quantile", "Cloess")
          )
        ),
        column(
          3,
          pickerInput(
            inputId = "imputation",
            label = "Choose imputation type",
            choices = c("None", "MinDet", "ImpSeq", "MissForest", "MinProb"),
            choicesOpt = list(subtext = c("", "", "", "Very slow (maxiter=1)", ""))
          )
        )
      ),
      br(),
      h5("By default the normalisation and imputation is performed on all samples.\n
         If this should be done by condition, please put selected samples (at least two per group) into respective groups.\n"),
      fluidRow(
        column(
          3,
          selectInput(inputId = "group1", choices = sample_info$Sample_name,
                      label = "Group_1", multiple = TRUE, selected = NULL)
        ),
        column(3, selectInput(inputId = "group2", choices = sample_info$Sample_name,
                              label = "Group_2", multiple = TRUE, selected = NULL))
      ),
      br(),
      # output in the form of boxplot and heatmap
      fluidRow(
        column(6, withSpinner(plotOutput("boxplot_data"))),
        column(6, withSpinner(plotOutput("heatmap_data")))
      ),
      br(),
      # output in the form of table than can be downloaded and saved as csv file
      fluidRow(DTOutput("matrix_log_norm_imp")),
      fluidRow(
      downloadButton("downloadData", label="Download"),
      tableOutput("table")
      )
    ),
    # new tab for statistical analysis
    tabPanel(
      title="Statistical analysis",
      fluidRow(
        br(),
        column(
          2,
          checkboxInput(
            inputId = "log",
            label = HTML("<b>Convert to log2</b>")
          )
        ),
        column(
          3,
          pickerInput(
            inputId = "normalisation",
            label = "Choose normalisation type",
            choices = c("None", "Quantile", "Cloess")
          )
        ),
        column(
          3,
          pickerInput(
            inputId = "imputation",
            label = "Choose imputation type",
            choices = c("None", "MinDet", "ImpSeq", "MissForest", "MinProb"),
            choicesOpt = list(subtext = c("", "", "", "Very slow (maxiter=1)", ""))
          )
        )
      ),
      br(),
      h5("By default the normalisation and imputation is performed on all samples.\n
         If this should be done by condition, please put selected samples (at least two per group) into respective groups.\n"),
      fluidRow(
        column(
          3,
          selectInput(inputId = "group1", choices = sample_info$Sample_name, 
                      label = "Group_1", multiple = TRUE, selected = NULL)
        ),
        column(3, selectInput(inputId = "group2", choices = sample_info$Sample_name, 
                              label = "Group_2", multiple = TRUE, selected = NULL))
      ),
      br(),
      # output in the form of volcano and heatmap
      fluidRow(
        column(6, withSpinner(plotOutput("volcano"))),
        column(6, withSpinner(plotOutput("heatmap")))
      ),
      br(),
      # output in the form of table than can be downloaded and saved as csv file
      fluidRow(DTOutput("matrix_stat")),
      downloadButton("download", "Download the data"),
      fluidRow(column(7, dataTableOutput("dto")))
    )
  )
)
