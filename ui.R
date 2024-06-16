#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Author: Katarzyna Zoltowska
# Aim: Testing proteomics data analysis workflows

ui <- fluidPage(
  # Application title
  titlePanel("Proteomics data analysis"),
  h4("Hobby project"),
  # create tabset
  tabsetPanel(
    # create first tab
    tabPanel(
      title = "Data overview", fluid = TRUE,
      sidebarLayout(
        sidebarPanel(
          width = 2,
          fileInput("file1", "Choose data file"),
          fileInput("file2", "Choose sample info file"),
          # create radiobuttons to select how to view the data
          radioButtons(
            inputId = "Description",
            label = "How would you like to view the data?",
            choices = c("Entire dataset", "Tabular summary", "Graphical summary", "Boxplot - raw intensities", "Boxplot - log2 intensities")
          ),
          checkboxInput(inputId="Keep_one",
                        label="Keep only one uniprot/geneID. Shortest uniprot ID is kept.",
                        value=TRUE),
        textInput(inputId="filter", label="Write regex matching input of gene column \n
                  to filter out selected genes", value="")
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
          selectInput(
            inputId = "group1", choices = NA,
            label = "Group_1_normalisation_imputation",
            multiple = TRUE, selected = NULL
          )
        ),
        column(3, selectInput(
          inputId = "group2", choices = NA,
          label = "Group_2_normalisation_imputation",
          multiple = TRUE, selected = NULL
        ))
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
        downloadButton("downloadData", label = "Download"),
        tableOutput("table")
      )
    ),
    # new tab for statistical analysis
    tabPanel(
      title = "Statistical analysis",
      h5("Please select groups for statistical comparison.\n"),
      fluidRow(
        column(
          2,
          selectInput(
            inputId = "group1_stat", choices = NA,
            label = "Group_1_statistics", multiple = TRUE, selected = NULL
          )
        ),
        column(2, selectInput(
          inputId = "group2_stat", choices = NA,
          label = "Group_2_statistics", multiple = TRUE, selected = NULL
        )),
        column(2, selectInput(
          inputId = "stat", 
          label = "Statistical test to use:",
          choices = c("ROTS", "limma"))),
        column(2, 
               radioButtons(inputId="pval", 
                            label="P value to show on volcano plot", 
                            choices=c("p value", "adjusted p value / FDR")))),
      # output in the form of volcano and heatmap
      fluidRow(
        withSpinner(plotOutput("volcano"))
      )
    )
  )
)
