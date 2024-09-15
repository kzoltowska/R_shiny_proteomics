#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Author: Katarzyna Zoltowska
# Aim: Testing proteomics data analysis workflows

ui <- fluidPage(
  position="fixed-top",
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
            choices = c(
              "Entire dataset", "Tabular summary", "Number of identifications & NA summary", 
              "PCA and correlation plots"
            
            )
          ),
          checkboxInput(
            inputId = "Keep_one",
            label = "Keep only one uniprot/geneID. Shortest uniprot ID is kept.",
            value = TRUE
          ),
          textInput(inputId = "filter", label = "Write regex matching input of gene column \n
                  to filter out selected genes", value = "")
        ),
        mainPanel(
          width = 10,
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
            condition = "input.Description == 'Number of identifications & NA summary'",
            plotOutput("barplot"),
            br(),
            br(),
            plotOutput("summary_plot")
          )
          ,
          # Show a graphical summary
          conditionalPanel(
            condition = "input.Description == 'PCA and correlation plots'",
            fluidRow(
            (column(5,plotOutput("pca"))),
            (column(5, plotOutput("cor_plot")))
          )
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
            choices = c("None", "Quantile", "Cloess", "Median")
          )
        ),
        column(
          3,
          pickerInput(
            inputId = "imputation",
            label = "Choose imputation type",
            choices = c("None", "MinDet", "ImpSeq", "MissForest", "MinProb", "KNN")
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
      fluidRow(
        column(3, withSpinner(plotOutput("dens_raw"))),
        column(3, withSpinner(plotOutput("dens_log"))),
        column(3, withSpinner(plotOutput("dens_norm"))),
        column(3, withSpinner(plotOutput("dens_imp")))
      ),
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
          choices = c("ROTS", "limma", "ttest-BH")
        )),
        column(2, actionButton(inputId="calc_stat","Calculate"))
      ),
      fluidRow(
        column(
          2,
          radioButtons(
            inputId = "pval",
            label = "P value to show on volcano plot",
            choices = c("p value", "adjusted p value / FDR")
          )),
        column(2, sliderInput(
          inputId = "slider_p",
          label = "P value / FDR cut-off for volcano plot", min = 0,
          max = 1, value = 0.05
        )),
        column(2, sliderInput(
          inputId = "slider_fc",
          label = "LogFC cut-off for volcano plot", min = 0, max = 4, value = 1, step=0.01
        )),
        column(3, numericInput(
          inputId = "pval_heat_tab", label = "FDR cut-off for heatmap and table",
          value = 0.05, min = 0, max = 1, step = NA
        )),
        column(3, numericInput(
          inputId = "logfc_heat_tab_min", label = "Lower logFC cut-off for heatmap and table",
          value = (-1), min = NA, max = NA, step = NA
        )),
        column(3, numericInput(
          inputId = "logfc_heat_tab_max", label = "Upper logFC cut-off for heatmap and table",
          value = 1, min = NA, max = NA, step = NA
        )),
        column(3, checkboxInput(inputId = "Rownames_heat", label = "Display rownames on heatmap", value = FALSE))
      ),
      br(),
      br(),
      # output in the form of volcano and heatmap
      fluidRow(
        column(6, withSpinner(plotOutput("volcano"))),
        column(6, withSpinner(plotOutput("heatmap_stat")))
      ),
      br(),
      br(),
      fluidRow(DTOutput("summary_num")),
      br(),
      br(),
      # output in the form of table than can be downloaded and saved as csv file
      fluidRow(
        column(2,DTOutput("stat_table"))),
      fluidRow(
        downloadButton("downloadData_stat", label = "Download"),
        tableOutput("table_stat")
      )
    ),
    # new tab for functional analysis
    tabPanel(
      title = "ORA",
      br(),
      br(),
      fluidRow(
        column(2, selectInput(inputId="GOset", label="Choose gene ontology", 
                              choices=c("BP","CC", "MF"))),
        column(2, radioButtons(inputId="organism", label="Species", choices=c("human", "mouse"))),
        column(2, actionButton(inputId="calc_ora",label="Calculate")),
        column(3,numericInput(label="Number of categories on dot plot", inputId = "ndot", value=10  )),
        column(3,numericInput(label="Number of categories on cnet plot", inputId = "ncnet", value=10 ))
      ),
      br(),
      fluidRow(column(6, withSpinner(plotOutput("dotplot"))),
               column(6, withSpinner(plotOutput("cnetplot")))),
    
    br(),
    br(),
    fluidRow(withSpinner(DTOutput("oratable"))),
    br(),
    
    fluidRow(column(6, withSpinner(plotOutput("simheat"))),
             column(6, withSpinner(plotOutput("simsquare")))),
    fluidRow(withSpinner(DTOutput("reducedTerms")))
  ),
  tabPanel(
    title = "DOSE",
    br(),
    br(),
    fluidRow(
      column(3,numericInput(label="Number of categories on dot plot", inputId = "ndot_DO", value=10  )),
      column(3,numericInput(label="Number of categories on cnet plot", inputId = "ncnet_DO", value=10 )),
      column(2, actionButton(inputId="calc_dose",label="Calculate"))
    ),
    br(),
    fluidRow(column(6, withSpinner(plotOutput("dotplot_DO"))),
             column(6, withSpinner(plotOutput("cnetplot_DO")))),
    
    br(),
    br(),
    fluidRow(withSpinner(DTOutput("DOtable"))),
    br(),
  )
)
)
