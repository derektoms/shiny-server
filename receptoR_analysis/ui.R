## receptoRA ui.R

ui<-fluidPage(
  
  #useShinyjs(),
  navbarPage("receptoR expression Analysis",
    
    theme = "sandstone.css",
    
  # Load Genes ----------------------------------------------------------------------------------  
    tabPanel("Load Genes",
      sidebarLayout(
        sidebarPanel(
          uiOutput("geneListsUI"),
          br(),
          uiOutput("geneUI")
        ),
        mainPanel(
          fluidRow(
            column(6,
              h4("Average Expression"),
              DT::dataTableOutput("genes")),
            column(6,
              h4("Gene Boxplot"),
              plotOutput("singleGenePlot"))
          )
        )
      )
    ),
    
  
  # Expression tab ------------------------------------------------------------------------------
  
    tabPanel("Expression plots", 
      sidebarLayout(
        sidebarPanel(
          checkboxGroupInput("tissues", label = "Select tissues to inclued",
              choices = c("photoreceptors","RPE","whole.retina"), selected = c("photoreceptors","RPE","whole.retina")),
          br(),
          checkboxInput("de_state", label = "Show differential expressed only", value = FALSE),
          uiOutput("de_choices"),
          br(),
          h3("Heatmap parameters"),
          checkboxInput("hm_probes", "Show probe-level", value = FALSE),
          checkboxInput("hm_rownames", "Show rownames", value = TRUE),
          checkboxInput("hm_col_cluster", "Cluster columns", value = TRUE),
          checkboxInput("hm_row_cluster", "Cluster rows", value = TRUE),
          h4("Select plot dimensions (px)"),
          numericInput("hm_width", "Width", value = 900, min = 100, max = 2400, step = 10),
          numericInput("hm_height", "Height", value = 1200, min = 100, max = 2400, step = 10)
        ),
        mainPanel(
          tabsetPanel(type = "tabs",
            tabPanel("Heatmap", uiOutput("heatmap_ui")),
            tabPanel("Summary boxplots", plotOutput("overallPlot", height = 600)),
            tabPanel("By-gene boxplots", plotOutput("byGenePlot", height = 600))
          )
        )
      )
    )

    # tabPanel("PLS-DA",
    #   sidebarLayout(
    #     sidebarPanel(
    #       checkboxGroupInput("pls_tissues", label = "Select tissues to inclued",
    #           choices = groups, selected = groups),
    #       checkboxInput("pls_probe", "Perform PLS-DA at probe level", value = FALSE),
    #       br(),
    #       h4("Gene contribution plot"),
    #       uiOutput("numGenesUI"),
    #       radioButtons("pls_ncomp", "Select component for gene contribution plot", choices = c(1,2)),
    #       br()
    #       # downloadButton("pls_download", "Download gene contribution data")
    #     ),
    #     mainPanel(
    #       plotOutput("indPlot", height = 800),
    #       plotOutput("varPlot", height = 800),
    #       plotOutput("contribPlot", height = 800)
    #       # DT::dataTableOutput("contribTable")
    #      )
    #    )
    #  )
  )
)