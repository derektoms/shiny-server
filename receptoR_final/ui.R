## 2018-04-02

## Javascript

# this codes for the 'enter/return' key as an action button
jscode <- '
$(function() {
  var $els = $("[data-proxy-click]");
  $.each(
    $els,
    function(idx, el) {
      var $el = $(el);
      var $proxy = $("#" + $el.data("proxyClick"));
      $el.keydown(function (e) {
        if (e.keyCode == 13) {
          $proxy.click();
        }
      });
    }
  );
});
'


## UI
ui <- fluidPage(
tags$head(tags$script(HTML(jscode))),
# tags$script(HTML("$('body').addClass('fixed);")),
shinyjs::useShinyjs(),
navbarPage("receptoR",
    theme = "sandstone.css",

    tabPanel("Search for datasets",
        # Search for GSE  ------------------------------------------------------------------------------              
        sidebarLayout(
        sidebarPanel(
        # Search for datasets ------------------------------------------------------
        h4("Search for GEO data series (GSE)"),
        radioButtons("gplSelection", "Choose species:", choices = c("Mouse (GPL1261)" = "mouse", "Human (GPL570)" = "human")),
        tagAppendAttributes(
        textInput("Key", "Enter search terms, separated by commas", value = ""),
        `data-proxy-click` = "Search"
        ),
        actionButton("Search", "Search"),
        hr(),
        # Define categories --------------------------------------------------------
        h4("Define the categories that you wish to assign each sample (GSM) for comparison."),
        textInput("cat1", "Define Category 1"),
        textInput("cat2", "Define Category 2"),
        textInput("cat3", "Define Category 3")
        ),
        mainPanel(
            tabsetPanel(
            tabPanel("Select GEO data series (GSE)",  #search GSE, and select which to include
                h4("2. Highlight the desired search results (GSE) and click 'Retrieve GSM' to proceed"),
                actionButton("getGSM", "Retrieve GSM"),
                helpText("Do not click 'finish' until all selections have been made. This button removes the unselected rows and generates a new table on the next page."),
                DT::dataTableOutput("filteredgse")
            ),
            # Assign samples to categories ------------------------------------------------------
            tabPanel("Assign samples (GSM) to categories", 
                h4("4. Highlight the desired search results and click 'assign' to assign them to the specificed category"),
                uiOutput("categorySelect"),
                actionButton("Assign", "Assign Categories"),
                actionButton("Remove", "Finalize selections and remove not included"),
                # ^ why aren't these two the same button?
                helpText("Do not click 'finish' until all selections have been made. This button removes the unselected rows and generates a new table on the next page."),
                DT::dataTableOutput("gsm_table")
            ),
            # This will be where the CEL files are downloaded (confirmation, etc) ------------
            tabPanel("Selection details", uiOutput("page4"), 
                DT::dataTableOutput("finishedtable"),
                actionButton("downloadCEL","Download CEL files"),
                tableOutput("CELdl")
            )
            )
        )
        )
    ),

    # Load Genes tab -------------------------------------
    tabPanel("Load data",
        sidebarLayout(
        sidebarPanel(
            selectInput(inputId="user_data",label="Select user data for analysis",choices=c("none"="none","Photoreceptors v RPE"="2018-04-13_app_data.rda"),selected="2018-04-13_app_data.rda"),
            br(),
            uiOutput("geneListsUI"),
            br(),
            uiOutput("geneUI")
        ),
        mainPanel(
            fluidRow(
                column(6, h4("Average Expression"), DT::dataTableOutput("genes")),
                column(6, h4("Gene Boxplot"), plotOutput("singleGenePlot"))
            )
        )
        )
    ),
    # Expression tab ------------------------------------------------------------------------------

    tabPanel("Expression plots", 
        sidebarLayout(
        sidebarPanel(
            checkboxGroupInput("tissues", label = "Select tissues to inclued",
            choices = c("photoreceptors","RPE","whole.retina"), selected = c("photoreceptors","RPE","whole.retina")
            ),
            br(),
            checkboxInput("de_state", label = "Show differential expressed only", value = FALSE),
            uiOutput("de_choices"),
            br(),
            h3("Heatmap parameters"),
            checkboxInput("hm_probes", "Show probe-level", value = FALSE),
            checkboxInput("hm_gsm", "Show GSM (column names)", value = TRUE),
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
    ),

    # Mixomics tab ---------------------------------------------
    tabPanel("PLS-DA",
        sidebarLayout(
        sidebarPanel(
            checkboxGroupInput("pls_tissues", label = "Select tissues to inclued",
            choices = groups, selected = groups),
            checkboxInput("pls_probe", "Perform PLS-DA at probe level", value = FALSE),
            br(),
            h4("Gene contribution plot"),
            uiOutput("numGenesUI"),
            radioButtons("pls_ncomp", "Select component for gene contribution plot", choices = c(1,2)),
            br()
            # downloadButton("pls_download", "Download gene contribution data")
        ),
        mainPanel(
            plotOutput("indPlot", height = 800),
            plotOutput("varPlot", height = 800),
            plotOutput("contribPlot", height = 800)
            # DT::dataTableOutput("contribTable")
        )
        )
    )
)
)