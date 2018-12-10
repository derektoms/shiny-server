## 2018-04-02
library(shiny)
library(shinythemes)
library(shinyjs)

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
    theme = shinytheme("spacelab"),
    
    tabPanel("Search for datasets",
        # Search for GSE  ------------------------------------------------------------------------------              
        tabsetPanel(
        tabPanel("Search for GEO data series (GSE)",
            h4("Highlight the desired search results (GSE) and click 'Retrieve GSM' to proceed"),
            # fluidRow(
                # column(6,
                radioButtons("gplSelection", "Choose species:", choices = c("Mouse (GPL1261)" = "mouse", "Human (GPL570)" = "human")),
                tagAppendAttributes(textInput("Key", "Enter search terms, separated by commas", value = ""),`data-proxy-click` = "Search"),
                # ),
                # column(6,
                actionButton("Search", "Search"),
                actionButton("getGSM", "Retrieve GSM"),
            #     )
            # ),
            hr(),
            DT::dataTableOutput("filteredgse")
        ),
        # Assign samples to categories ------------------------------------------------------
        tabPanel("Assign samples (GSM) to categories",
            h4("Define the categories that you wish to assign each sample (GSM) for comparison."),
            # uiOutput("categorySelect"
            # fluidRow(
            #     column(5,
                    textInput("cat1", label=NULL, placeholder="Category 1"),
                    textInput("cat2", label=NULL, placeholder="Category 2"),
                    textInput("cat3", label=NULL, placeholder="Category 3"),
                #     ),
                # column(7,
                    # actionButton("Remove", "Finalize selections and remove not included"),#))),
            
            # ^ why aren't these two the same button?
            # helpText("Do not click 'finish' until all selections have been made. This button removes the unselected rows and generates a new table on the next page."),
            
            hr(),
            DT::dataTableOutput("gsm_table"),
            # Category assignment panel
            absolutePanel(id="assignCat",class="panel panel-default",fixed=TRUE,draggable=TRUE,top=120,left="auto",right=20,bottom="auto",width="420",height=50,
            wellPanel(
                h4("Highlight samples, then click to Assign them to the specificed category."),
                uiOutput("categorySelect"),
                actionButton("Assign", "Assign GSM to Category")
                )
            )
        ),
            # This will be where the CEL files are downloaded (confirmation, etc) ------------
        tabPanel("Confirm sample categories", uiOutput("page4"),
            h4("Please check that your samples are appropriately categorized."),
            fluidRow(
                column(8,DT::dataTableOutput("finishedtable")),
                column(4,actionButton("downloadCEL","Download CEL files")))
        )
        )
    ),
    
    # Load Gene Expression Data tab -------------------------------------
    tabPanel("Load expression data",
        sidebarLayout(
        sidebarPanel(
            selectInput(inputId="user_data",label="Select user data for analysis",choices=c("none"="none","Photoreceptors v RPE"="2018-04-13_app_data.rda"),selected="none"),
            br(),
            uiOutput("geneListsUI"),
            br(),
            uiOutput("geneUI")
        ),
        mainPanel(
            tabsetPanel(type="tabs",
            tabPanel("Experimental design",h4("Category definitions and contrasts")),
            tabPanel("Gene-level expression",
                fluidRow(
                column(6, h4("Average Expression"), DT::dataTableOutput("genes")),
                column(6, h4("Gene Boxplot"), plotOutput("singleGenePlot"))
            )))
        )
        )
    ),
    
    # Magnitude expression tab ------------------------------------------------------------------------------

    tabPanel("Absolute expression", 
        sidebarLayout(
        sidebarPanel(
            style = "position:fixed",
            checkboxGroupInput("tissues", label = "Select tissues to inclued",
            choices = c("photoreceptors","RPE","whole.retina"), selected = c("photoreceptors","RPE","whole.retina")
            ),
            br(),
            checkboxInput("de_state", label = "Show differential expressed only", value = TRUE),
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
    tabPanel("Relative expression",
        sidebarLayout(
        sidebarPanel(
            checkboxGroupInput("pls_tissues", label = "Select tissues to inclued",
            choices = c("photoreceptors","RPE","whole.retina"), selected = c("photoreceptors","RPE","whole.retina")
            ),
            checkboxInput("pls_probe", "Perform PLS-DA at probe level", value = FALSE),
            br(),
            h4("Gene contribution plot"),
            uiOutput("numGenesUI"),
            radioButtons("pls_ncomp", "Select component for gene contribution plot", choices = c(1,2)),
            br()
            # downloadButton("pls_download", "Download gene contribution data")
        ),
        mainPanel(
            tabsetPanel(type = "tabs",
            tabPanel("PCA Analysis", plotOutput("indPlot", height = 800)),
            tabPanel("Circle variance", plotOutput("varPlot", height = 800)),
            tabPanel("Loadings plot", plotOutput("contribPlot", height = 800))
        )
        )
        )
    )
)
)