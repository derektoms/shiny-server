#                          _        ____
#  _ __ ___  ___ ___ _ __ | |_ ___ |  _ \
# | '__/ _ \/ __/ _ \ '_ \| __/ _ \| |_) |
# | | |  __/ (_|  __/ |_) | || (_) |  _ <
# |_|  \___|\___\___| .__/ \__\___/|_| \_\
#                   |_|
#
# June 2019 receptoR v 1.3
## Last update: 2019-06-22, Derek Toms
## ui.R

########################################
#$#$#$#$#$#$    HEADER     $#$#$#$#$#$#$
########################################

library(shiny)
library(shinythemes)
library(shinyjs)


########################################
#$#$#$#$#$#$  Javascript   $#$#$#$#$#$#$
########################################


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

jscode2 <- '
$(function() {
    $(#receptorMain li a[data-value="Gene-level Expression"]).hide();
    $(#receptorMain li a[data-value="Sample-level Expression"]).hide();
    
});
'


########################################
#$#$#$#$#$#$#$    UI     $#$#$#$#$#$#$#$
########################################

ui <- fluidPage(
tags$head(tags$script(HTML(jscode))),
tags$head(tags$script(HTML(jscode2))),
tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "receptor.CSS")),
tags$head(tags$link(rel = "stylesheet", href = "https://use.fontawesome.com/releases/v5.6.3/css/all.css",  integrity="sha384-UHRtZLI+pbxtHCWp1t77Bi1L4ZtiqrqD80Kn4Z8NTSRyMA2Fd33n5dQ8lWUE00s/", crossorigin="anonymous"),
tags$head(tags$style(
        type="text/css",
        "#QC img {max-width: 100%; width: 100%; height: auto}"
    ))

),
# tags$script(HTML("$('body').addClass('fixed);")),
shinyjs::useShinyjs(),
navbarPage("receptoR",
    id = "receptorMain",
    theme = shinytheme("spacelab"),

# Start page  ------------------------------------------------------------------------------

    tabPanel("Start here",
       value ="startPanel",
       h3("Welcome to receptoR!"),
       hr(),
       sidebarLayout(
           sidebarPanel(
               # h4("An automated hypothesis generation software to identify cellular signaling pathways from transcriptomics data"),
               p("This software allows you to browse and analyze public transcriptomics data. This is based on the idea that each cell type expresses a particular suite of cellular receptors that drive its behaviour."),
               tags$ol(tags$li("A cell transcribes mRNA that will be translated into functional receptor proteins."),tags$li("Isolated RNA from this cell is converted to labeled cDNA, which is hybridized to an oligonucleotide probe array to measure expression."),tags$li("Each array represents a snapshot of a specific transcriptome. Thousands of these have been digitized and made publicly available."),tags$li("By mining this data, we can predict which receptors are expressed by our cells or tissues of interest to direct bioengineering strategies.")),
               hr(),
               
               #div
               p("There are two ways to begin using receptor, either by ",
               actionLink("linkSearch","searching for expression data"),
               " to design your own experiment, or by ",
               actionLink("linkLoad","loading and analysing an existing experiment.")),
               # To proceed, click \'Search for datasets\', above"),
               hr(),
               p("code created by Derek Toms, Qing Yun Tong and Matthew Workentine"),
               p("Copyright (C) 2019, code licensed under GPLv3")
               #/div
               ),
           mainPanel(
               img(src="overview.png",width="100%")
               ))
               
        ),

# Search for GSM  ------------------------------------------------------------------------------

    tabPanel("Search Microarray Database",
       value = "searchPanel",
       includeCSS("www/receptor.CSS"),
       h3("Organize publicly available expression data"),
       hr(),
       sidebarLayout(
           
       sidebarPanel(
           # style = "position:fixed;width:30%",
           conditionalPanel(condition="input.searchpanel==1",
           h4("Search Expression Data"),
           p("Begin by searching the Gene Expression Omnibus (GEO) database for publicly available transcriptome data. Depending on availability, these may be available for specific tissues or isolated cell types. In the next step, these samples will be assigned to one of three categories to determine differential expression between sample types."),
           br(),
           radioButtons("gplSelection", "Choose species:", choices = c("Mouse (GPL1261)" = "mouse", "Human (GPL570)" = "human")),
           tagAppendAttributes(textInput("searchText", "Enter search terms:", value = ""),`data-proxy-click` = "searchButton"),
           helpText("Search for multiple keywords using the boolean operators 'AND','OR','NOT', and the wildcard '*'. For example: 'liver AND hepa* NOT brain'."),
           actionButton("searchButton", "Search for arrays"),
           hr(),
           # HTML(paste("These experiments, each containing multiple biological samples, are refered to as ",span("G",style="font-weight:bold"),"EO data ",span("se",style="font-weight:bold"),"ries (GSE). Each ",span("G",style="font-weight:bold"),"EO ",span("s",style="font-weight:bold"), "a",span("m",style="font-weight:bold"),"ple (GSM) represents a digitized transcriptional snapshot.",sep="")),
           p("Click \'Add array to experiment\' to retrieve array (GSM) information and then click on the \'Assign\' tab above to organize this data for analysis."),
          
           actionButton("addButton", "Add array to experiment")),
           
           conditionalPanel(condition="input.searchpanel==2",
           h4("Define the categories that you wish to assign each sample (GSM) for comparison."),
           p("Each sample of interest should be assigned to a category. In this way, experimental comparisons can be performed to determine differential expression between categories. A minimum of two and a maximum of three categories should be defined. If you are only interested in a single sample type it is recommended that this is compared to a 'background' sample to identify enriched receptor genes."),

           tags$div(class="inputWithIcon", textInput("cat1", label=NULL, placeholder="Category 1 (e.g., pancreatic endocrine cells)"), tags$span(style="color:#E41A1C",icon("circle",class="fa-2x"))),

           tags$div(class="inputWithIcon",textInput("cat2", label=NULL, placeholder="Category 2 (e.g., photoreceptors)"), tags$span(style="color:#377EB8",icon("circle",class="fa-2x"))),
           
           tags$div(class="inputWithIcon",textInput("cat3", label=NULL, placeholder="Category 3 (optional)"),
           tags$span(style="color:#4DAF4A",icon("circle",class="fa-2x"))),
                      
           hr(),
           h4("Highlight samples, then click to Assign them to the specificed category."),
           p("Using the table at right and the drop down menu below, click on samples and \'Assign\' them to different categories. Samples can be filtered using the search bar."),
           fluidRow(column(8,uiOutput("categorySelect")),
           column(4,actionButton("assignButton", "Assign")))
           ),
           
           conditionalPanel(condition="input.searchpanel==3",
               h4("Thank you for using receptoR!"),
               p(" Please enter your name and any comments/bugs/questions/requests in the box below, then click the \'Download and Process\' button to retrieve the raw files from the NCBI server and process them based on their assigned categories."),
               textAreaInput("comments","Comments",width="100%",height="100px",resize="vertical"),
               textInput("downloadId","Download ID"),
               downloadButton("report","Download Report"),
               actionButton("downloadCEL","Process")),
           hr(),
               # Help banner on the bottom -------------------------
           # h4("Help me!"),
           p("Click ",actionLink("linkReset","here "),"to start again.")
       ),
       
       mainPanel(
           # Search GSE based on species
        tabsetPanel(
        tabPanel("Search", value=1,
            h4("GEO microarrays (\'GSM\') matching your search query"), # return search here!
            DT::dataTableOutput("searchResultsGSM")
        ),
        # Assign samples to categories ------------------------------------------------------
        tabPanel("Assign", value=2,
            h4("Assign individual arrays (GSM) to categories of your choosing"),
            DT::dataTableOutput("gsm_table")
        ),
        # This will be where the CEL files are downloaded (confirmation, etc) ------------
        tabPanel("Process", value=3,
        h4("Please confirm samples are properly categorized before proceeding"),
        p("Expression samples annotated:"),
                DT::dataTableOutput("finishedtable")
        ),
        
        id = "searchpanel"
        )
        
        )
        )
        
    ),
    
    # Load Gene Expression Data tab -------------------------------------
    tabPanel("Load Expression Datasets",
        value="expressionPanel",
        h3("Pick from user-defined experiments to perform analyses"),
        hr(),
        sidebarLayout(
        sidebarPanel(
            h4("Load Experiment"),
            uiOutput("loadUserExperiments"),
            hr(),
            checkboxGroupInput("genelist", "Select a receptor type to analyze", 
                  choices = NULL),
            br(),
            selectInput("gene", "Select additional non-receptor coding gene(s) to include in the analysis.", choices = NULL, multiple = TRUE),
            helpText("Search by gene symbol; availability of a given gene is based on microarray probe annotations."),
            downloadButton("reportDEG","Download differential gene expression analysis"),
            helpText("This Microsoft Excel file (.XSLX) contains all differentially expressed genes among these tissues. It can be further used for downstream analyses including functional enrichment analysis.")
        ),
        mainPanel(
            tabsetPanel(type="tabs",selected="Gene-by-gene Expression",
            tabPanel("Quality control",
            uiOutput("QC")
        ),
            # tabPanel("Experimental design",h4("Category definitions and contrasts"),p("Coming soon!")),
            tabPanel("Gene-by-gene Expression",
                fluidRow(
                column(6, h4("Average Expression"), DT::dataTableOutput("genes")),
                column(6, h4("Gene Violin Plot"), plotOutput("singleGenePlot"))
            )))
        )
        )
    ),
    
    # Magnitude expression tab ------------------------------------------------------------------------------
    
    tabPanel("Gene-level Expression",
        h3("Compare genes based on absolute expression and differential expression between experimental groups"),
        hr(),
        sidebarLayout(
        sidebarPanel(
            h4("Gene expression"),
            checkboxGroupInput("tissues", label = "Select tissues to include",
            choices = NULL, selected = NULL),
            br(),
            checkboxInput("de_state", label = "Show differential expressed only", value = TRUE),
            checkboxGroupInput("de", label = "Choose comparison(s) to show", choices = NULL, selected = NULL),
            br(),
            conditionalPanel(condition="input.absexpanel==1",
                h5("Heatmap parameters"),
                checkboxInput("hm_probes", "Show probe-level", value = FALSE),
                checkboxInput("hm_gsm", "Show GSM (column names)", value = TRUE),
                checkboxInput("hm_rownames", "Show gene symbols (row names)", value = TRUE),
                checkboxInput("hm_col_cluster", "Cluster columns", value = TRUE),
                checkboxInput("hm_row_cluster", "Cluster rows", value = TRUE),
                numericInput("hm_width", "Plot width (px)", value = 900, min = 100, max = 2400, step = 10),
                numericInput("hm_height", "Plot height (px)", value = 1200, min = 100, max = 2400, step = 10))
        ),
        mainPanel(
            tabsetPanel(type = "tabs",
            tabPanel("Heatmap", value=1, h4("Cluster analysis and a heatmap representation of gene expression."), p("Genes with similar expression patterns will cluster as rows, while individual microarrays cluster as columns."), uiOutput("heatmap_ui")),
            tabPanel("Summary boxplots", h4("Boxplots of expression data by tissue."), plotOutput("overallPlot", height = 600)),
            tabPanel("By-gene boxplots", h4("Boxplots of expression data by gene."), plotOutput("byGenePlot", height = 600)),
            id = "absexpanel"
        )
        )
        )
    ),

    # Mixomics tab ---------------------------------------------
    tabPanel("Sample-level Expression",
        h3("Compare trends in samples based on experimental groups"),
        hr(),
        sidebarLayout(
        sidebarPanel(
            h4("Sample expression"),
            checkboxGroupInput("pls_tissues", label = "Select tissues to inclued",
            choices = NULL, selected = NULL),
            # checkboxInput("pls_probe", "Perform PLS-DA at probe level", value = FALSE),
            br(),
            h4("Gene contribution plot"),
            uiOutput("numGenesUI"),
            radioButtons("pls_ncomp", "Select component for gene contribution plot", choices = c(1,2)),
            br()
            # downloadButton("pls_download", "Download gene contribution data")
        ),
        mainPanel(
            tabsetPanel(type = "tabs",
            tabPanel("Discriminant Analysis", h4("Sparse Partial Least Squares Discriminant Analysis (sPLS-DA)  selects genes that are informative about a specific group. "), plotOutput("indPlot", height = 800)),
            
            tabPanel("Component loadings plot", h4("Gene contribution to each principle component."), p("The longer the bar (in either direction) the more that gene contributes to that component."), plotOutput("contribPlot", height = 800)),
            
            tabPanel("Circle variance", h4("Circle variance projections onto tissue."), p("Strongly correlated genes are projected in the same direction from the origin; the greater the distance the stronger the association."), plotOutput("varPlot", height = 800))
        ),
        position = c("right","left"),
        fluid = TRUE
        )
        )
    )
)
)