#                          _        ____
#  _ __ ___  ___ ___ _ __ | |_ ___ |  _ \
# | '__/ _ \/ __/ _ \ '_ \| __/ _ \| |_) |
# | | |  __/ (_|  __/ |_) | || (_) |  _ <
# |_|  \___|\___\___| .__/ \__\___/|_| \_\
#                   |_|
#
# March 2019 receptoR v 1.2
## Last update: 2019-04-14, Derek Toms
## ui.R

########################################
#$#$#$#$#$#$    HEADER     $#$#$#$#$#$#$
########################################

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


########################################
#$#$#$#$#$#$#$    UI     $#$#$#$#$#$#$#$
########################################

ui <- fluidPage(
tags$head(tags$script(HTML(jscode))),
tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "receptor.css")),
tags$head(tags$link(rel = "stylesheet", href = "https://use.fontawesome.com/releases/v5.6.3/css/all.css",  integrity="sha384-UHRtZLI+pbxtHCWp1t77Bi1L4ZtiqrqD80Kn4Z8NTSRyMA2Fd33n5dQ8lWUE00s/", crossorigin="anonymous")
),
# tags$script(HTML("$('body').addClass('fixed);")),
shinyjs::useShinyjs(),
navbarPage("receptoR",
    theme = shinytheme("spacelab"),

# Start page  ------------------------------------------------------------------------------

    tabPanel("Start here",
       h3("Welcome to receptoR!"),
       hr(),
       sidebarLayout(
           sidebarPanel(
               # h4("An automated hypothesis generation software to identify cellular signaling pathways from transcriptomics data"),
               p("This software allows you to browse and analyze public transcriptomics data. This is based on the idea that each cell type expresses a particular suite of cellular receptors that drive its behaviour."),
               tags$ol(tags$li("A cell transcribes mRNA that will be translated into functional receptor proteins."),tags$li("Isolating RNA from the cell and converting it to labeled cDNA allows us to hybridize it to an probe array to measure expression."),tags$li("Each sample represents a particular transcriptomic snapshot. Thousands of these have been digitized and made publicly available."),tags$li("By mining this data, we can predict which receptors are expressed by our samples of interest to direct tissue engineering strategies.")),
               hr(),
               
               #div
               p("There are two ways to begin using receptor, either by searching for expression data to design your own experiment, or by loading and analysing an existing experiment."),
               # To proceed, click \'Search for datasets\', above"),
               hr(),
               p("(C) 2019 Derek Toms"),
               p("License")
               #/div
               ),
           mainPanel(
               img(src="overview.png",width="100%")
               ))
               
        ),

# Search for GSM  ------------------------------------------------------------------------------

    tabPanel("Search Expression Data",
       h3("Organize publicly available expression data"),
       hr(),
       sidebarLayout(
           
       sidebarPanel(
           # style = "position:fixed;width:30%",
           conditionalPanel(condition="input.searchpanel==1",
           h4("Search Expression Data"),
           p("Begin by searching for experiments that expression data for your cell or tissue type of interest."),
           br(),
           radioButtons("gplSelection", "Choose species:", choices = c("Mouse (GPL1261)" = "mouse", "Human (GPL570)" = "human")),
           tagAppendAttributes(textInput("searchText", "Enter search terms:", value = ""),`data-proxy-click` = "searchButton"),
           actionButton("searchButton", "Search for arrays"),
           hr(),
           # HTML(paste("These experiments, each containing multiple biological samples, are refered to as ",span("G",style="font-weight:bold"),"EO data ",span("se",style="font-weight:bold"),"ries (GSE). Each ",span("G",style="font-weight:bold"),"EO ",span("s",style="font-weight:bold"), "a",span("m",style="font-weight:bold"),"ple (GSM) represents a digitized transcriptional snapshot.",sep="")),
           p("Click \'Add array to experiment\' to retrieve array (GSM) information and then click on the \'Assign\' tab above to organize this data for analysis."),
          
           actionButton("addButton", "Add array to experiment")),
           
           conditionalPanel(condition="input.searchpanel==2",
           h4("Define the categories that you wish to assign each sample (GSM) for comparison."),
           p("Each sample of interest should be assigned to a category. In this way, experimental comparisons can be performed to determine differential expression between categories."),

           tags$div(class="inputWithIcon",textInput("cat1", label=NULL, placeholder="Category 1")
           # tags$span(style="color:#E41A1C",icon("circle",class="fa-2x"))
           ),
           
           tags$div(class="inputWithIcon",textInput("cat2", label=NULL, placeholder="Category 2")
           # tags$span(style="color:#377EB8",icon("circle",class="fa-2x"))
           ),
           
           tags$div(class="inputWithIcon",textInput("cat3", label=NULL, placeholder="Category 3 (optional)")
           # tags$span(style="color:#4DAF4A",icon("circle",class="fa-2x"))
           ),

           ### https://www.aridhia.com/blog/the-sky-is-not-the-limit-embedding-raw-html-and-javascript-to-create-dynamic-ui-elements-in-shiny-applications/   
           ### ^ this should help with dynamically adding/subtracting categories
           
           hr(),
           h4("Highlight samples, then click to Assign them to the specificed category."),
           p("Using the table at right and the drop down menu below, click on samples and \'Assign\' them to different categories. Samples can be filtered using the search bar. \nPLEASE NOTE: once you have clicked the \'Assign\' button you will no longer be able to add arrays to your experiment."),
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
               hr()
               # Help banner on the bottom -------------------------
               # h4("Help me!"),
               # p("Turducken leberkas t-bone tongue, tail frankfurter corned beef strip steak buffalo picanha beef tri-tippork belly rump flank. Chicken cupim sausage, spare ribs prosciutto beef pork corned beef salami leberkas shankle.",style="color:#D8BFD8")
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
    tabPanel("Load Experiment",
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
            selectInput("gene", "Select gene(s) to show", choices = NULL, multiple = TRUE),
            downloadButton("reportDEG","Download differential gene expression analysis")
        ),
        mainPanel(
            tabsetPanel(type="tabs",selected="Gene-level expression",
            tabPanel("Quality control",
            uiOutput("QC"),
            plotOutput("degPlot")
        ),
            tabPanel("Experimental design",h4("Category definitions and contrasts"),p("Coming soon!")),
            tabPanel("Gene-level expression",
                fluidRow(
                column(6, h4("Average Expression"), DT::dataTableOutput("genes")),
                column(6, h4("Gene Boxplot"), plotOutput("singleGenePlot"))
            )))
        )
        )
    ),
    
    # Magnitude expression tab ------------------------------------------------------------------------------
    
    tabPanel("Absolute Expression",
        h3("Compare genes based on absolute expression"),
        hr(),
        sidebarLayout(
        sidebarPanel(
            h4("Absolute expression"),
            # p("Bacon ipsum dolor amet chuck tongue flank bresaola corned beef hamburger leberkas pig bacon pork loin. Andouille hamburger strip steak ground round, ham filet mignon swine kielbasa pork chop jerky.",style="color:#D8BFD8"),
            # style = "position:fixed",
            checkboxGroupInput("tissues", label = "Select tissues to inclued",
            choices = NULL, selected = NULL),
            br(),
            checkboxInput("de_state", label = "Show differential expressed only", value = TRUE),
            checkboxGroupInput("de", label = "Choose comparison(s) to show", choices = NULL, selected = NULL),
            br(),
            conditionalPanel(condition="input.absexpanel==1",
                h5("Heatmap parameters"),
                checkboxInput("hm_probes", "Show probe-level", value = FALSE),
                checkboxInput("hm_gsm", "Show GSM (column names)", value = TRUE),
                checkboxInput("hm_rownames", "Show rownames", value = TRUE),
                checkboxInput("hm_col_cluster", "Cluster columns", value = TRUE),
                checkboxInput("hm_row_cluster", "Cluster rows", value = TRUE),
                numericInput("hm_width", "Plot width (px)", value = 900, min = 100, max = 2400, step = 10),
                numericInput("hm_height", "Plot height (px)", value = 1200, min = 100, max = 2400, step = 10))
        ),
        mainPanel(
            tabsetPanel(type = "tabs",
            tabPanel("Heatmap", value=1, uiOutput("heatmap_ui")),
            tabPanel("Summary boxplots", plotOutput("overallPlot", height = 600)),
            tabPanel("By-gene boxplots", plotOutput("byGenePlot", height = 600)),
            id = "absexpanel"
        )
        )
        )
    ),

    # Mixomics tab ---------------------------------------------
    tabPanel("Relative Expression",
        h3("Compare genes based on relative expression between experimental groups"),
        hr(),
        sidebarLayout(
        sidebarPanel(
            h4("Relative expression"),
            # p("Bacon ipsum dolor amet chuck tongue flank bresaola corned beef hamburger leberkas pig bacon pork loin. Turducken leberkas t-bone tongue, tail frankfurter corned beef strip steak buffalo picanha beef tri-tip pork belly rump flank. Chicken cupim sausage, spare ribs prosciutto beef pork corned beef salami leberkas shankle. Andouille hamburger strip steak ground round, ham filet mignon swine kielbasa pork chop jerky.",style="color:#D8BFD8"),
            checkboxGroupInput("pls_tissues", label = "Select tissues to inclued",
            choices = NULL, selected = NULL),
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
        ),
        position = c("right","left"),
        fluid = TRUE
        )
        )
    )
)
)