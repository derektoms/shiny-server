## 2019-01-11
## Beta version to be completed today!


## 2018-04-02
library(shiny)
library(shinythemes)
library(shinyjs)
#  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  
# ( )( )( )( )( )( )( )( )( )( )( )( )( )( )( )( )( )( )( )( )( )
# \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/\ 
# (_)(_)(_)(_)(_)(_)(_)(_)(_)(_)(_)(_)(_)(_)(_)(_)(_)(_)(_)(_)(_)

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

# Start page  ------------------------------------------------------------------------------

    tabPanel("Start here",
       h3("Welcome to receptoR!"),
       hr(),
       sidebarLayout(
           sidebarPanel(
               h4("receptoR is an automated hypothesis generation software to identify cellular signaling pathways from transcriptomics data."),
               p("Bacon ipsum dolor amet chuck tongue flank bresaola corned beef hamburger leberkas pig bacon pork loin. Turducken leberkas t-bone tongue, tail frankfurter corned beef strip steak buffalo picanha beef tri-tip pork belly rump flank. Chicken cupim sausage, spare ribs prosciutto beef pork corned beef salami leberkas shankle. Andouille hamburger strip steak ground round, ham filet mignon swine kielbasa pork chop jerky. Andouille t-bone biltong bacon beef ribs boudin frankfurter ham hock pork loin capicola tail ground round brisket tenderloin tri-tip. Ham pork bacon strip steak, ball tip leberkas meatball capicola pork loin. Ribeye turducken tri-tip jowl filet mignon drumstick shank corned beef prosciutto spare ribs sausage leberkas cupim burgdoggen bacon."),
               p("(C) 2019 Derek Toms")),
           mainPanel(
               p("This software allows you to browse and analyze public transcriptomics data. To proceed, click \'Search for datasets\', above"),
               img(src="overview.png"),
               p("Spicy jalapeno bacon ipsum dolor amet brisket ribeye tri-tip tail meatloaf ground round salami fatback. Ball tip flank pork turkey prosciutto. Venison burgdoggen beef pork hamburger tongue shankle rump frankfurter kielbasa pastrami pork belly. T-bone fatback venison tenderloin salami biltong turkey chuck."),

p("Turkey tenderloin buffalo frankfurter, strip steak capicola filet mignon ribeye cow t-bone biltong. Tri-tip hamburger ham hock, shank landjaeger kielbasa pork loin spare ribs pastrami salami shoulder alcatra. Biltong tail brisket, turkey beef ribs hamburger prosciutto bacon tenderloin pancetta shank venison picanha cupim. Frankfurter porchetta turkey biltong corned beef, burgdoggen prosciutto jerky tail. Capicola beef ribs burgdoggen pancetta, frankfurter leberkas pig. Filet mignon jerky ground round cow burgdoggen, shoulder boudin ham hock pancetta."),

p("Sirloin salami strip steak burgdoggen pork chop ribeye, pastrami fatback short ribs corned beef. Beef shank shoulder, leberkas pancetta ground round jowl tenderloin ribeye hamburger shankle flank. Strip steak short loin kevin pork belly meatball swine sausage ham jowl pig tongue venison fatback t-bone beef. Spare ribs shoulder ham hock strip steak pastrami.")
               ))
               
        ),

# Search for GSE  ------------------------------------------------------------------------------

    tabPanel("Search for datasets",
       h3("Search for and categorize expression datasets"),
       hr(),
       sidebarLayout(
       sidebarPanel(
           conditionalPanel(condition="input.searchpanel==1",
           h4("Find expression data from publicly available datasets to begin defining your own experiment."),
           p("Begin by searching for ",span("G",style="font-weight:bold"),"EO data ",span("se",style="font-weight:bold"),"ries (GSE) containing expression data for your cell or tissue type of interest."),
           radioButtons("gplSelection", "Choose species:", choices = c("Mouse (GPL1261)" = "mouse", "Human (GPL570)" = "human")),
           tagAppendAttributes(textInput("Key", "Enter search terms, separated by commas", value = ""),`data-proxy-click` = "Search"),
           actionButton("Search", "Search"),
           hr(),
           p("When you have selected the appropriate datasets, click \'Retrieve GSM\' to collect sample information and then click on the \'Assign samples\' tab above."),
           actionButton("getGSM", "Retrieve GSM")),
           
           conditionalPanel(condition="input.searchpanel==2",
           h4("Define the categories that you wish to assign each sample (GSM) for comparison."),
           textInput("cat1", label=NULL, placeholder="Category 1"),
           textInput("cat2", label=NULL, placeholder="Category 2"),
           textInput("cat3", label=NULL, placeholder="Category 3 (optional)"),
           ### https://www.aridhia.com/blog/the-sky-is-not-the-limit-embedding-raw-html-and-javascript-to-create-dynamic-ui-elements-in-shiny-applications/   
           ### ^ this should help with dynamically adding/subtracting categories
           
           hr(),
           h4("Highlight samples, then click to Assign them to the specificed category."),
           uiOutput("categorySelect"),
           actionButton("Assign", "Assign GSM to Category")),
           
           conditionalPanel(condition="input.searchpanel==3",
               h4("Please check that your samples are appropriately categorized."),
               actionButton("downloadCEL","Download CEL files"))
       ),
       
       mainPanel(
        tabsetPanel(
        tabPanel("Search for GEO data series (GSE)", value=1,
            DT::dataTableOutput("filteredgse")
        ),
        # Assign samples to categories ------------------------------------------------------
        tabPanel("Assign samples (GSM) to categories", value=2,
            
            DT::dataTableOutput("gsm_table")
            
            # ,
            # # Category assignment panel
            # absolutePanel(id="assignCat",class="panel panel-default",fixed=TRUE,draggable=TRUE,top=120,left="auto",right=20,bottom="auto",width="420",height=50,
            # wellPanel(
            #     h4("Highlight samples, then click to Assign them to the specificed category."),
            #     uiOutput("categorySelect"),
            #     actionButton("Assign", "Assign GSM to Category")
            #     )
            # )
        ),
            # This will be where the CEL files are downloaded (confirmation, etc) ------------
        tabPanel("Confirm sample categories", value=3,
                DT::dataTableOutput("finishedtable")
        ),
        
        id = "searchpanel"
        )
        
        )
        )
        
    ),
    
    # Load Gene Expression Data tab -------------------------------------
    tabPanel("Load expression data",
        h3("Load experiments to perform analysis"),
        hr(),
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
        h3("Compare genes based on absolute expression"),
        hr(),
        sidebarLayout(
        sidebarPanel(
            # style = "position:fixed",
            checkboxGroupInput("tissues", label = "Select tissues to inclued",
            choices = c("photoreceptors","RPE","whole.retina"), selected = c("photoreceptors","RPE","whole.retina")
            ),
            br(),
            checkboxInput("de_state", label = "Show differential expressed only", value = TRUE),
            uiOutput("de_choices"),
            br(),
            conditionalPanel(condition="input.absexpanel==1",
                h3("Heatmap parameters"),
                checkboxInput("hm_probes", "Show probe-level", value = FALSE),
                checkboxInput("hm_gsm", "Show GSM (column names)", value = TRUE),
                checkboxInput("hm_rownames", "Show rownames", value = TRUE),
                checkboxInput("hm_col_cluster", "Cluster columns", value = TRUE),
                checkboxInput("hm_row_cluster", "Cluster rows", value = TRUE),
                h4("Select plot dimensions (px)"),
                numericInput("hm_width", "Width", value = 900, min = 100, max = 2400, step = 10),
                numericInput("hm_height", "Height", value = 1200, min = 100, max = 2400, step = 10))
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
    tabPanel("Relative expression",
        h3("Compare genes based on relative expression between experimental groups"),
        hr(),
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