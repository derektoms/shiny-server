## 2018-02-27

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

library(shiny)
library(shinythemes)
library(shinyjs)

## UI
ui <- fluidPage(
  tags$head(tags$script(HTML(jscode))),
  #creation of a navigation bar and mulitple pages
  navbarPage("receptoR",
  
              theme = "sandstone.css",
    # Search for GSE  ------------------------------------------------------------------------------              
            tabPanel("Select GEO data series (GSE)",  #search GSE, and select which to include
                sidebarLayout(
                    sidebarPanel(
                        radioButtons("gplSelection", "Choose species:", choices = c("Mouse (GPL1261)" = "mouse", "Human (GPL570)" = "human")),
                        tagAppendAttributes(
                            textInput("Key", "Enter search terms, separated by commas", value = ""),
                            `data-proxy-click` = "Search"
                        ),
                        actionButton("Search", "Search"),
                        br(),
                        hr(),
                        helpText("Define the categories that you wish to assign each sample (GSM) for comparison."),
                        textInput("cat1", "Define Category 1"),
                        textInput("cat2", "Define Category 2"),
                        textInput("cat3", "Define Category 3")
                    ),
                    mainPanel(
                        uiOutput("page1"), 
                        helpText("Highlight the desired search results (GSE) and click 'Retrieve GSM' to proceed"),
                        actionButton("GSE_GSM", "Retrieve GSM"),
                        helpText("Do not click 'finish' until all selections have been made. This button removes the unselected rows and generates a new table on the next page."),
                        DT::dataTableOutput("filteredgse"),
                        tableOutput("GSEtoGSMlist")
                    )
                )
            ),
              tabPanel("Assign samples to categories", uiOutput("page3"), 
                      helpText("Highlight the desired search results and click 'assign' to assign them to the specificed category"),
                      actionButton("Assign", "Assign Categories"),
                      verbatimTextOutput("selectedRows"), ## doesn't seem to be working
                      actionButton("Remove", "Finalize selections and remove not included"),
                      helpText("Do not click 'finish' until all selections have been made. 
                               This button removes the unselected rows and generates a new table on the next page."),
                      DT::dataTableOutput("gsm_table")
             ),
             tabPanel("Selection details", uiOutput("page4"), 
                      tableOutput("finishedtable")
             )
  )
)
