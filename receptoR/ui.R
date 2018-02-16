## Javascript
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
  #creation of a navigation bar and mulitple pages
  navbarPage("Bioinformatics Software",
             tabPanel("Search for GEO data series (GSE)",  
                      #search GSE, and select which to include
                      helpText("After searching, click on the second tab to proceed to the next page"),
                      radioButtons("gplSelection", "Choose species:", choices = c("Mouse (GPL1260)" = "mouse", "Human (GPL570)" = "human")),
                      textOutput("gplSelection"),
                      tagAppendAttributes(
                        textInput("Key", "Enter search terms, separated by commas", value = ""),
                        `data-proxy-click` = "Search"
                      ),
                      actionButton("Search", "Search")
             ),
             tabPanel("Select GEO data series (GSE)", uiOutput("page1"), 
                      helpText("Highlight the desired search results (GSE) and click 'Retrieve GSM' to proceed"),
                      actionButton("GSE_GSM", "Retrieve GSM"),
                      helpText("Do not click 'finish' until all selections have been made. 
                               This button removes the unselected rows and generates a new table on the next page."),
                      DT::dataTableOutput("filteredgse"),
                      tableOutput("GSEtoGSMlist")
                      
             ),
             tabPanel("Define categories for GEO samples (GSM)", uiOutput("page2"), 
                      helpText("Define the categories that you wish to compare. 
                               After this is complete, click on the third tab to proceed to the next page"),
                      textInput("cat1", "Define Category 1"),
                      textInput("cat2", "Define Category 2"),
                      textInput("cat3", "Define Category 3")
             ),
             tabPanel("Assign samples to categories", uiOutput("page3"), 
                      helpText("Highlight the desired search results and click 'assign' to assign them to the specificed category"),
                      actionButton("Assign", "Assign Categories"),
                      verbatimTextOutput("selectedRows"),
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
