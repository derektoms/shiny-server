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

  navbarPage("receptoR",
  
              theme = "sandstone.css",
              
    # Search for GSE  ------------------------------------------------------------------------------              
            tabPanel("Select GEO data series (GSE)",  #search GSE, and select which to include
                sidebarLayout(
                    sidebarPanel(
                        # Search for datasets ------------------------------------------------------
                        h4("1. Search for GEO data series (GSE)"),
                        radioButtons("gplSelection", "Choose species:", choices = c("Mouse (GPL1261)" = "mouse", "Human (GPL570)" = "human")),
                        tagAppendAttributes(
                            textInput("Key", "Enter search terms, separated by commas", value = ""),
                            `data-proxy-click` = "Search"
                        ),
                        actionButton("Search", "Search"),
                        hr(),
                        # Define categories --------------------------------------------------------
                        h4("3. Define the categories that you wish to assign each sample (GSM) for comparison."),
                        textInput("cat1", "Define Category 1"),
                        textInput("cat2", "Define Category 2"),
                        textInput("cat3", "Define Category 3")
                    ),
                    # Filtered GSE list -----------------------------------------------------------
                    mainPanel(
                        h4("2. Highlight the desired search results (GSE) and click 'Retrieve GSM' to proceed"),
                        actionButton("getGSM", "Retrieve GSM"),
                        helpText("Do not click 'finish' until all selections have been made. This button removes the unselected rows and generates a new table on the next page."),
                        DT::dataTableOutput("filteredgse")
                    )
                )
            ),
            # Assign samples to categories ------------------------------------------------------
              tabPanel("Assign samples (GSM) to categories", 
                      h4("4. Highlight the desired search results and click 'assign' to assign them to the specificed category"),
                      uiOutput("categorySelect"),
                      actionButton("Assign", "Assign Categories"),
                      actionButton("Remove", "Finalize selections and remove not included"),
                      # ^ why aren't these two the same button?
                      helpText("Do not click 'finish' until all selections have been made. 
                               This button removes the unselected rows and generates a new table on the next page."),
                      DT::dataTableOutput("gsm_table")
             ),
             # This will be where the CEL files are downloaded (confirmation, etc) ------------
             tabPanel("Selection details", uiOutput("page4"), 
                      DT::dataTableOutput("finishedtable"),
                      actionButton("downloadCEL","Download CEL files")
             )
  )
)