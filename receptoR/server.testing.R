########################################
#$#$#$#$#$#$    Shiny App  $#$#$#$#$#$#$
########################################


## SERVER
server <- function(input, output,session) {
  ## Change GPL
  output$gplSelection <- renderText({
    species <- switch(input$gplSelection,
                   mouse = 'GPL1260',
                   human = 'GPL570')
    paste("You chose", species)
  })
  
  ## Search functions
  Totalchar <- eventReactive(input$Search, {nchar(input$Key)})
  Commas <- eventReactive(input$Search, {which(strsplit(input$Key, "")[[1]]==",")})
  Ncommas <- eventReactive(input$Search, {length(Commas())})
  Commasstart <- eventReactive(input$Search, {Commas() + 1})
  Commasend <- eventReactive(input$Search, {Commas() - 1})
  
  Searchterms <- eventReactive(input$Search, {
    substring(input$Key, c(1, Commasstart()), c(Commasend(), Totalchar()))
  })
  
  filtered_gse <- eventReactive(input$Search, {
     
  #### This worked here (and made the checkbox selection work faster) but all the *.df variables are within this eventReactive scope, so can't be seen by the other functions.
  #### I definitely need to work on filtering the SQL database as much as possible, since this will be persistent, then slapping the results into a data frame and ultimately a RDA object for long term storage / archiving.
       
      ## Filter the SQL database *before* converting to data frames
      #gse.gpl.filter = gse %>% 
      #  left_join(gse_gpl, copy = TRUE) %>% 
      #  filter(gpl == species) %>% 
      #  collect()
        
      ## Convert SQLite to data frames
      withProgress(message = 'Loading GEO Database', value = 0, {
          as.data.frame.DataTable(gse) -> gse.df
          incProgress(1/3)
          as.data.frame.DataTable(gse_gsm) -> gse_gsm.df ## should be filtered first!
          incProgress(1/3)
          as.data.frame.DataTable(gsm) -> gsm.df ## should be filtered first!
          incProgress(1/3)
          gsm.df$category <- rep("Not yet assigned", nrow(gsm.df)) # I don't love this but it should get the job done
          gse_to_filter <- data.frame(gse="")
      })
      dplyr::filter(gse.df, str_detect(gse.df$title, Searchterms()))
  })
  
  ## Collect samples to use (GSE - GSM)
    # List of the GSM associated with the selected GSE
  
  gse_to_keep <- eventReactive(input$GSE_GSM, {
    filtered_gse()[input$filteredgse_rows_selected,c(3,7)]
  })
  
  # Use GSE to load GSM
  gsm_annotated <- eventReactive(input$GSE_GSM, {
     dplyr::filter(gsm.df,series_id %in% gse_to_keep()$gse)
  })
  
  ## Assign categories to each sample (GSM)
  # Assign categories
  rows <- reactiveValues() 
      observeEvent(input$Assign, {
          if (input$Assign == 1) {
            gsm_selected <- gsm_annotated()
            gsm_selected[input$gsm_table_rows_selected,"category"] <- input$selection
            rows$df <- gsm_selected
            gsm_selected <<- rows$df # '<<-' is necessary to get this to the enclosing environment
          }
          else
          {
            gsm_selected[input$gsm_table_rows_selected,"category"] <- input$selection
            rows$df <- gsm_selected
            gsm_selected <<- rows$df
          }
      })

  finishedtable <- eventReactive(input$Remove, {
    dplyr::filter(rows$df, category %in% c(input$cat1, input$cat2, input$cat3))
  })
 
  ## Outputs
  output$page3 <- renderUI(
    fluidRow(
      column(3,
             selectInput("selection", "Select a Category",
                         c("category1" <- {input$cat1},
                           "category2" <- {input$cat2},
                           "category3" <- {input$cat3},
                           "category4" <- "Not included"))
      )
    )
  )
  
  output$filteredgse <- DT::renderDataTable({
      filtered_gse()[,c(1,2,7)]}, options=list(searching=TRUE, pageLength=20))
 
  output$GSEtoGSMlist <- renderTable(
    if (input$GSE_GSM == 0)
      return ()
    else
      return (filter(gse_gsm.df,gse %in% gse_to_keep()$gse)))
  
  output$gsm_table <- DT::renderDataTable({
    if (input$Assign == 0)
      return (gsm_annotated()[,c(2,3,27,33)])
    else
      return (gsm_annotated()[,c(2,3,27,33)])}, options=list(searching=FALSE))

  output$finishedtable <- renderTable({finishedtable()})
      
  ## Kill shinyApp when session closes
  session$onSessionEnded(stopApp)

}
