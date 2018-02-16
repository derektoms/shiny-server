# 2017-07-31

########################################
#$#$#$#$#$#$    HEADER     $#$#$#$#$#$#$
########################################
# Bioinformatics packages installed via biocLite:
#source("https://bioconductor.org/biocLite.R")
#biocLite(c('limma','annotate','genefilter','ComplexHeatmap','pheatmap','cowplot','GEOmetadb','mouse4302.db','hgu133plus2.db'))

library(limma)
library(annotate)
library(genefilter)
library(ComplexHeatmap)
library(pheatmap)
library(cowplot)
library(GEOmetadb)

#biocLite(c('MergeMaid','GEOquery','inSilicoMerging','affy','sva','Rtsne','metaArray','testthat'))
library(MergeMaid)
library(GEOquery)
library(testthat)
library(metaArray)
#library(inSilicoMerging) not available for this version of R
library(Rtsne)
library(sva)
library(affy)

# Microarray platform annotations:
# Equivalent human platform is GPL570 with 127 514 samples
# HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array

#library(mouse4302.db) 
#library(hgu133plus2.db)  Both aren't available

# Other packages:
#install.packages(c('dplyr','dbplyr','tidyr','ggplot2','RColorBrewer','readr','stringr','shiny','shinythemes','shinyjs','DT'))

library(dplyr)
library(dbplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(readr)
library(stringr)
library(DT)
library(shiny)
library(shinythemes)
library(shinyjs)

########################################
#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$
########################################

# Trying out getting and analyzing the data from GEO

#if(!file.exists('./../data/GEOmetadb.sqlite')) getSQLiteFile()

db = src_sqlite('data/GEOmetadb.sqlite')
src_tbls(db)
gse = tbl(db, 'gse')
gse_gpl = tbl(db, 'gse_gpl')
gpl = tbl(db, 'gpl') 
gsm = tbl(db, 'gsm')
gse_gsm = tbl(db, 'gse_gsm')


########################################
#$#$#$#$#$#$    Shiny App  $#$#$#$#$#$#$
########################################


## SERVER
server <- function(input, output,session) {
  ## Change GPL
  #gplSelection <- switch(input$gplSelection,
  #                 mouse = 'GPL1261',
  #                 human = 'GPL570')
  output$gplSelection <- renderText({
    paste("You chose", input$gplSelection)
  })
  
  ## Convert SQLite to data frames
  as.data.frame.DataTable(gse) -> gse.df
  as.data.frame.DataTable(gse_gsm) -> gse_gsm.df
  as.data.frame.DataTable(gsm) -> gsm.df
  gsm.df$category <- rep("Not yet assigned", nrow(gsm.df)) # I don't love this but it should get the job done
  gse_to_filter <- data.frame(gse="")
  
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
