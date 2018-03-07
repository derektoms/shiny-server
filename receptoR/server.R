# 2018-03-01

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

## Replace the following line with periodic updates of the db from GEO
#if(!file.exists('./../data/GEOmetadb.sqlite')) getSQLiteFile()

GEOdb = src_sqlite('../../data/GEOmetadb.sqlite')
#src_tbls(GEOdb)
gse = tbl(GEOdb, 'gse')
gse_gpl = tbl(GEOdb, 'gse_gpl')
#gpl = tbl(GEOdb, 'gpl') 
gsm = tbl(GEOdb, 'gsm')
gse_gsm = tbl(GEOdb, 'gse_gsm')


# Filter GEOdb tables based on GPL that can be used for subsequent searches based on the species selected (this greatly speeds up the actual run time). I thought it best to keep the whole table, then select specific columns as needed farther along.

## GPL570
gseGPL570 <- gse %>% 
left_join(gse_gpl, copy = TRUE) %>% 
filter(gpl == "GPL570") %>% 
collect()

gsmGPL570 <- gsm %>% 
filter(gpl == "GPL570") %>% 
collect()

## GPL1261
gseGPL1261 <- gse %>% 
left_join(gse_gpl, copy = TRUE) %>% 
filter(gpl == "GPL1261") %>% 
collect()

gsmGPL1261 <- gsm %>% 
filter(gpl == "GPL1261") %>% 
collect()


########################################
#$#$#$#$#$#$    Shiny App  $#$#$#$#$#$#$
########################################


## SERVER
server <- function(input, output,session) {
    
  ## Change GPL
  output$gplSelection <- renderText({
    species <- switch(input$gplSelection,
                   mouse = 'GPL1261',
                   human = 'GPL570')
    paste("You chose", species)
  })
  
#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$
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
      if(input$gplSelection=='human'){
          dplyr::filter(gseGPL570, str_detect(gseGPL570$title, Searchterms()))
      } else {
          dplyr::filter(gseGPL1261, str_detect(gseGPL1261$title, Searchterms()))
      }
      
  })

#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$
  
  ## Collect samples to use (GSE - GSM)
    # List of the GSM associated with the selected GSE
  
  gse_to_keep <- eventReactive(input$GSE_GSM, {
    filtered_gse()[input$filteredgse_rows_selected,]
  })
  
  # Use GSE to load GSM
  gsm_annotated <- eventReactive(input$GSE_GSM, {
      if(input$gplSelection=='human'){
          dplyr::filter(gsmGPL570,series_id %in% gse_to_keep()$gse)
      } else {
          dplyr::filter(gsmGPL1261,series_id %in% gse_to_keep()$gse)
      }
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
  
 # output$filteredgse <- DT::renderDataTable({
 #     filtered_gse()[,c(1,2,7)]}, options=list(searching=TRUE, pageLength=20))
 
  output$filteredgse <- DT::renderDataTable({
          filtered_gse()}, options=list(searching=TRUE, pageLength=20))
 
 
  output$GSEtoGSMlist <- renderTable(
    if (input$GSE_GSM == 0)
      return ()
    else
      return (filter(gse_gsm,gse %in% gse_to_keep()$gse)))
  
  output$gsm_table <- DT::renderDataTable({
    if (input$Assign == 0)
      return (gsm_annotated())
    else
      return (finishedtable())}, options=list(searching=FALSE))

  output$finishedtable <- renderTable({finishedtable()})
      
  ## Kill shinyApp when session closes
  session$onSessionEnded(stopApp)

}
