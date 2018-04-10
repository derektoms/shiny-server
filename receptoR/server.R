# 2018-04-02
# Hoping to work with the filtered table to integrate 
# 1. saving table data (e.g. search results) for reproducibility
# 2. displaying sorted table tibble and associated categories
# 3. .rda of CEL to get that can be loaded (thinking along the lines of a !exist escape if there is no processed datafile)
# 4. aforementioned processed datafile, saved for each user for some length of time

########################################
#$#$#$#$#$#$    HEADER     $#$#$#$#$#$#$
########################################

# App structural packages:
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

# Bioinformatics packages installed via biocLite:
#source("https://bioconductor.org/biocLite.R")
#biocLite(c('limma','annotate','genefilter','ComplexHeatmap','pheatmap','cowplot','GEOmetadb','mouse4302.db','hgu133plus2.db'))

library(GEOmetadb)
library(GEOquery)

library(affy)

library(limma)
library(annotate)
library(genefilter)
library(ComplexHeatmap)
library(pheatmap)
library(cowplot)

#biocLite(c('MergeMaid','GEOquery','inSilicoMerging','affy','sva','Rtsne','metaArray','testthat'))
library(MergeMaid)

library(testthat)
library(metaArray)
#library(inSilicoMerging) ## package ‘inSilicoMerging’ is not available (for R version 3.4.3) 
library(Rtsne)
library(sva)

# Microarray platform annotations:
# Equivalent human platform is GPL570 with 127 514 samples
# HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array

library(mouse4302.db) 
library(hgu133plus2.db)



########################################
#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$
########################################

# Trying out getting and analyzing the data from GEO

## Replace the following line with periodic updates of the db from GEO
#if(!file.exists('./../data/GEOmetadb.sqlite')) getSQLiteFile()

#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$
### UPDATE March 2018
### process and save GEOdb as .rda files to save time loading
### using .rda instead of saveRDS and .rds to maintain the original names

### This should all be replaced by a function that can periodically update the GEOmetadb and subsequently update the .rda files
#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$

#GEOdb = src_sqlite('../../data/GEOmetadb.sqlite')
#src_tbls(GEOdb)
#gse = tbl(GEOdb, 'gse')
#gse_gpl = tbl(GEOdb, 'gse_gpl')
#gpl = tbl(GEOdb, 'gpl') 
#gsm = tbl(GEOdb, 'gsm')
#gse_gsm = tbl(GEOdb, 'gse_gsm')


# Filter GEOdb tables based on GPL that can be used for subsequent searches based on the species selected (this greatly speeds up the actual run time). I thought it best to keep the whole table, then select specific columns as needed farther along.

## GPL570
#gseGPL570 <- gse %>% 
#left_join(gse_gpl, copy = TRUE) %>% 
#filter(gpl == "GPL570") %>% 
#collect()
#save(gseGPL570, file = "../../data/gseGPL570.rda")

#gsmGPL570 <- gsm %>% 
#filter(gpl == "GPL570") %>% 
#collect()
#save(gsmGPL570, file = "../../data/gsmGPL570.rda")

## GPL1261
#gseGPL1261 <- gse %>% 
#left_join(gse_gpl, copy = TRUE) %>% 
#filter(gpl == "GPL1261") %>% 
#collect()
#save(gseGPL1261, file = "../../data/gseGPL1261.rda")

#gsmGPL1261 <- gsm %>% 
#filter(gpl == "GPL1261") %>% 
#collect()
#save(gsmGPL1261, file = "../../data/gsmGPL1261.rda")

### not the best place to put this, but it should work for now
load("../../data/gseGPL570.rda")
load("../../data/gsmGPL570.rda")
load("../../data/gseGPL1261.rda")
load("../../data/gsmGPL1261.rda")
source("../receptoR/functions.R")
## GSM downloads
CELtoDownload <- c("gsm", "category", "timestamp")
#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$

########################################
#$#$#$#$#$#$    Shiny App  $#$#$#$#$#$#$
########################################

## SERVER
server <- function(input, output, session) {
      
#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$
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
          dplyr::filter(gseGPL1261, str_detect(gseGPL1261$summary, Searchterms()))
      }
  })

  output$filteredgse <- DT::renderDataTable({
          filtered_gse()}, options=list(searching=TRUE, pageLength=6, columnDefs=list(list(
              targets = c(8,9,12),
              render = JS(
                  "function(data, type, row, meta) {",
                      "return type === 'display' && typeof data === 'string' && data.length > 100 ?",
                      "'<span title=\"' + data + '\">' + data.substr(0, 100) + '...</span>' : data;",
                      "}") 
                      )))) ## typeof data needs to be a string, as a "NA" converted to JS "NULL" breaks things

 
#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$
  ## Collect samples to use (GSE - GSM)
    # List of the GSM associated with the selected GSE
  
  gse_to_keep <- eventReactive(input$getGSM, {
    filtered_gse()[input$filteredgse_rows_selected,]
  })
  
  # Use GSE to load GSM from the prefiltered lists
  gsm_annotated <- eventReactive(input$getGSM, {
      withProgress(message='Collecting GSM',{
      if(input$gplSelection=='human'){
          dplyr::filter(gsmGPL570,series_id %in% gse_to_keep()$gse)
      } else {
          dplyr::filter(gsmGPL1261,series_id %in% gse_to_keep()$gse)
      }
      })
  })

  ## ^ these two things should be condensed, so that there is one action on the button click


#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$

  ## Assign categories to each sample (GSM)
  
  output$gsm_table <- DT::renderDataTable({

       if(input$Assign==0){
          return (gsm_annotated()[,c(-5:-7,-11,-12,-14:-26,-28:-32)])
       } else {
          return (samples$df[,c(-5:-7,-11,-12,-14:-26,-28:-32)])
       }
  }, options=list(searching=FALSE, columnDefs=list(list(
              targets = "_all",
              render = JS(
                  "function(data, type, row, meta) {",
                      "return type === 'display' && typeof data === 'string' && data.length > 100 ?",
                      "'<span title=\"' + data + '\">' + data.substr(0, 100) + '...</span>' : data;",
                      "}") 
                      ))))
                      
  proxy.gsm = dataTableProxy('gsm_table')
  observeEvent(input$Assign,{
      proxy.gsm %>% selectRows(NULL)
  })
  
  ## UI output

    output$categorySelect <- renderUI(
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

  ## Assign categories
  samples <- reactiveValues()
  samples$df <- data.frame()
  
  observeEvent(input$Assign, {
      if (input$Assign == 1) {
        gsm_selected <- gsm_annotated()
        gsm_selected$category <- rep("Not yet assigned", nrow(gsm_selected))
        gsm_selected[input$gsm_table_rows_selected,"category"] <- input$selection
        samples$df <<- gsm_selected
      }
      else
      {
        samples$df[input$gsm_table_rows_selected,"category"] <<- input$selection
      }
  })      
  
  # ^ don't love this... would like to have the category set without a button click (maybe change to this tab), but it's working for the moment
#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$

## Finished table, to ultimately lead to CEL download

  finishedtable <- eventReactive(input$Remove, {
    dplyr::filter(samples$df, category %in% c(input$cat1, input$cat2, input$cat3))
  })
 
  output$finishedtable <- DT::renderDataTable({finishedtable()[,c(2,3,4,10,31,32,33)]})
    
  proxy.finishedtable = dataTableProxy('finishedtable')
  observeEvent(input$downloadCEL, {
      proxy.finishedtable %>% selectRows(2) %>% selectColumns(2) # selectColumn doesn't work. Instead what I would prefer is a PDF of the final table and a R list of CEL files to download
      saveData(formData()) # success, this will save the tables as variables I can use to download the CEL files
      ####
      #### Move to SQL 
      ####
      withProgress(
          message = "Downloading and processing GSM",
          processData(CELtoDownload$gsm))
      
  })
 
  formData <- eventReactive(input$Remove, {
      time <- strftime(Sys.time(),"%Y.%m.%d %H:%M")
      stamp <- data.frame("timestamp"=rep(time,nrow(finishedtable())))
      CELdl <- c(finishedtable()[,3],finishedtable()[,33],stamp)
      CELdl
  })
 
 output$CELdl <- renderTable(formData()) 
#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$  
  ## Kill shinyApp when session closes
  session$onSessionEnded(stopApp)

}
