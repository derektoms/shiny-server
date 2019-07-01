#                          _        ____
#  _ __ ___  ___ ___ _ __ | |_ ___ |  _ \
# | '__/ _ \/ __/ _ \ '_ \| __/ _ \| |_) |
# | | |  __/ (_|  __/ |_) | || (_) |  _ <
# |_|  \___|\___\___| .__/ \__\___/|_| \_\
#                   |_|
#
# June 2019 receptoR v 1.3
## Last update: 2019-06-22, Derek Toms
## server.R


########################################
#$#$#$#$#$#$    HEADER     $#$#$#$#$#$#$
########################################

# App structural packages:
#install.packages(c('dplyr','dbplyr','tidyr','ggplot2','RColorBrewer','readr','stringr','shiny','shinythemes','shinyjs','DT'))

library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(readr)
library(stringr)
library(shiny)
library(shinythemes)
library(shinyjs)
library(dbplyr)
library(DT)

# Bioinformatics packages installed via biocLite:
#source("https://bioconductor.org/biocLite.R")
#biocLite(c('limma','annotate','genefilter','ComplexHeatmap','pheatmap','cowplot','GEOmetadb','mouse4302.db','hgu133plus2.db'))
#biocLite(c('mixOmics','MergeMaid','GEOquery','inSilicoMerging','affy','sva','Rtsne','metaArray','testthat'))

library(GEOmetadb)
library(GEOquery)
library(affy)

library(limma)
library(annotate)
library(pheatmap)
library(mixOmics)
library(cowplot)

library(pool)
library(writexl)

# Microarray platform annotations:
library(mouse4302.db)
library(hgu133plus2.db)
library(hgu133plus2cdf) # 2019-04-14 required
library(mouse4302cdf)

source("functions.R")
## 2019-03-27 Ran this to get the latest database
# if(!file.exists('./../data/GEOmetadb.sqlite')) getSQLiteFile()
load("./../2019-04_genelists.rda")

#------------------------------------------------------------------------------------+
### for local work
# load("~/Documents/Retina/CNIB_TuckMacPhee/Bioinformatics/2018-12_genelists.rda")
# poolGEO <- dbPool(
#   drv = RSQLite::SQLite(),
#   dbname = "/Volumes/ULTRA/across_array/GEOmetadb.sqlite"
# )
#
# poolUserData <- dbPool(
#   drv = RSQLite::SQLite(),
#   dbname = "~/Documents/Retina/CNIB_TuckMacPhee/Bioinformatics/2019-06-15 v1.3 Update/receptoRUserData.sqlite"
# )
#------------------------------------------------------------------------------------+

# 2019-03-04
## Connection to GEO Metadata DB
poolGEO <- dbPool(
  drv = RSQLite::SQLite(),
  dbname = "./data/GEOmetadb.sqlite"
)

poolUserData <- dbPool(
  drv = RSQLite::SQLite(),
  dbname = "./data/receptoRUserData.sqlite"
)

onStop(function() {
  poolClose(poolGEO)
  poolClose(poolUserData)
})

## Initialize user experiments to load
global <- reactiveValues (DatasetTable = loadUserDatasets(poolUserData))

########################################
#$#$#$#$#$#$#    SERVER    #$#$#$#$#$#$#
########################################

server <- function(input, output, session) {
    
#   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-
#  / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \
# `-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'
#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$
## This is the database search begins
#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$

# Quick link from the main page
observeEvent(input$linkSearch, {
  updateNavbarPage(session, "receptorMain", selected="searchPanel")
})


# Set up colour environment
#_,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,_
  catCol <- brewer.pal(3, "Set1")
  rowCol <-desat(catCol)
  userID <- NULL
 
# Search functions 
#_,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,_
  
  searchGSM <- eventReactive(input$searchButton, {
      if(input$gplSelection=='human'){
          sql<-"SELECT * FROM appgsm WHERE description MATCH ?id1 AND gpl LIKE 'GPL570';"
      } else {
          sql<-"SELECT * FROM appgsm WHERE description MATCH ?id1 AND gpl LIKE 'GPL1261';"
      }
      query<-sqlInterpolate(poolGEO,sql,id1=input$searchText)
      queryGSM<-dbGetQuery(poolGEO,query)
      shinyjs::disable("gplSelection")
      return(queryGSM)
  })

  output$searchResultsGSM <- DT::renderDataTable({
          searchGSM()}, extensions = 'Buttons', options=list(
              dom = 'Bfrtip',
              buttons = list(list(extend = 'colvis')),
              searching=TRUE, 
              paging=FALSE,
              scrollX=TRUE, 
              scrollY='60vh', 
              scrollCollapse=TRUE,
              fixedHeader=TRUE,
              autoWidth=TRUE,
              columnDefs=list(list(
              targets = "_all",
              render = JS(
                  "function(data, type, row, meta) {",
                      "return type === 'display' && typeof data === 'string' && data.length > 100 ?",
                      "'<span title=\"' + data + '\">' + data.substr(0, 100) + '...</span>' : data;",
                      "}") 
                      )))) ## typeof data needs to be a string, as a "NA" converted to JS "NULL" breaks things

# Set up tables to store user-selected data
#_,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,_

  proxy.search = dataTableProxy('searchResultsGSM')

  ## Set up reactive table to store experimental samples
  userSamples <- reactiveValues()
  userSamples$df <- data.frame()
    
# Add sample (array) record to the current experiment 
#_,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,_
 
  observeEvent(input$addButton, {
      gsm_selected <- searchGSM()[input$searchResultsGSM_rows_selected,]
      gsm_selected$category <- rep("Not yet assigned", nrow(gsm_selected))
      userSamples$df <<- rbind(userSamples$df,gsm_selected)
      proxy.search %>% selectRows(NULL)
      updateTabsetPanel(session = session, inputId = "searchpanel", selected = "2")
  })

# Assign categories to each sample (GSM)
#_,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,_
  
  observeEvent(input$assignButton, {
        userSamples$df[input$gsm_table_rows_selected,"category"] <<- input$selection
  })      
   
  output$gsm_table <- DT::renderDataTable({
      if(input$assignButton == 0){
         return (datatable(userSamples$df, extensions = 'Buttons', options=list(
               dom = 'Bfrtip',
               buttons = list(list(extend = 'colvis')),
               searching=TRUE, 
               paging=FALSE,
               scrollX=TRUE, 
               scrollY='60vh', 
               scrollCollapse=TRUE,
               fixedHeader=TRUE,
               autoWidth=TRUE,
               columnDefs=list(list(
             targets = "_all",
             render = JS(
                 "function(data, type, row, meta) {",
                     "return type === 'display' && typeof data === 'string' && data.length > 100 ?",
                     "'<span title=\"' + data + '\">' + data.substr(0, 100) + '...</span>' : data;",
                     "}")
                     )))))
      } else {
         return (datatable(userSamples$df, extensions = 'Buttons', options=list(
             dom = 'Bfrtip',
             buttons = list(list(extend = 'colvis')),
             searching=TRUE, 
             paging=FALSE,
             scrollX=TRUE, 
             scrollY='60vh', 
             scrollCollapse=TRUE,
             fixedHeader=TRUE,
             autoWidth=TRUE,
             columnDefs=list(list(
             targets = "_all",
             render = JS(
                 "function(data, type, row, meta) {",
                     "return type === 'display' && typeof data === 'string' && data.length > 100 ?",
                     "'<span title=\"' + data + '\">' + data.substr(0, 100) + '...</span>' : data;",
                     "}")
                     )))) %>%
                     formatStyle('category', target="row", backgroundColor=styleEqual(c(input$cat1, input$cat2, input$cat3), c(rowCol[1], rowCol[2], rowCol[3]))))
      }
  })
                 
  proxy.gsm = dataTableProxy('gsm_table')
 
  observeEvent(input$assignButton,{
      proxy.gsm %>% selectRows(NULL)
  }) 
  
  
  ## UI output

    output$categorySelect <- renderUI(
      fluidRow(
        column(12,
               selectizeInput("selection", "Select a Category",
                           c("category1" <- {input$cat1},
                             "category2" <- {input$cat2},
                             "category3" <- {input$cat3},
                             "category4" <- "Not included")
                             # , options = list(create=TRUE, plugins = list("remove_button")))  ### <- "remove_button" isn't what I thought it was. I would also like the "create" option but I will need to link this to the table as cat1-3 are linked (otherwise new variables are not coloured or sent along for processing)
        )
      )     ### 2018-12-10 I'd like to have a button to add category 3
    )
    )  


# Finished table, to ultimately lead to CEL download
#_,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,_
userSamples$finishedtable <- NULL

observeEvent(input$assignButton, {
    userSamples$finishedtable <<- dplyr::filter(userSamples$df, category %in% c(input$cat1, input$cat2, input$cat3))
})
 
  output$finishedtable <- DT::renderDataTable({
      if(!is.null(userSamples$finishedtable)){
      datatable(userSamples$finishedtable,
      options=list(
          searching=FALSE, 
          paging=FALSE,
          scrollX=TRUE, 
          scrollY='60vh', 
          scrollCollapse=TRUE,
          fixedHeader=TRUE,
          autoWidth=TRUE,
          columnDefs=list(list(
          targets = "_all",
          render = JS(
              "function(data, type, row, meta) {",
                  "return type === 'display' && typeof data === 'string' && data.length > 100 ?",
                  "'<span title=\"' + data + '\">' + data.substr(0, 100) + '...</span>' : data;",
                  "}")
          )))) %>%
      formatStyle('category',target="row",
      backgroundColor=styleEqual(c(input$cat1,input$cat2,input$cat3),c(rowCol[1],rowCol[2],rowCol[3]))
  )}})

rv <- reactiveValues(download_flag = 0)


  output$report <- downloadHandler(
      filename = function(){paste(input$downloadId,"GSM_report.csv",sep="_")},
      content = function(file){
          write.csv(userSamples$finishedtable,file)
          rv$download_flag <- rv$download_flag + 1
      })

# Modal confirming CEL download, and processing function
#_,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,_     
observeEvent(input$downloadCEL, {
    userSamples$finishedtable %>% group_by(category) %>% summarise(n.gse = n_distinct(series_id)) -> gse.check
    warning <- "Please click below to begin processing the data."
    numCat <- length(gse.check$category)>1
    if(length(which(gse.check$n.gse==1))!=0){
        catAlert <- paste(gse.check$category[which(gse.check$n.gse==1)], collapse = ", ")
        warning <- paste("WARNING: The following categories contain samples from a single experiment (GSE) and as such they will be confounded by batch effects: ",catAlert,".<br>Please proceed with caution or cancel and select additional samples to add to these categories.",sep="")
    }
    if(!numCat){
        showModal(modalDialog(title="Error! A minimum of two categories are needed.","Experimental samples need to be organized into 2 or 3 categories for appropriate downstream analysis. If you are interested in only one type of sample, we suggest choosing samples to act as 'background', which will allow for differential analysis to identify which receptor genes are enriched or depleted in your sample of interest.",
        easyClose = TRUE,
        footer = tagList(
            modalButton("Cancel")))) 
    } else {
        showModal(modalDialog(title="Important! Downloading raw .CEL files from the NCBI server.",HTML(paste("June 20th, 2019<br>",warning)),
        easyClose = TRUE,
        footer = tagList(
            modalButton("Cancel"),
            actionButton("process","Proceed"))))      
    }
  })

observeEvent(input$process, {
    shinyjs::disable("process")
    userID <<- processData(userSamples$finishedtable, input$downloadId, input$comments, input$gplSelection, poolUserData)
    global$DatasetTable <<- loadUserDatasets(poolUserData)
    removeModal()
    showModal(modalDialog(title="Your dataset was successfully processed!","Analyse your data in the 'Load Expression Datasets' tab. You can also download a report from this page.",
    easyClose = TRUE,
    footer = tagList(
        modalButton("OK"))))# modal
  })

# Reset button, modal confirmation
#_,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,_
  observeEvent(input$linkReset, {
      showModal(modalDialog(title="Important! Are you sure you want to reset everything?","All searches and categorized samples will be lost. This can not be undone.",
      footer = tagList(
          modalButton("Cancel"),
          actionButton("buttonReset","Yes, reset."))))# modal
      # confirm reset (all categories, sample search, gone)
      observeEvent(input$buttonReset, {
          shinyjs::enable("gplSelection")
          userSamples$df <<- userSamples$df[0,]
          reset("searchText")
          reset("cat1")
          reset("cat2")
          reset("cat3")
          reset("downloadId")
          replaceData(proxy.search, NULL)
          replaceData(proxy.gsm, NULL)
          userSamples$finishedtable <<- NULL
          removeModal()
          updateTabsetPanel(session = session, inputId = "searchpanel", selected = "1")
        })

    
  })
  


#   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-
#  / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \
# `-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'

#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$
## This is where the analysis begins
#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$

# Quick link from the main page
observeEvent(input$linkLoad, {
  updateNavbarPage(session, "receptorMain", selected="expressionPanel")
})

# Conditional nav tabs
hideTab(inputId = "receptorMain", target = "Gene-level Expression")
hideTab(inputId = "receptorMain", target = "Sample-level Expression")

observeEvent(input$user_data,{
    if(input$user_data!="none"){
        showTab(inputId = "receptorMain", target = "Gene-level Expression")
        showTab(inputId = "receptorMain", target = "Sample-level Expression")
    } else {
        hideTab(inputId = "receptorMain", target = "Gene-level Expression")
        hideTab(inputId = "receptorMain", target = "Sample-level Expression")

    }
})

# Load dataset
#_,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,_
output$loadUserExperiments = renderUI({
    selectizeInput(inputId="user_data",label="Select an experiment for analysis",choices=c("none"="none",split(global$DatasetTable$desc, global$DatasetTable$species)),selected="none")
})


observeEvent(input$user_data,{
    id <- NULL
    datasetToLoad <- NULL
   
   if(input$user_data=="none"){
        mapped_probes<<-NULL
        eset<<-NULL
        de_choices<<-NULL
        sig_genes_lfc<<-NULL
    }else{         
        id <- global$DatasetTable$userID[which(global$DatasetTable$desc == input$user_data)]
        assign(
              x = "species", value = global$DatasetTable$species[which(global$DatasetTable$desc == input$user_data)], envir = .GlobalEnv
           )
        datasetToLoad <- paste("./data/app_data_", id, ".rda", sep='')
        withProgress(message="Loading dataset",value=0.2,{
            load(datasetToLoad,envir=.GlobalEnv)
            
            incProgress(0.3, message = "Loading contrasts")
            
            updateCheckboxGroupInput(session, "tissues", choices = groups, selected = groups)
            updateCheckboxGroupInput(session, "pls_tissues", choices = groups, selected = groups)
            updateCheckboxGroupInput(session, "de", choices = de_choices, selected = de_choices[1])
            
            incProgress(0.3, message ="Loading genelists")
            updateCheckboxGroupInput(session, "genelist", label = NULL, choices = names(gene_lists[[species]]), selected = NULL, inline = FALSE)
            
            incProgress(0.2, message = "Loading gene names")
            updateSelectInput(session, "gene", choices = all_genes[species])
        })
        
    }
    
})

# Download DEG
#_,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,_

rvDEG <- reactiveValues(download_flag = 0)

  # proxy.finishedtable = dataTableProxy('finishedtable')
  output$reportDEG <- downloadHandler(
      filename = paste(input$user_data,"DEG_report.xlsx",sep="_"),
      # filename = paste(input$user_data,"DEG_report.csv",sep="_"),
      content = function(file){
          write_xlsx(sig_genes_lfc, path=file)
          rvDEG$download_flag <- rvDEG$download_flag + 1
      })

# QC output
#_,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,_

output$QC = renderUI({
  validate(
    need(input$user_data!="none","No dataset selected")
  )

  id <- global$DatasetTable$userID[which(global$DatasetTable$desc == input$user_data)]

  fluidRow(
      inlineCSS(list(
          "#norm" = c("max-width:100%","width=100%"))),
      h4("Expression normalization (array intensity, before and after)"), 
      tags$div(class="norm",
          tags$img(src=paste("array_normalization_", id, ".png", sep=''))
          ),
      h4("RNA degradation plot (probe position along transcript vs intensity)"),
      tags$img(src=paste("probe_degradation_", id, ".png", sep='')))
})

      
# Load genes tab
#_,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,_

  geneList = reactive({
    if (is.null(input$genelist) && is.null(input$gene)) {
      return(NULL)
    }
    
    genes = c()

    if (!is.null(input$genelist)) {
      for (gene in input$genelist) {
        genes = c(genes, gene_lists[[species]][[gene]])
      }
    }

    if (!is.null(input$gene)) {
      genes = c(genes, input$gene)
    }
    
    return(unname(genes))
  })
  
 summary_gene_data = reactive({
     validate(
       need(input$user_data!="none","No dataset selected. Please select an experiment for analysis."),
       need(geneList(), "No genes selected. Please select receptor type(s) to analyse.")
     )
   get_expression_summary(eset, geneList())
 })


# Gene outputs
#_,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,_
  
  output$genes = DT::renderDataTable({
      validate(
        need(input$user_data!="none","No dataset selected. Please select an experiment for analysis in 'Load Expression Data'."),
        need(geneList(), "No genes selected. Please select receptor type(s) to analyse in 'Load Expression Data'.")
      )
    
     summary_gene_data() %>% datatable() %>% 
      formatRound(2:4)
    
  })
  
  # single gene plot
 output$singleGenePlot = renderPlot({
     validate(
       need(input$user_data!="none","No dataset selected."),
       need(geneList(), "No genes selected."),
       need(input$genes_rows_selected >= 1, "Please select one or more genes from the 'Average Expression' table to inspect expression by tissue type.")
     )
    
    rows = as.integer(input$genes_rows_selected)
    genes_to_plot = summary_gene_data()$Symbol[rows]
    
    gene_data = get_gene_data(eset, genes_to_plot)
    by_gene_violplot(gene_data,tissues=groups)
    
    
  })

# Expression tab
#_,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,_
  observe({
    toggle("de_choices", anim = TRUE, condition = input$de_state )
  })
  
  genesToPlot = reactive({
    validate(
      need(input$user_data!="none","No dataset selected. Please select an experiment for analysis in 'Load Expression Data'."),
      need(geneList(), "No genes selected. Please select receptor type(s) to analyse in 'Load Expression Data'.")
    )

    genes = geneList()
    
    if(input$de_state) {
      selected_de = input$de
      de_lists = lapply(selected_de, function(x) { as.character(get_de_genes(genes, x, sig_genes_lfc)$Symbol) })
      genes = Reduce(union, de_lists)
    } 
   
    return(genes) 
  }) 

  
  
# Heatmap plot
#_,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,_
  output$expressionPlot = renderPlot({
      validate(
          need(input$user_data!="none","No dataset selected. Please select an experiment for analysis in 'Load Expression Data'."),
          
          need(geneList(), "No genes selected. Please select receptor type(s) to analyse in 'Load Expression Data'."),
          
          need(input$tissues, "No tissues selected. Please choose at least one tissue to plot receptor heatmap."),
          
          need(length(genesToPlot())>10, if(input$de_state){paste("Based on the genes selected in 'Load Data', ", length(genesToPlot())," genes were differentially expressed in these tissues (",paste(input$tissues, collapse = ", "), "); try unselecting that option in the side menu.",sep="")}else{paste("No genes to plot as a heatmap (minimum = 10). Try including more receptor types in 'Load Data'.")})
    )
       
    selected_tissues = input$tissues
    sub_eset = eset[, eset$tissue %in% selected_tissues]
    genes = gene2probe(genesToPlot(), mapped_probes)
    
    cat(file=stderr(), "Preparing heatmap:\n Tissues:", paste(input$tissues, collapse = ", "), "\n gene list: ",paste(genesToPlot(),collapse=", "),"\n genes: ", paste(genes,collapse=", "),"\n")
    
    gene_heatmap(sub_eset, genes, scale = "row",
                  probe_level = input$hm_probes,
                  gsm_show = input$hm_gsm,
                  show_rownames = input$hm_rownames,
                  cluster_rows = input$hm_row_cluster,
                  cluster_cols = input$hm_col_cluster,
                  border_color = NA)
    })
  
  output$heatmap_ui = renderUI({
    plotOutput("expressionPlot", height = input$hm_height, width = input$hm_width)
  })

# Overall expression
#_,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,_
  output$overallPlot = renderPlot({
    validate(
        need(input$user_data!="none","No dataset selected. Please select an experiment for analysis in 'Load Expression Data'."),
        need(geneList(), "No genes selected. Please select receptor type(s) to analyse in 'Load Expression Data'."),
    need(input$tissues, "No tissues selected. Please choose at least one tissue to plot receptor heatmap."),
    need(length(genesToPlot())>1, if(input$de_state){paste("Based on the genes selected in 'Load Data', ", length(genesToPlot())," genes were differentially expressed in these tissues (",paste(input$tissues, collapse = ", "), "); try unselecting that option in the side menu.",sep="")}else{paste("No genes to plot. Try including more receptor types in 'Load Data'.")}) 
    )
    
    gene_data = get_gene_data(eset, genesToPlot())
    overall_expression_boxplot(gene_data, tissues = input$tissues)
    
  })
  

# By gene boxplots ----------------------------------------------------------------------------

  output$byGenePlot = renderPlot({
      validate(
          need(input$user_data!="none","No dataset selected. Please select an experiment for analysis in 'Load Expression Data'."),
          need(geneList(), "No genes selected. Please select receptor type(s) to analyse in 'Load Expression Data'."),
      need(input$tissues, "No tissues selected. Please choose at least one tissue to plot receptor heatmap."),
      need(length(genesToPlot())>1, if(input$de_state){paste("Based on the genes selected in 'Load Data', ", length(genesToPlot())," genes were differentially expressed in these tissues (",paste(input$tissues, collapse = ", "), "); try unselecting that option in the side menu.",sep="")}else{paste("No genes to plot. Try including more receptor types in 'Load Data'.")}) 
      )    
    gene_data = get_gene_data(eset, genesToPlot())
    by_gene_boxplot(gene_data, tissues = input$tissues)
  })
  


  plsdaData = reactive({
    
    selected_tissues = input$pls_tissues
    if(length(selected_tissues) < 2) {
      return(NULL)
    }
    
    
    sub_eset = eset[, eset$tissue %in% selected_tissues]
    genes = gene2probe(geneList(), mapped_probes)
    
    # probe = input$pls_probe
    
    get_plsda(sub_eset, genes, probe = FALSE) 
    
  })

# PCA plot ----------------------------------------------------------------------------
  output$indPlot = renderPlot({
    validate(
        need(input$user_data!="none","No dataset selected. Please select an experiment for analysis in 'Load Expression Data'."),
        need(geneList(), "No genes selected. Please select receptor type(s) to analyse in 'Load Expression Data'."),
       need(length(input$pls_tissues) >= 2, "Please select at least two tissues for a PLS-DA plot.")
    )
    
    plotIndiv(plsdaData()$result, ind.names = FALSE, group = factor(plsdaData()$tissue_grps), pch = 16, 
              col.per.group = brewer.pal(3, "Set1")[1:length(input$pls_tissues)], legend = TRUE, cex = 2, ellipse=TRUE, title="Plot of individual arrays",style="graphics")
  })

# Correlation Circle plot ----------------------------------------------------------------------------  
  output$varPlot = renderPlot({
      validate(
          need(input$user_data!="none","No dataset selected. Please select an experiment for analysis in 'Load Expression Data'."),
          need(geneList(), "No genes selected. Please select receptor type(s) to analyse in 'Load Expression Data'."),
         need(length(input$pls_tissues) >= 2, "Please select at least two tissues for a Correlation circle plot.")
      )
      comp = as.integer(input$pls_ncomp)
    plotVar(plsdaData()$result, var.names = list(plsdaData()$varNames), comp.select=comp, cex = 1, overlap=FALSE, col="grey",title="Correlation circle between genes and discriminant components", style="graphics")
    
  })

  output$numGenesUI = renderUI({
    numericInput("pls_num_genes", "Select number of genes to show contributions for", 
                 value = 25, min = 1, max = length(geneList()), step = 1)
  })
  
# Loadings plot ----------------------------------------------------------------------------
  output$contribPlot = renderPlot({
      validate(
          need(input$user_data!="none","No dataset selected. Please select an experiment for analysis in 'Load Expression Data'."),
          need(geneList(), "No genes selected. Please select receptor type(s) to analyse in 'Load Expression Data'."),
         need(length(input$pls_tissues) >= 2, "Please select at least two tissues for a Loadings plot.")
      )

    
    grps = plsdaData()$result$names$colnames$Y
    ndisplay = input$pls_num_genes
    comp = as.integer(input$pls_ncomp)
    plotLoadings(plsdaData()$result, name.var = plsdaData()$varNames, ndisplay = ndisplay, comp = comp, contrib='max', method='mean',legend.color = catCol[1:length(grps)],title=paste("Weight of the top ", ndisplay, " genes contributing to discriminant component ", comp, sep=""),size.title=1)
     
  })
  
#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$  
  ## Kill shinyApp when session closes
  session$onSessionEnded(stopApp)

}
