#                          _        ____
#  _ __ ___  ___ ___ _ __ | |_ ___ |  _ \
# | '__/ _ \/ __/ _ \ '_ \| __/ _ \| |_) |
# | | |  __/ (_|  __/ |_) | || (_) |  _ <
# |_|  \___|\___\___| .__/ \__\___/|_| \_\
#                   |_|
#
# March 2019 receptoR v 1.2
## Last update: 2019-04-14, Derek Toms
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

## 2018-12-02 not currently needed:
    # library(genefilter)
   #  library(ComplexHeatmap)
   #
   #  library(MergeMaid)
   #  library(testthat)
   #  library(metaArray)
   #  library(Rtsne)
   #  library(sva)

# Microarray platform annotations:
library(mouse4302.db)
library(hgu133plus2.db)
library(hgu133plus2cdf) # 2019-04-14 required
library(mouse4302cdf)

source("functions.R")
## 2019-03-27 Ran this to get the latest database
# if(!file.exists('./../data/GEOmetadb.sqlite')) getSQLiteFile()
load("./../2018-12_genelists.rda")

### for local work
# load("~/Documents/Retina/CNIB_TuckMacPhee/Bioinformatics/gseGPL570.rda")
# load("~/Documents/Retina/CNIB_TuckMacPhee/Bioinformatics/gsmGPL570.rda")
# load("~/Documents/Retina/CNIB_TuckMacPhee/Bioinformatics/gseGPL1261.rda")
# load("~/Documents/Retina/CNIB_TuckMacPhee/Bioinformatics/gsmGPL1261.rda")
# load("~/Documents/Retina/CNIB_TuckMacPhee/Bioinformatics/2018-12_genelists.rda")
# poolGEO <- dbPool(
#   drv = RSQLite::SQLite(),
#   dbname = "/Volumes/ULTRA/across_array/GEOmetadb.sqlite"
# )
#

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

userDatasetTable <- loadUserDatasets(poolUserData)

onStop(function() {
  poolClose(poolGEO)
  poolClose(poolUserData)
})

########################################
#$#$#$#$#$#$#    SERVER    #$#$#$#$#$#$#
########################################

server <- function(input, output, session) {

# Set up colour environment
#_,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,_
  catCol <- brewer.pal(3, "Set1")
  rowCol <-desat(catCol)
  groups <- NULL
  # groups <- c("group1","group2","group3") ## Use these in all following code! They should have a "name" variable for user-assigned names 2018-12-10
  # groups<-c("photoreceptors","RPE","whole.retina") ## what is has to be for the moment
  userID <- NULL
  
# Search functions 
#_,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,_
  ### 2019-03-04 UPDATE to SQL searching directly
  
  searchGSM <- eventReactive(input$searchButton, {
      if(input$gplSelection=='human'){
          sql<-"SELECT * FROM appgsm WHERE description MATCH ?id1 AND gpl LIKE 'GPL570';"
      } else {
          sql<-"SELECT * FROM appgsm WHERE description MATCH ?id1 AND gpl LIKE 'GPL1261';"
      }
      query<-sqlInterpolate(poolGEO,sql,id1=input$searchText)
      queryGSM<-dbGetQuery(poolGEO,query)
      return(queryGSM)
  })

  output$searchResultsGSM <- DT::renderDataTable({
          searchGSM()}, options=list(searching=TRUE, pageLength=50, scrollY='60vh', columnDefs=list(list(
              targets = c(8),
              render = JS(
                  "function(data, type, row, meta) {",
                      "return type === 'display' && typeof data === 'string' && data.length > 100 ?",
                      "'<span title=\"' + data + '\">' + data.substr(0, 100) + '...</span>' : data;",
                      "}") 
                      )))) ## typeof data needs to be a string, as a "NA" converted to JS "NULL" breaks things

# Add sample (array) record to the current experiment 
#_,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,_
  proxy.search = dataTableProxy('searchResultsGSM')
  testTable <- NULL
  gsm_annotated <- eventReactive(input$addButton, {
      testTable <<- rbind(testTable,searchGSM()[input$searchResultsGSM_rows_selected,])
      proxy.search %>% selectRows(NULL)
      return(testTable)
  })
  
  observeEvent(input$addButton, {
      updateTabsetPanel(session = session, inputId = "searchpanel", selected = "2")
  })

# Assign categories to each sample (GSM)
#_,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,_
  ## Set up reactive table to store category data
  samples <- reactiveValues()
  samples$df <- data.frame()
  
  observeEvent(input$assignButton, {
  
      groups <<- c(input$cat1,input$cat2,input$cat3) ## Use these in all following code! They should have a "name" variable for user-assigned names 2018-12-10
 
      if (input$assignButton == 1) {
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
  
  ## ^ don't love this... would like to have the category set without a button click (maybe change to this tab), but it's working for the moment
   
  output$gsm_table <- DT::renderDataTable({
      if(input$assignButton == 0){
         return (datatable(gsm_annotated(),options=list(searching=TRUE, pageLength=50, scrollY='60vh',## 2018-12-10 Pick which columns are necessary ^
             columnDefs=list(list(
             targets = "_all",
             render = JS(
                 "function(data, type, row, meta) {",
                     "return type === 'display' && typeof data === 'string' && data.length > 100 ?",
                     "'<span title=\"' + data + '\">' + data.substr(0, 100) + '...</span>' : data;",
                     "}")
                     )))))
      } else {
         return (datatable(samples$df,options=list(searching=TRUE, pageLength=50, scrollY='60vh',
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
  

  # outputOptions(output, "searchResultsGSM", suspendWhenHidden = FALSE)
  # outputOptions(output, "gsm_table", suspendWhenHidden = FALSE)

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
  finishedtable <- eventReactive(input$assignButton, {
    dplyr::filter(samples$df, category %in% c(input$cat1, input$cat2, input$cat3))
  })
 
  output$finishedtable <- DT::renderDataTable({datatable(finishedtable(),
      options=list(searching=FALSE,pageLength=100, scrollY='60vh')) %>%
      formatStyle('category',target="row",
      backgroundColor=styleEqual(c(input$cat1,input$cat2,input$cat3),c(rowCol[1],rowCol[2],rowCol[3]))
  )})


rv <- reactiveValues(download_flag = 0)

  # proxy.finishedtable = dataTableProxy('finishedtable')
  output$report <- downloadHandler(
      filename = paste(input$downloadId,userID,"GSM_report.csv",sep="_"),
      content = function(file){
          write.csv(finishedtable(),file)
#           tempReport <- file.path(tempdir(),"report.Rmd")
#           file.copy("report.Rmd",tempReport,overwrite=TRUE)
#           params <- list(annotatedGSM = finishedtable())
#
#           rmarkdown::render(tempReport,output_file = file,
#               params = params,
#               envir = new.env(parent=globalenv())
#               )
rv$download_flag <- rv$download_flag + 1
      })
      
observeEvent(input$downloadCEL, {
    
    showModal(modalDialog(title="Important! Downloading raw .CEL files from the NCBI server.","April 11th, 2019: App should be working now. Please click below to begin processing the data.",
    footer = tagList(
        modalButton("Cancel"),
        actionButton("process","Proceed"))))      
  })


  observeEvent(input$process, {
      removeModal()
   })


  observeEvent(input$process, {
      withProgress(
          message = "Downloading and processing GSM",
          {userID<<-processData(finishedtable(),input$downloadId,input$comments,input$gplSelection,poolUserData)})
  })


#   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-
#  / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \
# `-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'

#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$
## This is where the analysis part of the application begins
#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$


# Load dataset
#_,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,_
output$loadUserExperiments = renderUI({
    selectizeInput(inputId="user_data",label="Select an experiment for analysis",choices=c("none"="none",split(userDatasetTable$desc, userDatasetTable$species)),selected="none")
})


observeEvent(input$user_data,{
   if(input$user_data=="none"){
        mapped_probes<<-NULL
        eset<<-NULL
        de_choices<<-NULL
        sig_genes_lfc<<-NULL
        groups <<-NULL
    }else{         
        datasetToLoad <- paste("./data/app_data_",userDatasetTable$userID[which(userDatasetTable$desc == input$user_data)], ".rda", sep='')
        cat(file=stderr(), "attempting to load dataset", datasetToLoad, "based on the id", input$user_data, "\n")
        withProgress(message="Dataset loading",value=0.4,{load(datasetToLoad)})
        mapped_probes<<-mapped_probes
        eset<<-eset
        de_choices<<-de_choices
        sig_genes_lfc<<-sig_genes_lfc
        groups <<- groups
    }
    
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
        genes = c(genes, gene_lists[[gene]])
      }
    }

    if (!is.null(input$gene)) {
      genes = c(genes, input$gene)  
    }
    
    return(unname(genes))
  })
  
  ## gene list UI
  output$geneListsUI = renderUI({
    checkboxGroupInput("genelist", "Select a receptor type to analyze", 
          choices = names(gene_lists))
  })
  
  ## single gene UI
  output$geneUI = renderUI({
    withProgress(message="Loading gene lists",value=0.6,{selectInput("gene", "Select gene(s) to show", choices = all_genes, multiple = TRUE)})
  })
#### This was key to loading the output before we get to this page. All that remains now is either loading both human and mouse, or loading just one depending on the species button. I think loading both at the beginning will help it be snappier overall...
  outputOptions(output, "geneUI", suspendWhenHidden = FALSE)
  
 summary_gene_data = reactive({
   validate(
      need(geneList(), "No genes selected"),
      need(!is.null(eset),"No dataset selected")
    )
   get_expression_summary(eset, geneList())
 })

# QC output
#_,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,_
  
 output$QC = renderUI({
    validate(
      need(input$user_data!="none","No dataset selected")
    )
    fluidRow(h4("Expression normalization (array intensity, before and after)"), tags$img(src="array-processing.png",width="100%"), h4("RNA degradation plot (probe position vs intensity)"),
    tags$img(src="RNA-deg.png",width="100%"))
  })
  
  
  output$genes = DT::renderDataTable({
    validate(
      need(geneList(), "No genes selected")
    )
    
     summary_gene_data() %>% datatable() %>% 
      formatRound(2:4)
    
  })
  
  # single gene plot
 output$singleGenePlot = renderPlot({
  # output$singleGenePlot = renderTable({
    validate(
      need(input$genes_rows_selected >= 1, "No genes selected")
    )
    
    rows = as.integer(input$genes_rows_selected)
    genes_to_plot = summary_gene_data()$Symbol[rows]
    
    gene_data = get_gene_data(eset, genes_to_plot)
    by_gene_violplot(gene_data,tissues=groups)
    
    
  })

  # DE choices UI
  output$de_choices = renderUI({
    checkboxGroupInput("de", "Choose comparison(s) to show", choices = de_choices, selected = de_choices[1])
  })

# Expression tab
#_,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,__,.-'~'-.,_
  observe({
    toggle("de_choices", anim = TRUE, condition = input$de_state )
  })
  
  genesToPlot = reactive({
    validate(
      need(geneList(), "No genes selected")
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
      need(genesToPlot(), "No genes selected"),
      need(input$tissues, "No tissues selected")
    )
   
    selected_tissues = input$tissues
    sub_eset = eset[, eset$tissue %in% selected_tissues]
    genes = gene2probe(genesToPlot(), mapped_probes)
    
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
      need(genesToPlot(), "No genes selected"),
      need(input$tissues, "No tissues selected")
    )
    
    gene_data = get_gene_data(eset, genesToPlot())
    overall_expression_boxplot(gene_data, tissues = input$tissues)
    
  })
  

# By gene boxplots ----------------------------------------------------------------------------

  output$byGenePlot = renderPlot({
    validate(
      need(genesToPlot(), "No genes selected"),
      need(input$tissues, "No tissues selected")
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
    
    probe = input$pls_probe
    #ncomp = input$pls_ncomp
    
    get_plsda(sub_eset, genes, probe) 
    
  })

# PCA plot ----------------------------------------------------------------------------
  output$indPlot = renderPlot({
    validate(
      need(plsdaData(), "No PLS-DA to plot"),
      need(length(input$pls_tissues) >= 2, "Please select at least two tissues")
    )
    
    plotIndiv(plsdaData()$result, ind.names = FALSE, group = plsdaData()$tissue_grps, pch = 16, 
              col.per.group = brewer.pal(3, "Set1")[1:length(input$pls_tissues)], legend = TRUE, cex = 2, ellipse=TRUE)
  })

# Correlation Circle plot ----------------------------------------------------------------------------  
  output$varPlot = renderPlot({
     validate(
      need(plsdaData(), "No PLS-DA to plot")
    )

    plotVar(plsdaData()$result, var.names = list(plsdaData()$varNames), cex = 3,overlap=FALSE)
    
  })

  output$numGenesUI = renderUI({
    numericInput("pls_num_genes", "Select number of genes to show contributions for", 
                 value = 10, min = 1, max = length(geneList()), step = 1)
  })
  
# Loadings plot ----------------------------------------------------------------------------
  output$contribPlot = renderPlot({
    validate(
      need(plsdaData(), "No PLS-DA to plot"),
      need(input$pls_num_genes, "")
    )
    
    grps = plsdaData()$result$names$Y
    cols = brewer.pal(3, "Set1")[1:length(grps)]
     
    ndisplay = input$pls_num_genes
    comp = as.integer(input$pls_ncomp)
    plotLoadings(plsdaData()$result, name.var = plsdaData()$varNames, ndisplay = ndisplay,
                comp = comp, legend.color = c(1:2))
     
  })
  
#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$  
  ## Kill shinyApp when session closes
  session$onSessionEnded(stopApp)

}
