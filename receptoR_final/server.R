# October 2018 receptoR v 1.0
## Last update: 2018-11-23
## server.R
### Integrating both applications to a final shiny executable

#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$

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
library(tidyr)
library(ggplot2)
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
# library(hgu133plus2.db)

########################################
#$#$#$#$#$#$    Shiny App  $#$#$#$#$#$#$
########################################

source("functions.R")

### not the best place to put this, but it should work for now

# load("../../data/gseGPL570.rda")
# load("../../data/gsmGPL570.rda")
load("../gseGPL1261.rda")
load("../gsmGPL1261.rda")

load("../2018-12_genelists.rda")
### I'm going to try and not have this loaded to start
# load("2018-04-13_app_data.rda")

## SERVER
server <- function(input, output, session) {
  catCol <- brewer.pal(3, "Set1")
  rowCol <-desat(catCol)
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
          filtered_gse()}, options=list(searching=TRUE, pageLength=50, columnDefs=list(list(
              targets = c(8,9,12),
              render = JS(
                  "function(data, type, row, meta) {",
                      "return type === 'display' && typeof data === 'string' && data.length > 100 ?",
                      "'<span title=\"' + data + '\">' + data.substr(0, 100) + '...</span>' : data;",
                      "}") 
                      )))) ## typeof data needs to be a string, as a "NA" converted to JS "NULL" breaks things

#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$
### 2018-10-28 Disable platform selection to get it working with mice
shinyjs::disable("gplSelection")
shinyjs::disable("downloadCEL")

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
          return (datatable(gsm_annotated()[,c(-5:-7,-11,-12,-14:-26,-28:-32)],options=list(searching=FALSE, pageLength=50,
              columnDefs=list(list(
              targets = "_all",
              render = JS(
                  "function(data, type, row, meta) {",
                      "return type === 'display' && typeof data === 'string' && data.length > 100 ?",
                      "'<span title=\"' + data + '\">' + data.substr(0, 100) + '...</span>' : data;",
                      "}") 
                      )))))
       } else {
          return (datatable(samples$df[,c(-5:-7,-11,-12,-14:-26,-28:-32)],options=list(searching=FALSE, pageLength=50,
              columnDefs=list(list(
              targets = "_all",
              render = JS(
                  "function(data, type, row, meta) {",
                      "return type === 'display' && typeof data === 'string' && data.length > 100 ?",
                      "'<span title=\"' + data + '\">' + data.substr(0, 100) + '...</span>' : data;",
                      "}") 
                      ))))%>% formatStyle('category',target="row",backgroundColor=styleEqual(c(input$cat1,input$cat2,input$cat3),c(rowCol[1],rowCol[2],rowCol[3]))))
       }
  })
                      
  proxy.gsm = dataTableProxy('gsm_table')
  observeEvent(input$Assign,{
      proxy.gsm %>% selectRows(NULL)
  }) 
  
  ## UI output

    output$categorySelect <- renderUI(
      fluidRow(
        column(12,
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

  finishedtable <- eventReactive(input$Assign, {
    dplyr::filter(samples$df, category %in% c(input$cat1, input$cat2, input$cat3))
  })
 
  output$finishedtable <- DT::renderDataTable({datatable(finishedtable()[,c(2,3,4,10,31,32,33)],, options=list(pageLength=100, scrollY=220)) %>% formatStyle('category',target="row",backgroundColor=styleEqual(c(input$cat1,input$cat2,input$cat3),c(rowCol[1],rowCol[2],rowCol[3])))})
    
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
 
  formData <- eventReactive(input$Assign, {
      time <- strftime(Sys.time(),"%Y.%m.%d %H:%M")
      stamp <- data.frame("timestamp"=rep(time,nrow(finishedtable())))
      CELdl <- c(finishedtable()[,3],finishedtable()[,33],stamp)
      CELdl
  })
 
 output$CELdl <- renderTable(formData()) 
 
    groups<-c("photoreceptors","RPE","whole.retina")

#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$
## This is where the analysis part of the application begins
#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$



observeEvent(input$user_data,{
   if(input$user_data=="none"){
        mapped_probes<<-NULL
        eset<<-NULL
        de_choices<<-NULL
        sig_genes_lfc<<-NULL
    }else{
        withProgress(message="Dataset loading",value=0.4,{load("../2018-04-13_app_data.rda",envir=.GlobalEnv)})
    }
    
})

# Load genes tab ------------------------------------------------------------------------------

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
  
  # gene list UI
  output$geneListsUI = renderUI({
    checkboxGroupInput("genelist", "Select a gene list", 
          choices = names(gene_lists))
  })
  
  # single gene UI
  output$geneUI = renderUI({
    withProgress(message="Dataset loading",value=0.4,{selectInput("gene", "Select gene(s) to show", choices = all_genes, multiple = TRUE)})
  })
  
 summary_gene_data = reactive({
   validate(
      need(geneList(), "No genes selected"),
      need(!is.null(eset),"No dataset selected")
    )
   get_expression_summary(eset, geneList())
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
    by_gene_violplot(gene_data,tissues=c("photoreceptors","RPE","whole.retina"))
    
    
  })

  # DE choices UI
  output$de_choices = renderUI({
    checkboxGroupInput("de", "Choose comparison(s) to show", choices = de_choices, selected = de_choices[1])
  })

# Expression tab ------------------------------------------------------------------------------
  
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
  
  
# heatmap plot --------------------------------------------------------------------------------
  
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

# Overall expression --------------------------------------------------------------------------

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
