## temporary packages from the first part, necessary for this to work
#####---------------------
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
# source("https://bioconductor.org/biocLite.R")
# biocLite(c('limma','annotate','genefilter','ComplexHeatmap','pheatmap','cowplot','GEOmetadb','mouse4302.db','hgu133plus2.db'))

library(GEOmetadb)
library(GEOquery)

library(affy)

library(limma)
library(annotate)
library(genefilter)
library(ComplexHeatmap)
library(pheatmap)
library(cowplot)

#biocLite(c('mixOmics','MergeMaid','GEOquery','inSilicoMerging','affy','sva','Rtsne','metaArray','testthat'))

library(mixOmics)
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
#####---------------------


## receptoRA server.R
source("../receptoR_final/functions.R")
load("../genelists.rda")
load("../2018-04-13_app_data.rda")

server<-function(input, output, session) {
  
    ## load up the groups (this information should be saved, but isn't)
    groups<-c("photoreceptors","RPE","whole.retina")
    
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
    selectInput("gene", "Select gene(s) to show", choices = all_genes, multiple = TRUE)
  })
  
 summary_gene_data = reactive({
   validate(
      need(geneList(), "No genes selected")
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
    by_gene_boxplot(gene_data,tissues=c("photoreceptors","RPE","whole.retina"))
    
    
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
  
  # contribData = reactive({
  #
  #   ndisplay = input$pls_num_genes
  #   comp = as.integer(input$pls_ncomp)
  #   contrib = plotLoadings(plsdaData()$result, name.var = plsdaData()$varNames, ndisplay = ndisplay,
  #               comp = comp, plot = TRUE)
  #
  #   contrib$contrib %>% dplyr::select(-Contrib) %>% add_rownames("Gene") %>%
  #     mutate(Gene = getSYMBOL(Gene, "mouse4302.db"))
  # })
  
  # output$contribTable = renderDataTable({
#     validate(
#       need(contribData(), "No data")
#     )
#      contribData()
#   })
  
  # output$pls_download = downloadHandler(
#     filename = 'gene_contribution_data.csv',
#     content = function(file) {
#       write_csv(contribData(), file)
#     }
#
#   )
    session$onSessionEnded(stopApp)
}

  
