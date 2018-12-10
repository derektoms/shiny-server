# October 2018 receptoR v 1.0
## 2018-11-23 
### Developing data management for multiple users and polishing the app
## 2018-10-27
## functions.R
### Integrating both applications to a final shiny executable

#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$
desat = function(cols, sat=0.5) {
    X <- diag(c(1, sat, 1)) %*% rgb2hsv(col2rgb(cols))
    hsv(X[1,], X[2,], X[3,])
}


#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$
 #$# Data processing 9 April 2018
 ##  Updates 2018-12-10 to switch between user input and CEL downloads
 
 processData = function(finished_table){
 gsm_to_fetch <- finished_table$gsm
 ## timestamp
 timeStamp <- strftime(Sys.time(),"%Y%m%d-%H%M")

 # get_files = TRUE
 get_files = FALSE
 
 if (get_files) {
   # get raw CEL files

   setDir<-paste(getwd(),'/',timeStamp,sep='')

   rawFilePaths = lapply(gsm_to_fetch, function(x) {
       dir.create(file.path(setDir), showWarnings = FALSE)
       getGEOSuppFiles(x,baseDir=setDir)
   })
 
gsm_dirs = list.files(path=setDir, pattern = "GSM", full.names=TRUE)

gsm_files = lapply(gsm_dirs, list.files, pattern = "[Cc][Ee][Ll].gz", full.names = TRUE)

 # New approach: Normalize together ------------------------------------------------------------
 all_data = ReadAffy(filenames = unlist(gsm_files))
 all_eset = rma(all_data)
 all_pData = pData(all_eset)

 gsm = str_match(rownames(all_pData), "(GSM[[:alnum:]]+)")[,2]

 ## missing CEL files (2018-04-13) required the fix below
 gsm <- tibble(gsm)
 dlSamples <- inner_join(gsm,finished_table)

 all_pData<-tibble("tissue"=dlSamples$category)
 colnames(all_eset)<-dlSamples$gsm
 rownames(all_pData)<-dlSamples$gsm # Warning: setting row names on a tibble is deprecated
 ## end fix

 all_eset_final<-all_eset
 pData(all_eset_final)<-all_pData
 pData(all_eset_final) %>% View

 identical(colnames(exprs(all_eset_final)), rownames(pData(all_eset_final)))

 # run_tsne(t(exprs(all_eset_final)), pData(all_eset_final))
 ##^ couldn't do this


 # ok that looks good let's save this now
 ## 2018-04-13
 ## just before saving, there is a problem with the annotation (see below at the call to contrast_matrix)
 ## quick and dirty fix:

 # pData(all_eset_final)<-pData(all_eset_final)%>%mutate(tissue=str_replace(tissue,"whole retina","whole.retina"))

 save(all_eset_final, file = paste("final_processed_data_",timeStamp,".rda",sep='')) # filename should include timestamp

 ### DE

 ID = featureNames(all_eset_final)
 # Symbol = getSYMBOL(ID, "hgu133plus2.db")
 Symbol = getSYMBOL(ID, "mouse4302.db")
 fData(all_eset_final) = data.frame(Symbol = Symbol)

 eset = all_eset_final

 tissue = as.factor(pData(eset)$tissue)
 design = model.matrix(~0 + tissue)
 colnames(design) = levels(tissue)
 
 fit = lmFit(eset, design)
 matrices<-t(combn(levels(tissue),2))
 contrasts <- paste(matrices[,1],matrices[,2],sep='-')
 
 ## 2018-04-13
 ## Ran into more issues at the contrast matrix because the levels need to by syntatically allowed (i.e. no spaces)
 ## adding code higher up to avoid this
 contrast_matrix = makeContrasts(contrasts=contrasts, levels = design)
 
 fit2 = contrasts.fit(fit, contrast_matrix)
 efit = eBayes(fit2)
 tfit = treat(fit2, lfc = 1)

 results = decideTests(efit)
 results_lfc = decideTests(tfit)

 # I like these... they should return something in the UI!
 vennDiagram(results, include = c("up", "down"))
 vennDiagram(results_lfc)
 
  coefs = colnames(contrast_matrix)
  sig_genes = lapply(coefs, function(x) {
    topTable(efit, coef = x, number = Inf, p.value = 0.01, sort.by = "p")
  })

  sig_genes_lfc = lapply(coefs, function(x) {
    topTreat(tfit, coef = x, number = Inf, p.value = 0.01, sort.by = "p")
  })
  names(sig_genes) = coefs
  names(sig_genes_lfc) = coefs

  sapply(sig_genes, nrow)
  sapply(sig_genes_lfc, nrow)
  
   # reverse map for symbol to probe conversion
   mapped_probes = as.list(revmap(mouse4302SYMBOL))

   de_choices = names(sig_genes_lfc)

   all_genes = featureData(eset)@data[["Symbol"]] %>% as.character() %>% unique()

 save(mapped_probes, eset, de_choices, sig_genes_lfc, file = paste("app_data_",timeStamp,".rda",sep='')) # filename should include timestamp
   # save(all_genes, mapped_probes, eset, de_choices, sig_genes_lfc, file = "2018-04-13_app_data.rda")
   #  save(all_genes, gene_lists, file = " 2018-12_genelists.rda")
   return(timeStamp)
} else {
    save(finished_table, file = paste("annotated_gsm_",timeStamp,".rda",sep=''))
    return(timeStamp)
}

}

 # end processing
 #$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$
 

 ## Functions for processing array expression data ##
 
 
#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$

 # Convert gene symbols back to probes --------------------------------------------------------

 gene2probe = function(gene_list, mapped_probes) {
   na.omit(unlist(mapped_probes[gene_list]))
 }


 # Get de genes from topTable output that match a genelist -------------------------------------

 get_de_genes = function(gene_list, de_choice, sig_genes_lfc) {
   sig_genes_lfc[[de_choice]] %>% tibble::rownames_to_column("probe") %>% filter(Symbol %in% gene_list)
 }


 # Get summary expression table ----------------------------------------------------------------

 get_expression_summary = function(eset, gene_list) {
   
   get_gene_data(eset, gene_list) %>% 
    group_by(Symbol, tissue) %>% 
    summarise(expression = mean(expression)) %>% 
    spread(tissue, expression) # %>% rowwise() %>%
    # mutate(mean = mean(tissue2,tissue1,tissue3)) %>%
    # ^ this code was never right; only photoreceptor mean was calculated
    # arrange(desc(mean))
  
 }

 #' Gene Heatmap
 #' 
 #' Plot a heatmap with provided expression set and gene list
 #' @param eset Expression set
 #' @param subset_probes subset of probes from eset to plot
 #' @param anno_col Annotation column, must be a column of pData(eset)
 #' @param probe_level Plot expression of probes instead of aggregating to genes
 #' @param ... additional arguments to pheatmap
 #' @export
 gene_heatmap = function(eset, subset_probes, anno_col = 'tissue', probe_level = FALSE, gsm_show = TRUE, ...) {

   mat = exprs(eset)[subset_probes,]
   #mat = t(scale(t(mat)))
   row_labs = paste(getSYMBOL(rownames(mat), "mouse4302.db"),rownames(mat),sep=":")
   rownames(mat) = getSYMBOL(rownames(mat), "mouse4302.db")
   mat = mat[order(rownames(mat)),]
   
   if (!probe_level) {
     mat = aggregate(mat, list(genes = rownames(mat)), mean)
     rownames(mat) = mat$genes
     mat = as.matrix(mat[,-1])
     row_labs=rownames(mat)
   }
  
   anno_df = pData(eset)[anno_col]
   # another fix for Tibble deprecating row names
   anno_df <- data.frame(anno_df)
   rownames(anno_df) <- colnames(mat)
   # end fix
   vars = as.character(unique(anno_df[[anno_col]]))
   anno_colours = brewer.pal(9, "Set1")[1:length(vars)]
   names(anno_colours) = vars
   anno_colours = list(anno_colours)
   names(anno_colours) = anno_col 
  
   pheatmap(mat, show_colnames = gsm_show, annotation_col = anno_df, annotation_colors = anno_colours, labels_row = row_labs,...)
  
 }


 # Get expression data from eset for given genes -----------------------------------------------

 get_gene_data = function(eset, gene_list) {
   ph = pData(eset) %>% tibble::rownames_to_column("Sample")
   ph$Sample <- colnames(exprs(eset)) # 2018-04-17 fix missing row names
   exprs(eset) %>% 
     as.data.frame() %>% 
     tibble::rownames_to_column("probe") %>% 
     mutate(Symbol = getSYMBOL(probe, "mouse4302.db")) %>% 
     filter(Symbol %in% gene_list) %>% 
     gather(Sample, expression, starts_with("GSM")) %>% 
     left_join(ph, by = "Sample")
 }

 # Create a top N gene list by abs LFC ---------------------------------------------------------

 top_genes = function(de_data, top = 25) {
   de_data %>% 
     group_by(Symbol) %>% 
     summarise(meanLFC = mean(logFC)) %>% 
     top_n(top, abs(meanLFC))
 }


 # Make a boxplot faceted by gene --------------------------------------------------------------

 by_gene_boxplot = function(gene_data, tissues = groups) {
   gene_data %>% 
     filter(tissue %in% tissues) %>%
     ggplot(aes(x = tissue, y = expression)) +
       geom_boxplot(aes(fill = tissue)) +
       facet_wrap(~Symbol) +
       theme_bw() + theme(axis.text.x = element_blank()) +
       scale_fill_brewer(palette = "Set1")
 }
 
 by_gene_violplot = function(gene_data, tissues = groups) {
    gene_data %>% 
      filter(tissue %in% tissues) %>%
      ggplot(aes(x = tissue, y = expression)) +
        geom_violin(aes(fill = tissue)) +
        facet_wrap(~Symbol) +
        theme_bw() + theme(axis.text.x = element_blank()) +
        scale_fill_brewer(palette = "Set1")
  }

 # Make a box plot for the expression of given genes, on plot per tissue -----------------------

 overall_expression_boxplot = function(gene_data, tissues = groups ) {

   tissue_cols = brewer.pal(3, "Set1")[1:length(tissues)]
   names(tissue_cols) = tissues

   gene_plots = gene_data %>% 
     filter(tissue %in% tissues) %>% 
     group_by(tissue, Symbol) %>% mutate(meanExp = mean(expression)) %>% 
     group_by(tissue) %>% 
     do(plot = mutate(., Symbol = factor(Symbol, levels = unique(Symbol[order(meanExp, decreasing = TRUE)]))) %>% 
       ggplot(aes(x = Symbol, y = expression)) +
         geom_boxplot(fill = tissue_cols[unique(.$tissue)]) +
         theme_bw() +
         theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
         labs(title = unique(.$tissue), x = "")
     )
   plot_grid(plotlist = gene_plots$plot, ncol = 1)    
    
 }


 # Perform sPLS-DA on selected genes -----------------------------------------------------------

 get_plsda = function(eset, genes, probe) {
    
     exp = exprs(eset)[genes,]
     tissue = factor(pData(eset)$tissue)
     tissue_grps = pData(eset)$tissue
    
     if (probe == FALSE) {
       exp = aggregate(exp, list(genes = rownames(exp)), mean)
       rownames(exp) = exp$genes
       exp = as.matrix(exp[,-1])
     }
    
     exp = t(exp)
     return(list(
       result = splsda(exp, tissue, ncomp = 2),
       tissue_grps = tissue_grps,
       varNames = getSYMBOL(colnames(exp), "mouse4302.db")
     ))
 } 