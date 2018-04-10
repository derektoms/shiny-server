
## Functions for photoreceptor app ##


# Converte gene symbols back to probes --------------------------------------------------------

gene2probe = function(gene_list, mapped_probes) {
  na.omit(unlist(mapped_probes[gene_list]))
}


# Get de genes from topTable output that match a genelist -------------------------------------

get_de_genes = function(gene_list, de_choice, sig_genes_lfc) {
  sig_genes_lfc[[de_choice]] %>% add_rownames("probe") %>% filter(Symbol %in% gene_list)
}


# Get summary expression table ----------------------------------------------------------------

get_expression_summary = function(eset, gene_list) {
  
  get_gene_data(eset, gene_list) %>% 
    group_by(Symbol, tissue) %>% 
    summarise(expression = mean(expression)) %>% 
    spread(tissue, expression) %>% rowwise() %>%  
    mutate(mean = mean(photoreceptors,retina,RPC)) %>% 
    arrange(desc(mean))
  
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
gene_heatmap = function(eset, subset_probes, anno_col = 'tissue', probe_level = FALSE, ...) {

  mat = exprs(eset)[subset_probes,]
  #mat = t(scale(t(mat)))
  rownames(mat) = getSYMBOL(rownames(mat), "mouse4302.db") 
  mat = mat[order(rownames(mat)),]
  
  if (!probe_level) {
    mat = aggregate(mat, list(genes = rownames(mat)), mean)
    rownames(mat) = mat$genes
    mat = as.matrix(mat[,-1])
  }
  
  anno_df = pData(eset)[anno_col]
  vars = as.character(unique(anno_df[[anno_col]]))
  anno_colours = brewer.pal(9, "Set1")[1:length(vars)]
  names(anno_colours) = vars
  anno_colours = list(anno_colours)
  names(anno_colours) = anno_col 
  
  pheatmap(mat, show_colnames = FALSE, annotation_col = anno_df, annotation_colors = anno_colours, ...)
  
}


# Get expression data from eset for given genes -----------------------------------------------

get_gene_data = function(eset, gene_list) {
  
  ph = pData(eset) %>% add_rownames("Sample")
  exprs(eset) %>% 
    as.data.frame() %>% 
    add_rownames("probe") %>% 
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

by_gene_boxplot = function(gene_data, tissues = c("retina", "photoreceptors", "RPC")) {
  gene_data %>% 
    filter(tissue %in% tissues) %>%
    ggplot(aes(x = tissue, y = expression)) +
      geom_boxplot(aes(fill = tissue)) +
      facet_wrap(~Symbol) +
      theme_bw() + theme(axis.text.x = element_blank()) +
      scale_fill_brewer(palette = "Set1")
}


# Make a box plot for the expression of given genes, on plot per tissue -----------------------

overall_expression_boxplot = function(gene_data, tissues = c("retina", "photoreceptors", "RPC") ) {

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
