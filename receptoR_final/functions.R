#                          _        ____
#  _ __ ___  ___ ___ _ __ | |_ ___ |  _ \
# | '__/ _ \/ __/ _ \ '_ \| __/ _ \| |_) |
# | | |  __/ (_|  __/ |_) | || (_) |  _ <
# |_|  \___|\___\___| .__/ \__\___/|_| \_\
#                   |_|
#
# June 2019 receptoR v 1.3
## Last update: 2019-06-22, Derek Toms
## functions.R


#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$
desat = function(cols, sat=0.5) {
    X <- diag(c(1, sat, 1)) %*% rgb2hsv(col2rgb(cols))
    hsv(X[1,], X[2,], X[3,])
}


#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$
 #$# Data processing 9 April 2018
 ##  Updates 2018-12-10 to switch between user input and CEL downloads
 
processData = function(finished_table,datasetID,userComments,gpl,userDB){
    ## categories need to have R-safe names so we'll store the user input as a new column
    finished_table$category.labels <- finished_table$category
    safeNames <- make.names(levels(factor(finished_table$category)),unique=TRUE)
    finished_table$category <- factor(finished_table$category, levels = levels(factor(finished_table$category)), labels=safeNames)
    ## colours
    catCol <- factor(brewer.pal(9,"Set1")[1:length(levels(factor(finished_table$category)))])
    finished_table$colours <- finished_table$category
    finished_table$colours <- factor(finished_table$colours, levels=levels(factor(finished_table$category)), labels=catCol)
    
    ## timestamp
    timeStamp <- strftime(Sys.time(),"%Y%m%d-%H%M")

    PATH = "./data/"
    get_files = TRUE
    # get_files = FALSE

    gsm_to_fetch <- finished_table$gsm
    
    if (get_files) {
        # get raw CEL files

        setDir<-paste(PATH,timeStamp,sep='')
          
withProgress(
   message = "Downloading and processing GSM", value = 0,
   {    
# Download expression files ------------------------------------------------------------
      rawFilePaths = lapply(gsm_to_fetch, function(x) {
        incProgress(0.5/length(gsm_to_fetch), message = "Downloading expression data")
        dir.create(file.path(setDir), showWarnings = FALSE)
        getGEOSuppFiles(x,baseDir=setDir)    
      })
        
        gsm_dirs = list.files(path=setDir, pattern = "GSM", full.names=TRUE)
        gsm_files = lapply(gsm_dirs, list.files, pattern = "[Cc][Ee][Ll].gz", full.names = TRUE)

# Normalize together ------------------------------------------------------------
        incProgress(0.1, message = "Performing normalization")
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

        ## QC – save PNG files for degradation and normalization
        ### save pictures of processed data
        
        finished_table <- finished_table[order(match(all_pData$tissue,finished_table$category)),] # order the two tables the same way
                        
        ## Array normalization
        mat<-matrix(c(1,2,3),3)
        png(filename = paste("./www/array_normalization_", timeStamp, ".png", sep=''), width=2400, height=1800, units="px")
        layout(mat,widths=c(1,1,1),heights=c(1,2,3))
        par(mar=c(1,3,1,1))
                ## Experimental samples (NULL plot, legend only)
        plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
        legend("topleft", legend =levels(factor(finished_table$category.labels)), pch=15, pt.cex=3, cex=1.5, bty='n', col = levels(finished_table$colours))
        par(mar=c(1,3,1,1))
        boxplot(all_data,las=2,main="Raw expression data",xaxt="n",col=paste(finished_table$colours))
        par(mar=c(7,3,1,1))
        boxplot(exprs(all_eset),las=2,main="Normalized expression data",col=paste(finished_table$colours))
        dev.off()
    
        ## Probe degradation
        all_data_deg <- AffyRNAdeg(all_data)
        ggsave(file = paste("./www/probe_degradation_", timeStamp, ".png", sep=''), plot =  plotAffyRNAdeg(all_data_deg, cols=paste(finished_table$colours)))
        
        ## Change of object name
        all_eset_final<-all_eset
        pData(all_eset_final)<-all_pData

        ## Save data files 
        save(all_eset_final, all_data, finished_table, file = paste(PATH, "final_processed_data_", timeStamp, ".rda", sep=''))

# Differentially Expressed Gene (DEG) Analysis ------------------------------------------------------------
        incProgress(0.1, message = "Determining differential gene expression")

        ## Set gene symbols based on species
        ID = featureNames(all_eset_final)
        if(gpl =='human'){
            Symbol = getSYMBOL(ID, "hgu133plus2.db")
            mapped_probes = as.list(revmap(hgu133plus2ALIAS2PROBE))
            # reverse map for symbol to probe conversion

        } else {
            Symbol = getSYMBOL(ID, "mouse4302.db")
            mapped_probes = as.list(revmap(mouse4302SYMBOL))
        }


        fData(all_eset_final) = data.frame(Symbol = Symbol)

        eset = all_eset_final

        tissue = as.factor(pData(eset)$tissue)
        design = model.matrix(~0 + tissue)
        colnames(design) = levels(tissue)
        
        incProgress(0.1, message = "Fitting model")
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

        # get list of DEG
        de_choices = names(sig_genes_lfc)
        # set groups
        groups = levels(tissue)

# Save user-generated experiments -----------------------------------   
        incProgress(0.1, message = "Saving processed data")
        db <- poolCheckout(userDB)
        data <- data.frame(userID = timeStamp, desc = datasetID, comments = userComments, species = gpl)
        dbWriteTable(conn=db, name="userData", data, append=T, row.names=F)
        poolReturn(db)

        save(mapped_probes, eset, de_choices, sig_genes_lfc, groups, file = paste(PATH,"app_data_",timeStamp,".rda",sep=''))
        
        incProgress(0.1, message = "Success!")
}) # end withProgress
    
        return(timeStamp)
    
    } else {
        save(userComments,finished_table, file = paste("annotated_gsm_",timeStamp,".rda",sep=''))
        return(timeStamp)
    }

}

loadUserDatasets <- function(userDB) {
     # Connect to the database
     db <- poolCheckout(userDB)
     # Construct the fetching query
     query <- sprintf("SELECT * FROM userData")
     # Submit the fetch query and disconnect
     data <- dbGetQuery(db, query)
     poolReturn(db)
     data
 }

 # end processing
 #$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$
 

 ## Functions for processing array expression data ##
 
 
#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$

 # Convert gene symbols back to probes --------------------------------------------------------

 gene2probe = function(gene_list, mapped_probes) {

     if(species=='mouse'){
         return(na.omit(unlist(mapped_probes[gene_list])))
     } else {
         return(names(mapped_probes)[sapply(mapped_probes, function(pr) any(gene_list %in% pr))])
     }
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
    spread(tissue, expression)
  
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

   cat(file=stderr(), "incoming ", length(subset_probes), " probes\n")
   #mat = t(scale(t(mat)))
   if(species == 'mouse'){
       mat = exprs(eset)[c(subset_probes),]
       cat(file=stderr(), "attempting to load ", length(rownames(mat)), " mouse genes\n")
       row_labs = paste(getSYMBOL(rownames(mat), "mouse4302.db"),rownames(mat),sep=":")
       rownames(mat) = getSYMBOL(rownames(mat), "mouse4302.db")
   } else {
       mat = exprs(eset)[c(subset_probes),]
       cat(file=stderr(), "attempting to load ", length(rownames(mat)), " human genes\n")
       row_labs = paste(getSYMBOL(rownames(mat), "hgu133plus2.db"),rownames(mat),sep=":")
       rownames(mat) = getSYMBOL(rownames(mat), "hgu133plus2.db")
       row_labs = row_labs[!is.na(rownames(mat))]
       mat = mat[which(!is.na(rownames(mat))),]
       
   }
   
   # debug
    cat(file=stderr(), "These are the genes to plot: ", rownames(mat),"\n")
   
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
   if(species=='mouse'){
       exprs(eset) %>% 
         as.data.frame() %>% 
         tibble::rownames_to_column("probe") %>% 
         mutate(Symbol = getSYMBOL(probe, "mouse4302.db")) %>% 
         filter(Symbol %in% gene_list) %>% 
         gather(Sample, expression, starts_with("GSM")) %>% 
         left_join(ph, by = "Sample")
   } else {
     exprs(eset) %>% 
       as.data.frame() %>% 
       tibble::rownames_to_column("probe") %>% 
       mutate(Symbol = getSYMBOL(probe, "hgu133plus2.db")) %>% 
       filter(Symbol %in% gene_list) %>% 
       gather(Sample, expression, starts_with("GSM")) %>% 
       left_join(ph, by = "Sample")
   }
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
     
     if(species=='mouse'){
         exp = exprs(eset)[c(genes),]
     } else {
         exp = exprs(eset)[genes,]
     }
     
     tissue = factor(pData(eset)$tissue)
     tissue_grps = pData(eset)$tissue
    
     if (probe == FALSE) {
       exp = aggregate(exp, list(genes = rownames(exp)), mean)
       rownames(exp) = exp$genes
       exp = as.matrix(exp[,-1])
     }
    
     exp = t(exp)
     if(species=='mouse'){
         return(list(
           result = splsda(exp, tissue, ncomp = 2),
           tissue_grps = tissue_grps,
           varNames = getSYMBOL(colnames(exp), "mouse4302.db")
         ))
     } else {
         return(list(
           result = splsda(exp, tissue, ncomp = 2),
           tissue_grps = tissue_grps,
           varNames = getSYMBOL(colnames(exp), "hgu133plus2.db")
         ))
     }
 } 