#' Perform GSEA on Differential Expression Output
#'
#' \code{cluster_wrapper} will perform GSEA quickly and efficiently on your differential expression output. All you need is a dataframe of genes and their associated log fold change. Future updates will allow for multi-level GSEA (eg: cluster marker analysis for scRNA)
#'
#' @param data An object of either class data.frame or a named numeric vector in one one of two formats. If data is a data.frame, you must supply gene_col and stat_col (see below). If it is a named vector, the names should be the gene names and the vector should be the stat value
#'
#' @param gene_col A character string of length one that is the name of the column in data that contains gene symbols. If NULL, will attempt to default to column 1
#'
#' @param stat_col A character string of length one that is the name of the column in data that contains the numeric values for ranking, can be logFC, pval*sign of log, etc. If NULL, will attempt to default to column 2
#'
#' @param species character string that is the species of interest for GSEA (eg: Mus musculus, Homo sapiens)
#'
#' @param category a character string that fits with msigdbr::, tells the program what to use
#'
#' @param n_max a numeric of length one that tells the max number of GSEA plots to return (does not effect the dataframe). Defaults to 10
#'
#' @param inputType a character string of length one that tells the naming convention that is used. Default is SYMBOL, acceptable inputs are also ENTREZID and ENSEMBL
#'
#'
#' @return Returns an object of class 'list' that contains all GSEA plots of significantly enriched sets based on the given results file. At the moment, this is limited to A vs B comparisons.
#'
#'
#' @examples
#' cluster_wrapper(data = df, gene_col = "SYMBOL",stat_col = "avg_logFC")
#' cluster_wrapper(data = df, gene_col = "SYMBOL",stat_col = "avg_logFC", category = "C2")
#' cluster_wrapper(data = df, species = "Human",gene_col = "SYMBOL",stat_col = "avg_logFC", category = "H")
#'

#' @export
#' @import org.Mm.eg.db
#' @import org.Hs.eg.db
#' @import clusterProfiler
#' @import dplyr

cluster_wrapper <- function(data, gene_col = NULL, stat_col = NULL,
                            species = "Homo sapiens",
                            category = NULL, subcategory = NULL,
                            inputType = "SYMBOL", exponent = 1, maxGSSize = 2000,
                            n_max = 10, pval_cutoff = .99, nPerm = 50000){
  
  output <- list()
  ## First, pick your org.db
  if(species %in% c("Mus musculus","mouse","Mouse")){
    db = org.Mm.eg.db
  }
  else if(species %in% c("Homo sapiens","human","Human")){
    db = org.Hs.eg.db
  }
  else{
    warning("no recognized Org.db for input species, see documentation")
  }
  # ========================================================================
  # ========================================================================
  
  # See if input is data.frame or named number
  if(class(data) == "data.frame"){
    # Trim dataset to desired columns
    if(is.null(gene_col) | is.null(stat_col)){
      warning("Data.frame input require both gene_col and stat_col be provided, please see documentation for help")
    }
    else{
      data = data[,c(gene_col,stat_col)]
      colnames(data) <- c("gene","stat")
      geneList = as.numeric(data[,"stat"])
      # If not entrez id already, convert here
      if(inputType != "ENTREZID"){
        names(geneList) = clusterProfiler::bitr(data$gene,
                                                fromType = inputType, toType = "ENTREZID",
                                                OrgDb = db)$ENTREZID
      }
      else{
        names(geneList) = data$gene
      }
      geneList = sort(geneList, decreasing = T)
    }
  }
  else{
    geneList = sort(data, decreasing = T)
    if(inputType != "ENTREZID"){
      names(geneList) = clusterProfiler::bitr(names(geneList),
                                              fromType = inputType, toType = "ENTREZID",
                                              OrgDb = db)$ENTREZID
    }
    else{
      names(geneList) = data$gene
    }
  }
  
  
  # Pull gene sets of interest from msigdb
  m_df <- msigdbr(species = species, category = category, subcategory = subcategory) %>%
    dplyr::select(gs_name, entrez_gene)
  
  
  # remove NAs from values and names (NA names can cause a "duplicate gene name" warning)
  to_remove = sum(is.na(names(geneList)))
  if(to_remove > 0){
    warning(paste0("Warning: ", to_remove," entries failed to map from input to ENTREZID."))
  }
  geneList = na.omit(geneList)
  geneList = geneList[-which(is.na(names(geneList)))]
  
  # Run GSEA here
  set.seed(1234)
  output[['GSEA_df']] <- clusterProfiler::GSEA(geneList, exponent = exponent, TERM2GENE = m_df,
                                               pvalueCutoff = pval_cutoff, maxGSSize = maxGSSize,
                                               nPerm = nPerm,
                                               seed = TRUE)
  
  if(nrow(output[['GSEA_df']]) == 0){
    warning("No gene sets found under p-value, look into your p-value cutoff or shrinking your search space via category argument")
    return(output)
  }
  else{
    print(paste0("Number of significant sets: ",length(output[['GSEA_df']]$Description)))
    for(j in 1:min(length(output[['GSEA_df']]$Description), n_max)){
      set <- output[['GSEA_df']]$Description[j]
      output[[set]] <- enrichplot::gseaplot2(output[['GSEA_df']],
                                             geneSetID = j,
                                             title = output[['GSEA_df']]$Description[j],
                                             pvalue_table = TRUE,
                                             subplots = 1:2)
    }
    output$GSEA_df <- as.data.frame(output$GSEA_df@result)
    rownames(output$GSEA_df) <- NULL
  }
  return(output)
}