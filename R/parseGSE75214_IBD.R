#' Parse IBD Dataset
#' 
#' \code{parseGSE75214} parses the IBD data fileand allows it to be ready for use by \code{AGED}.
#'
#' @export
#' @import GEOquery
#' @import readr


#-----------------------------GSE75214-----------------------------------

parseGSE75214 <- function() {
  IBD_dataset <- list()
  #Load expression from supplement
  IBD_dataset$ex <- read_delim(system.file("extdata/entire_dataset.txt", package = "aged"), 
                                                "\t", escape_double = FALSE, trim_ws = TRUE) %>%
    tibble::column_to_rownames("GeneSymbol") %>% as.matrix()
  ## Use GEOquery for sample info
  gse <- getGEO("GSE75214", 
                GSEMatrix = TRUE, 
                getGPL = FALSE)$GSE75214_series_matrix.txt.gz
  
  #parse sampInfo
  IBD_dataset$sampInfo <- data.frame(gse@phenoData@data[,c("title","characteristics_ch1","characteristics_ch1.1",
                                                           "characteristics_ch1.2","geo_accession")],
                                     row.names = NULL)
  colnames(IBD_dataset$ex) = IBD_dataset$sampInfo$title
  
  IBD_dataset$sampInfo$characteristics_ch1 = stringr::str_split_fixed(IBD_dataset$sampInfo$characteristics_ch1,": ",2)[,2]
  IBD_dataset$sampInfo$characteristics_ch1.1 = stringr::str_split_fixed(IBD_dataset$sampInfo$characteristics_ch1.1,": ",2)[,2]
  IBD_dataset$sampInfo$characteristics_ch1.1[which(IBD_dataset$sampInfo$characteristics_ch1.1 == "CD")] = "Crohn's"
  IBD_dataset$sampInfo$characteristics_ch1.1[which(IBD_dataset$sampInfo$characteristics_ch1.1 == "Crohn's disease")] = "Crohn's"
  IBD_dataset$sampInfo$characteristics_ch1.2 = stringr::str_split_fixed(IBD_dataset$sampInfo$characteristics_ch1.2,": ",2)[,2]
  IBD_dataset$sampInfo$cohort = paste(IBD_dataset$sampInfo$characteristics_ch1,
                                     IBD_dataset$sampInfo$characteristics_ch1.1,
                                     IBD_dataset$sampInfo$characteristics_ch1.2,
                                     sep = "_")
  
  IBD_dataset$sampInfo = as.data.frame(apply(IBD_dataset$sampInfo,2,as.factor))
  colnames(IBD_dataset$sampInfo)[2:4] = c("tissue","disease","activity")
  
  IBD_dataset$featInfo <- data.frame(SYMBOL = rownames(IBD_dataset$ex))
  IBD_dataset$metadata <- list(log.transformed = FALSE,
                                  reference = "Vancamelbeke M, Vanuytsel T, Farré R, Verstockt S et al. Genetic and Transcriptomic Bases of Intestinal Epithelial Barrier Dysfunction in Inflammatory Bowel Disease. Inflamm Bowel Dis 2017 Oct;23(10):1718-1729. PMID: 28885228",
                                  accession = "GSE75214",
                                  description = "Mucosal biopsies from colon of 97 UC patients, 8 CD patients, 11 controls and mucosal biopsies from illeum of 67 CD patients, 11 controls",
                                  survivalA = "None",
                                  survivalB = "None")
  str(IBD_dataset)
  save(IBD_dataset, file = "./data/IBD_dataset.RData", compress = T)
  return(NULL)
}

