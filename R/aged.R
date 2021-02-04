#' Automatic Gene Expression Deconvolution
#'
#' \code{AGED} will perform gene expression deconvolution on a matrix or data frame of gene counts using non-negative matrix factorization. The results of the gene expression deconvolution will then be interpreted using GSEA.
#' 
#' @param data Gene expression target data, a matrix-like object. The rows should represent genes, and each row must have a unique row name. Each column should represent a different sample.
#' 
#' @param rank The factorization rank (number of factors) to be used during NMF. This function argument should usually be a positive integer value.
#' 
#' @param n The number of barcode genes desired per metagene
#' 
#' @param nrun The desired number of NMF runs
#' 
#' @param nmf_seed The desired seed to be used for NMF
#' 
#' @param .options This argument is used to set runtime options. See \link[NMF]{nmf} for detailed information.
#' 
#' @param .pbackend This argument is used in accordance with the .options parameter. See \link[NMF]{nmf} for detailed information.
#' 
#' @param species The species corresponding to the dataset (human, mouse, etc.).
#' 
#' @param exponent The weight of each step between differentially-ranked genes in GSEA.
#' 
#' @param gsea_barcodes A boolean value that determines which values are passed to GSEA. If \code{TRUE}, only the barcode genes and their corresponding numerical values will be used as input for GSEA for each metagene. If \code{FALSE}, all genes and their corresponding numeric values will be input into GSEA for each metagene.
#' 
#' @param max_geneset_size The maximum number of genes in potential gene sets in gene set enrichment analysis.
#' 
#' @param category MSigDB collection abbreviation, such as H, C1, C2, C3, C4, C5, C6, C7.
#' 
#' @param subcategory MSigDB sub-collection abbreviation, such as CGP or BP.
#' 
#' @param input The input type for gene names such as Entrez, Ensembl, etc.
#' 
#' @param n_max Maximum number of plots to return.
#' 
#' @param pval_cutoff The p-value cutoff for which GSEA results are to be considered statistically significant.
#' 
#' @param nperm The number of permutations to be performed by gene set enrichment analysis.
#' 
#' @param clear_low_variance A boolean variable that determines whether rows with var < 1 are removed or not. This is done before any type of transformation is performed.
#' 
#' @param transformation_type A string variable that determines whether or not a log or VST transformation should be done on the original dataset. If this argument is used, it should be "vst" or "log" only. If no transformation is to be performed, the default value or any string other than "vst" or "log" can be used. For NMF, untransformed data should be log-transformed or VST-transformed.
#' 
#' @param blind If a VST is to be done, this boolean value determines whether it is blind or not.
#' 
#' @return A list containing barcode genes for each metagene, gene set enrichment analysis results for each metagene, and the raw W and H matrix returned by NMF.
#' 
#' @export
#' 
#' @import NMF
#' @import org.Mm.eg.db
#' @import org.Hs.eg.db
#' @import clusterProfiler
#' @import dplyr
#' @import msigdbr
#' @import DESeq2

aged <- function(data, rank, n = 25, nrun = 30, nmf_seed = 123456, .options = "p4", .pbackend = "", species = "Homo sapiens", exponent = 0, gsea_barcodes = TRUE, max_geneset_size = 200, category = NULL, subcategory = NULL, input = "SYMBOL", n_max = 10, pval_cutoff = 0.05, nperm = 50000, clear_low_variance = FALSE, transformation_type = "", blind = TRUE) {
   if (is.null(rownames(data))) {
      stop("The dataset must have row names for AGED to run properly. Please verify that your dataset has proper row names before continuing.")
   }
   if (clear_low_variance == TRUE) {
      print("Clearing low variance...")
      data <- data[apply(data, 1, var) > 1,]
   }
   if (transformation_type == "vst") {
      print("Applying a variance-stabilizing transformation...")
      data <- DESeq2::varianceStabilizingTransformation(data, blind = blind)
      detach("package:DESeq2")
      detach("package:SummarizedExperiment")
      detach("package:DelayedArray")
   } else if (transformation_type == "log") {
      print("Applying a log transformation...")
      data <- log1p(data)
   }
   
   print(paste("Performing NMF with rank ",rank,"...", sep = ""))
   if (.options == "" && .pbackend == "") {
      nmf_object <- NMF::nmf(data, rank = rank, nrun = nrun, seed = nmf_seed)
   } else if (.options != "" && .pbackend == "") {
      nmf_object <- NMF::nmf(data, rank = rank, nrun = nrun, seed = nmf_seed, .options = .options)
   } else if (.options == "" && .pbackend != "") {
      nmf_object <- NMF::nmf(data, rank = rank, nrun = nrun, seed = nmf_seed, .pbackend = .pbackend)
   } else {
      nmf_object <- NMF::nmf(data, rank = rank, nrun = nrun, seed = nmf_seed, .options= .options, .pbackend = .pbackend)
   }
   h <- nmf_object@fit@H
   w <- nmf_object@fit@W
   col_sums <- colSums(w)
   lambda_h <- matrix(,nrow=nrow(h),ncol=ncol(h))
   lambda_w <- matrix(,nrow=nrow(w), ncol=ncol(w))
   row.names(lambda_w) <- row.names(w)

   # Normalization
   print("Normalizing NMF results...")
   for (i in 1:rank) {
      desired_col <- w[,i]
      reciprocalColSum <- 1 / col_sums[i]
      lambda_w[,i] <- w[,i] * reciprocalColSum
      lambda_h[i,] <- h[i,] * (1 / reciprocalColSum)
   }

   # Calculate the barcode genes for each metagene
   print("Calculating barcode genes for each metagene...")
   if (rank == 2) {
      lambda_w <- cbind(lambda_w, lambda_w[,1] - lambda_w[,2])
      lambda_w <- cbind(lambda_w, lambda_w[,2] - lambda_w[,1])
      lambda_w <- lambda_w[,-1:-2]
      differenceMatrix <- lambda_w
   } else {
      differenceMatrix <- matrix(,nrow=nrow(lambda_w), ncol=ncol(lambda_w))
      row.names(differenceMatrix) <- row.names(w)
      for (i in 1:nrow(lambda_w)) {
         ro <- lambda_w[i,]
         for (j in 1:ncol(lambda_w)) {
            ro2 <- ro[-j]
            max_value <- max(ro2)
            differenceMatrix[i,j] <- lambda_w[i,j] - max_value
         }
      }
   }
   metagene_list <- list()
   print("Starting GSEA...")
   for (i in 1:rank) {
      metagene_name <- paste("metagene", i, sep="")
      assign(metagene_name, list())
      differenceMatrix <- differenceMatrix[order(differenceMatrix[,i], decreasing = TRUE),]
      rn <- row.names(differenceMatrix)
      row_values <- differenceMatrix[,i]
      names(row_values) <- rn
      lst <- list()
      lst[[1]] <- row_values[1:n]
      names(lst)[[1]] <- paste("barcode_genes")
      if (gsea_barcodes == TRUE) {
         row_values <- row_values[row_values >= 0]
      }
      lst[[2]] <- aged::cluster_wrapper(row_values, exponent = exponent, maxGSSize = max_geneset_size, species = species, category = category, subcategory = subcategory, input = input, n_max = n_max, pval_cutoff = pval_cutoff, nPerm = nperm)
      names(lst)[[2]] <- "gsea"
      metagene_list[[i]] <- lst
      names(metagene_list)[[i]] <- paste("metagene",i,sep="") 
   }
   metagene_list[[rank + 1]] <- w
   names(metagene_list)[[rank + 1]] <- "w"
   metagene_list[[rank + 2]] <- h
   names(metagene_list)[[rank + 2]] <- "h"
   return(metagene_list)
}
