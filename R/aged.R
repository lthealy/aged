#' Automatic Gene Expression Deconvolution Using FaStaNMF
#'
#' \code{AGED} will perform gene expression deconvolution on a matrix or data frame of gene expression data using FaStaNMF.
#' 
#' @param data Gene expression target data, a matrix-like object. The rows should represent genes, and each row must have a unique row name. Each column should represent a different sample.
#' 
#' @param rank The factorization rank (number of factors) to be used during NMF. This function argument should be a positive integer value.
#' 
#' @param n The number of barcode genes desired per metagene
#' 
#' @param nrun The desired number of NMF runs to be run on the initially reduced dataset.
#' 
#' @param nmf_seed The desired seed to be used for NMF on the initially reduced dataset.
#' 
#' @param mvg A numerical argument determining how many of the most variable genes to look at during the first steps of FaStaNMF.
#' 
#' @param clv A numerical value \code{x} that reduces the dataset by removing genes with variance < \code{x} across all samples. Our recommended value is to set this parameter to 1 if genes expression low variance across samples is desired. These genes will not be considered at all for the deconvolution. This is done before any type of transformation or other reduction is performed.
#' 
#' @param transformation A numerical value that determines whether or not a log or VST transformation should be done on the original dataset. A value of 0 indicates no transformation, a value of 1 indicates a log transformation using \link[base]{log1p}, a value of 2 indicates a VST transformation using \link[DESeq2]{varianceStabilizingTransformation} If this argument is used, it should be "0", "1" or "2" only. Any other value will assume no transformation. For FaStaNMF, untransformed data should be log-transformed or VST-transformed.
#' 
#' @param blind If a VST is to be done using the \code{transformation} parameter, this boolean value determines whether it is blind or not.
#' 
#' @param ... Other arguments to be passed to FaStaNMF.
#' 
#' @return A list containing barcode genes for each metagene, and the raw W and H matrix returned by FaStaNMF.
#' 
#' @export
#' 
#' @import NMF
#' @import DESeq2

aged <- function(data, rank, n = 25, nrun = 200, nmf_seed = 123456, mvg = 1000, clv = 0, transformation = 0, blind = TRUE, ...) {
   
   # The pipeline does not start if the matrix-like object is without row names.
   if (is.null(rownames(data))) {
      stop("In order for genes to be identifiable, the dataset must have row names for AGED to run properly. Please verify that your dataset has proper row names before continuing.")
   }
   
   # Clear low variance if desired.
   if (clv != 0) {
      print("Clearing low variance...")
      data <- data[apply(data, 1, var) > clv,]
   }
   
   # Apply a transformation if desired for untransformed data.
   if (transformation == 2) {
      print("Applying a variance-stabilizing transformation...")
      data <- DESeq2::varianceStabilizingTransformation(data, blind = blind)
      
      # Prevents conflicts between the two different seed generics from NMF and SummarizedExperiment.
      detach("package:DESeq2")
      detach("package:SummarizedExperiment")
      detach("package:DelayedArray")
   } else if (transformation == 1) {
      print("Applying a log transformation...")
      data <- log1p(data)
   }
   
   # Start NMF
   print(paste("Starting FaStaNMF with rank ",rank,"...", sep = ""))
   nmf_object <- FaStaNMF::fastanmf(data, rank = rank, nrun = nrun, nmf_seed = nmf_seed, mvg = mvg, ...)
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
   metagene_list[[1]] <- nmf_object
   names(metagene_list)[[1]] <- "nmf_object"
   for (i in 1:ncol(differenceMatrix)) {
      differenceMatrix <- differenceMatrix[order(differenceMatrix[,i], decreasing = TRUE),]
      metagene_list[[i + 1]] <- differenceMatrix[,i][1:n]
      names(metagene_list)[[i + 1]] <- paste0("metagene", as.character(i))
   }
   return(metagene_list)
}