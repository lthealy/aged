#' Cophenetic Correlation Coefficient Plot Generator
#'
#' \code{cophenetic_generator} will run non-negative matrix factorization to determine the cophenetic correlation coefficient for each rank of factorization in a desired range of ranks decided by the user. The cophenetic correlation coefficient can be helpful for the user in deciding what rank to use when running NMF. The raw cophenetic correlation coefficient value, the elbow method, or any other applicable approach can help determine a desirable rank for NMF. The higher the cophenetic correlation coefficient is, the more stable and reproducible the NMF results are. In the plot returned by this graph, the rank with the highest cophenetic correlation coefficient will be highlighted in red. If the input vector for \code{rank_range} is continuous, the rank directly before the biggest drop in cophenetic correlation coefficient, before any positive slopes, will be highlighted in cyan. If these two points are the same, the point will be highlighted in magenta. In the extremely rare event of a tie in numerical values, the first index is selected. However, it is ultimately up to the user to decide what rank is best fit for NMF runs. 
#' 
#' @param data Gene expression target data, a matrix-like object. The rows should represent genes, and each row must have a unique row name. Each column should represent a different sample.
#' 
#' @param rank_range Any numeric vector containing ranks of factorization to try (does not need to be continuous). Duplicates are removed, and the vector will be sorted in increasing order before use. All values should be positive and greater than 1.
#' 
#' @param nrun The desired number of NMF runs. For simply determing the cohpenetic correlation coefficient for each rank, it is not entirely necessary to perform a high number of runs or as many runs as normal when running NMF. This function defaults to 12, but any number of runs can be used.
#' 
#' @param nmf_seed The desired seed to be used for NMF
#' 
#' @param .options This argument is used to set runtime options. See \link[NMF]{nmf} for detailed information.
#' 
#' @param .pbackend This argument is used in accordance with the .options parameter. See \link[NMF]{nmf} for detailed information.
#' 
#' @param colors A boolean argument determining whether or not the specified points in the documentation (maximum value, point preceding the largest drop) should be highlighted in color. If TRUE, the points will be highlighted. If false, no points will be highlighted.
#' 
#' @param clear_low_variance A boolean variable that determines whether or not rows with variance less than 1 should be removed from the original dataset.
#' 
#' @param transformation_type A string variable that determines whether or not a log or VST transformation should be done on the original dataset. If this argument is intended to be used to perform a transformation, it should be "vst" or "log" only.
#' 
#' @param blind If a VST is to be done, this boolean value determines whether it is blind or not.
#'
#' @return A line graph that displays the cophenetic coefficients for the values in the range of ranks you selected. This function also plots it for you. 
#' 
#' @export
#' 
#' @import NMF
#' @import ggplot2
#' @import ggpubr
#' @import DESeq2

cophenetic_generator <- function(data, rank_range = 2:20, nrun = 12, nmf_seed = 123456, .options = "p4", .pbackend = "", colors = TRUE, clear_low_variance = FALSE, transformation_type = "", blind = TRUE) {
  
  # Validate data
  if (is.null(rownames(data))) {
    stop("The dataset must have row names for AGED to run properly. Please verify that your dataset has proper row names before continuing.")
  }
  
  # Clear low variance if desired
  if (clear_low_variance == TRUE) {
    print("Clearing low variance...")
    data <- data[apply(data, 1, var) > 1,]
  }
  
  # Perform desired transformation
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
  
  # Perform NMF
  if (.options == "" && .pbackend == "") {
    nmf_results <- NMF::nmfEstimateRank(data, range = rank_range, nrun = nrun, seed = nmf_seed)
  } else if (.options != "" && .pbackend == "") {
    nmf_results <- NMF::nmfEstimateRank(data, range = rank_range, nrun = nrun, seed = nmf_seed, .options = .options)
  } else if (.options == "" && .pbackend != "") {
    nmf_results <- NMF::nmfEstimateRank(data, range = rank_range, nrun = nrun, seed = nmf_seed, .pbackend = .pbackend)
  } else {
    nmf_results <- NMF::nmfEstimateRank(data, range = rank_range, nrun = nrun, seed = nmf_seed, .options= .options, .pbackend = .pbackend)
  }
  
  # Pull NMF results to create CCC plot
  df2 <- as.data.frame(cbind(nmf_results$measures$cophenetic, rank_range))
  colnames(df2) = c("cophenetic", "rank")
  if (colors == FALSE) {
    p <- ggplot(data = df2, aes(x = rank, y = cophenetic)) + xlab("Rank") + ylab("Cophenetic Correlation Coefficient") + ggtitle("Cophenetic Correlation Coefficient Plot") + geom_line() + geom_point() + theme_pubr()
    plot(p)
    return(p)
  } else {
    max_ccc <- which(df2$cophenetic==max(df2$cophenetic))
    if (length(max_ccc) > 1) {
      max_ccc <- max_ccc[1]
    }
    if (all(diff(rank_range) == 1)) {
      first_pos <- which(diff(df2$cophenetic) > 0)[1]
      if (is.na(first_pos)) {
        subset_coph <- diff(df2$cophenetic)
      } else {
        subset_coph <- diff(df2$cophenetic)[1:first_pos - 1]
      }
      subset_coph <- subset_coph * -1
      slope_drop_index <- which(subset_coph==max(subset_coph))
      if (length(slope_drop_index) > 1) {
        slope_drop_index <- slope_drop_index[1]
      }
      if (max_ccc == slope_drop_index) {
        p <- ggplot(data = df2, aes(x = rank, y = cophenetic)) + xlab("Rank") + ylab("Cophenetic Correlation Coefficient") + ggtitle("Cophenetic Correlation Coefficient Plot") + geom_line() + geom_point() + theme_pubr() + geom_point(x = df2$rank[max_ccc], y = df2$cophenetic[max_ccc], size = 5, color = "magenta2", alpha = 0.2)
      } else {
        p <- ggplot(data = df2, aes(x = rank, y = cophenetic)) + xlab("Rank") + ylab("Cophenetic Correlation Coefficient") + ggtitle("Cophenetic Correlation Coefficient Plot") + geom_line() + geom_point() + theme_pubr() + geom_point(x = df2$rank[max_ccc], y = df2$cophenetic[max_ccc], size = 5, color = "red", alpha = 0.2) + geom_point(x = df2$rank[slope_drop_index], y = df2$cophenetic[slope_drop_index], size = 5, color = "cyan", alpha = 0.2)
      }
    } else {
      p <- ggplot(data = df2, aes(x = rank, y = cophenetic)) + xlab("Rank") + ylab("Cophenetic Correlation Coefficient") + ggtitle("Cophenetic Correlation Coefficient Plot") + geom_line() + geom_point() + theme_pubr() + geom_point(x = df2$rank[max_ccc], y = df2$cophenetic[max_ccc], size = 5, color = "red", alpha = 0.2)
    }
    return(p)
  }
}
