#' AGED Scatterplot Generator
#'
#' \code{scatterplot_generator} will pull information from \code{AGED} results to plot gene expression for each sample for a selected gene. The value in the original dataset will be plotted on the y-axis, and the user's choice of data taken from either the W matrix, H matrix, or W * H will be plotted on the y-axis. The column or row taken from the H or W matrix is chosen by which metagene the selected gene is the most prevalent in. 
#'
#' @param aged_results The results of a successful call to \code{AGED}
#'
#' @param data The original dataset that was plugged into the initial call of \code{AGED}
#'
#' @param gene A string value perfectly corresponding to a row value of the original dataset, representing a gene name.
#' 
#' @param clv A numerical value \code{x} that reduces the dataset by removing genes with variance < \code{x} across all samples. Our recommended value is to set this parameter to 1 if genes expression low variance across samples is desired. These genes will not be considered at all for the deconvolution. This is done before any type of transformation or other reduction is performed.
#' 
#' @param transformation A numerical value that determines whether or not a log or VST transformation should be done on the original dataset. A value of 0 indicates no transformation, a value of 1 indicates a log transformation using \link[base]{log1p}, a value of 2 indicates a VST transformation using \link[DESeq2]{varianceStabilizingTransformation} If this argument is used, it should be "0", "1" or "2" only. Any other value will assume no transformation. For FaStaNMF, untransformed data should be log-transformed or VST-transformed.
#' 
#' @param blind If a VST is to be done, this boolean value determines whether it is blind or not.
#' 
#' @param x_axis A string value representing what data should be plotted. The value should be one of: "h" or "wh" only.
#' 
#' @param color A string vector representation of size \code{n} where \code{n} equals the number of samples. This vector should detail how each sample should be categorized and colored within each scatterplot. Using a column from a dataset's sample information dataframe is common. For example, one might want to separate points in the scatterplots by treatment or batch. Each treatment, batch, or category will be represented by a different color.
#' 
#' @param shape A string vector representation of size \code{n} where \code{n} equals the number of samples. This vector should detail how each sample should be categorized and represented within each scatterplot. Using a column from a dataset's sample information dataframe is common. For example, one might want to separate points in the scatterplots by treatment or batch. Each treatment, batch, or category will be represented by a different shape.
#' 
#' @param reg A boolean value representing whether a solid, linear regression line should be added to the plot.
#' 
#' @param reg_color If a linear regression line is added to the scatterplot, this string value should represent the color of that line.
#' 
#' @param xy A boolean value representing whether or not a dotted x = y should be added to the plot. This is typically not applicable to "h". 
#' 
#' @param xy_color If an x = y line is added to the scatterplot, this string value should represent the color of that line.
#' 
#' @param ellipse A boolean value representing whether or not ellipses should be added to the scatterplot. These ellipses will create normal confidence ellipses based on the "color" vector parameter. 
#' 
#' @return Returns a scatterplot containing the desired gene expression values. This function also plots it for you. 
#'
#' @export
#' @import ggplot2
#' @import ggbeeswarm

scatterplot_generator <- function(aged_results, data, gene, clv = 0, transformation = 0, blind = TRUE, x_axis = "wh", color = NULL, shape = NULL, reg = TRUE, reg_color = "black", xy = TRUE, xy_color = "gray", ellipse = TRUE) {
  
  # Validate data
  if (is.null(rownames(data))) {
    stop("The dataset must have row names for AGED to run properly. Please verify that your dataset has proper row names before continuing.")
  }
  
  # Clear low variance if desired
  if (clv > 0) {
    print("Clearing low variance...")
    data <- data[apply(data, 1, var) > clv,]
  }
  
  # Perform desired transformation
  if (transformation == 2) {
    print("Applying a variance-stabilizing transformation...")
    data <- DESeq2::varianceStabilizingTransformation(data, blind = blind)
    if("DESeq2" %in% (.packages())){
      detach("package:DESeq2", unload=TRUE) 
    }
    if("SummarizedExperiment" %in% (.packages())){
      detach("package:SummarizedExperiment", unload=TRUE) 
    }
    if("DelayedArray" %in% (.packages())){
      detach("package:DelayedArray", unload=TRUE) 
    }
  } else if (transformation == 1) {
    print("Applying a log transformation...")
    data <- log1p(data)
  }
  
  # Pull NMF results and find proper metagene for selected gene
  rank <- length(aged_results) - 1
  w <- aged_results$nmf_object@fit@W
  h <- aged_results$nmf_object@fit@H
  wh <- (w %*% h)
  if (!any(row.names(data) == gene)) {
    stop(paste("No gene",gene,"was found inside of the dataset.", sep = " "))
  }
  gene_row = data[gene,]
  lst <- vector()
  for (i in 1:rank) {
    w <- w[order(w[,i], decreasing = TRUE),]
    rn <- row.names(w)
    gene_index <- which(rn==gene)
    lst[i] = gene_index
  }
  
  # Create scatterplot
  min_index <- which(lst==min(lst))
  df2 <- as.numeric(data[gene,])
  x_lab = ""
  if (x_axis == "h") {
    df2 <- rbind(df2, as.numeric(h[min_index,]))
    x_lab = paste(gene,"Values, H Matrix, Metagene",min_index, sep = " ")
  } else if (x_axis == "wh") {
    df2 <- rbind(df2, as.numeric(wh[gene,]))
    x_lab = paste(gene,"Values, W * H Matrix (Expected Values)", sep = " ")
  } else {
    stop("The x_axis parameter must be one of \"h\" or \"wh\" only.")
  }
  df2 <- t(df2)
  df2 <- as.data.frame(df2)
  colnames(df2) <- c("original", "chosen")
  
  g <- qplot(x = chosen, y = original, data = df2, xlab = x_lab, ylab = paste(gene,"Values, Original Dataset", sep = " "), shape = shape, color = color, geom = c("point")) + theme_pubr()
  if (reg == TRUE) {
    g <- g + geom_smooth(color = reg_color, method='lm', formula=y~x, aes(group = 1), se = FALSE)
  }
  if (xy == TRUE) {
    g <- g + geom_abline(show.legend = T, color = xy_color, slope=1, lty=2) 
  }
  if (ellipse == TRUE) {
    g <- g + stat_ellipse(aes(group = color))
  }
  return(g)
}