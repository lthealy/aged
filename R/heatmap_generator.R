#' AGED Heatmap Generator
#'
#' \code{heatmap_generator} will pull metagene information from \code{AGED} results to create a heatmap. By default space for four components will be allocated in the plot. At the top left will be the color key, the top right will have the column dendrogram, the bottom left will have the row dendrogram, and the bottom right will have the image plot.
#'
#' @param aged_results The results of a successful call to \code{AGED}
#'
#' @param data The original data that was plugged into the initial call of \code{AGED}
#'
#' @param samp_info A matrix-like object containing sampling information regarding the samples used in the original AGED call. A proper samp_info object will have \code{n} rows corresponding to the \code{n} samples, and columns comprised of information pertaining to the samples such as extraction date, other characteristics, etc. 
#' 
#' @param batches A string vector that details the tracks of characteristics to outline on the heatmap. Each element of the string vector must correspond to a column name of samp_info.
#' 
#' @param clear_low_variance A boolean variable that determines whether rows with var < 1 are removed or not. This is done before any type of transformation is performed.
#' 
#' @param transformation_type A string variable that determines whether or not a log or VST transformation should be done on the original data. If this argument is used, it should be "vst" or "log" only. If no transformation is to be performed, the default value or any string other than "vst" or "log" can be used.
#' 
#' @param blind If a VST is to be done, this boolean value determines whether it is blind or not.
#' 
#' @param pearson A boolean value indicating whether or not pearson (row) clustering should be performed.
#' 
#' @param specific_order A vector of strings representing how the sample tracks (columns) should be sorted. Each element of this vector should represent a column name in \code{batches}. The sample tracks will be nested sorted by the first element first, then the second element, etc. 
#' 
#' @param legend A boolean value indicating whether or not a legend should be plotted alongside the heatmap. WARNING: Adding a row dendrogram will interfere with the legend if plotted concurrently.
#' 
#' @param hmap_color A string value or a vector of string values indicating the color(s) for the heatmap's scale. Low values on the scale are indicated by the first value in the vector, and high values on the scale will be indicated by the last value in the vector. Intermediate value(s) will be the color(s) indicated by this string or vector value.
#'
#' @param dendrogram A character string indicating whether to draw "none", "row", "column" or "both" dendrograms on the heatmap. Defaults to "column". This argument should only be one of these strings. WARNING: Adding a row dendrogram will interfere with the legend if plotted concurrently.
#' 
#' @param trace A character string indicating whether a solid trace line should be drawn across the "row"s, down the "column"s, "both", or "none". The distance of the line from the center of each color-cell is proportional to the size of the measurement. Defaults to "none". This argument should only be one of these strings.
#' 
#' @param scale A character string indicating if the values should be centered and scaled in the "row" direction, the "column" direction, or "none". The default is "none". This argument should only be one of these strings.
#' 
#' @param cexRow A positive number, used as 'cex.axis' for column axis labeling.
#' 
#' @param key A boolean value indicating whether or not a color key should be drawn.
#' 
#' @param lhei A numerical vector indicating the relative row heights of the rows of the plot.
#' 
#' @param lwid A numerical vector indicating the relative column widths of the rows of the plot.
#' 
#' @param legend_size A numerical value indicating the size of the legend. The default value is 0.75. Increasing this value will increase the text and block size of the legend.
#' 
#' @param legend_space A numerical value indicating how far spaced apart different elements in the legend should be. The default value is 1. Increasing this value will increase the space between different elements in the legend.
#' 
#' @param ... Optional arguments that can be passed to \code{heatmap.3}. For a full list of arguments, check \href{https://github.com/obigriffith/biostar-tutorials/blob/master/Heatmaps/heatmap.3.R}{here}. 
#' 
#' @return Plots the requested heatmap
#'
#' @export
#' 
#' @import DESeq2
#' @import bioDist

heatmap_generator <- function(aged_results, data, samp_info, batches = names(samp_info), clear_low_variance = FALSE, transformation_type = "", blind = TRUE, pearson = FALSE, specific_order = NULL, legend = TRUE, hmap_color = "skyblue", dendrogram = "none", trace = "none", scale = "row", cexRow = 0.5, key = FALSE, lhei = c(1,3), lwid = c(2,3), legend_size = 0.75, legend_space = 1, ...) {
  
  # Verify and prepare data
  if (is.null(rownames(data))) {
    stop("The data must have row names for AGED to run properly. Please verify that your data has proper row names before continuing.")
  }
  
  # Clear low variance
  if (clear_low_variance == TRUE) {
    print("Clearing low variance...")
    data <- data[apply(data, 1, var) > 1,]
  }
  
  # Requested transformation
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
  
  # Set up "featInfo"
  rank <- length(aged_results) - 1
  rn <- row.names(data)
  df <- data.frame(matrix(ncol=2,nrow=length(rn)))
  df[,1] <- rn
  colnames(df) <- c("genes", "metagenes")
  for (i in 1:rank) {
    nms <- names(aged_results[[i]][[1]])
    for (m in 1:(length(nms))) {
      index <- match(nms[m], df[,1])
      df[index, 2] <- i
    }
  }
  
  # For loop in place for adaptability of plotting multiple AGED results on same heatmap in future, for now should just iterate once
  for(j in names(df)[-1]){
    metagenes = unique(which(!is.na(df[,j])))
    data <- data[metagenes,]
    df <- df[metagenes,]
    
    # Remove any genes with SD of 0 (can't cluster)
    no_var = which(apply(data,1,sd) == 0)
    if(length(no_var) != 0){
      data <- data[-no_var,]
      df <- df[-no_var,]
    }
    
    # Remove any samples with SD of 0 (can't cluster)
    no_var = which(apply(data,2,sd) == 0)
    if(length(no_var) != 0){
      data <- data[,-no_var]
      df <- df[-no_var,]
    }
    
    # Prepare colors
    colorlists <- rep(list(c("gray94", "blue", "green",
                             "yellow", "orange", "red","black")),
                      length(batches))
    ColSideColors <- aged::get_side_colors(sampInfo = samp_info,
                                           sampleTracks = batches,
                                           drop.levels = TRUE,
                                           colorlists = colorlists,
                                           displaynames = batches)
    df$metagenes <- as.factor(df$metagenes)
    colorlists <- rep(list(c("gray94", "blue", "green",
                             "yellow", "orange", "red","black")),
                      rank)
    RowSideColors <- aged::get_side_colors(sampInfo = df,
                                           sampleTracks = "metagenes", 
                                           drop.levels = TRUE,
                                           colorlists = colorlists,
                                           displaynames = strsplit(j,"_")[[1]][1])
    
    # Adjusts track (column) sorting depending on user preference
    if (!is.null(specific_order)) {
      specific_order <- paste(specific_order, collapse = ", ")
      if (grepl(":", specific_order, fixed = TRUE) == TRUE) {
        specific_order = paste('\`', specific_order, '\`', sep = "")
      }
      sampleOrder <- with(samp_info, eval(parse(text = paste("order(", specific_order, ")"))))
      Colv <- aged::convert_order_to_dendrogram(sampleOrder)
    } else {
      Colv = as.dendrogram( hclust( bioDist::cor.dist(x = t(as.matrix(data))) ) )
    }
    
    # Adjusts row sorting depending on user preference
    geneOrder <- order(df[,j])
    if (pearson == TRUE) {
      Rowv = as.dendrogram( hclust( bioDist::cor.dist(x = as.matrix(data)) ) )
    } else {
      Rowv <- aged::convert_order_to_dendrogram(geneOrder)
    }
    
    if (length(hmap_color) > 1) {
      hmap_color = unlist(hmap_color)
    }
    
    # Establishes default values
    l <- list(...)
    if (is.null(l$Colv)) l$Colv <- Colv
    if (is.null(l$Rowv)) l$Rowv <- Rowv
    if (is.null(l$dendrogram)) l$dendrogram <- dendrogram
    if (is.null(l$trace)) l$trace <- trace
    if (is.null(l$scale)) l$scale <- scale
    if (is.null(l$labRow)) l$labRow <- NA
    if (is.null(l$col)) l$col <- colorRampPalette(c(hmap_color))(n = 299)
    if (is.null(l$ColSideColors)) l$ColSideColors <- ColSideColors$SideColors
    if (is.null(l$ColSideColorsSize)) l$ColSideColorsSize <- dim(ColSideColors$SideColors)[2]*1.2
    if (is.null(l$RowSideColors)) l$RowSideColors <- t(RowSideColors$SideColors)
    if (is.null(l$RowSideColorsSize)) l$RowSideColorsSize <-dim(RowSideColors$SideColors)[2]
    if (is.null(l$cexRow)) l$cexRow <- cexRow
    if (is.null(l$key)) l$key <- key
    if (is.null(l$lhei)) l$lhei <- lhei
    if (is.null(l$lhei)) l$lwid <- lwid
    breaks <- (((0:299)/299)-0.5)*4
    
    # Plots heatmap
    aged::heatmap.3(x = data,
                    Colv = l$Colv, 
                    Rowv = l$Rowv,
                    dendrogram = l$dendrogram,
                    trace = l$trace,
                    scale = l$scale,
                    labRow = l$labRow,
                    col = l$col,
                    ColSideColors = l$ColSideColors,
                    ColSideColorsSize = l$ColSideColorsSize,
                    RowSideColors = l$RowSideColors,
                    RowSideColorsSize = l$RowSideColorsSize,
                    cexRow = l$cexRow,
                    key = l$key,
                    lhei = l$lhei,
                    lwid = l$lwid,
                    breaks = breaks,
                    ...)
    if (legend == TRUE) {
      legend(xy.coords(x=0,y=1),
            legend=c("Legend","",ColSideColors$text,"Metagene Tags","",RowSideColors$text),
            fill=c("white","white",ColSideColors$colors,"white","white",RowSideColors$colors),
            border=FALSE, bty="n",
            y.intersp = legend_space, cex = legend_size)
    }
  }
}