#' Transform Data (DEPRECATED, FUNCTIONALITY BUILT IN TO AGED.R ARGS)
#'
#' \code{verify_and_transform_data} is a helper function used in multiple core functions to check input, clear low variance and apply log or variance stabilizing transformation.
#'
#' @param data The original data that was plugged into the initial call of \code{AGED}
#
#' @param clear_low_variance A boolean variable that determines whether rows with var < 1 are removed or not. This is done before any type of transformation is performed.
#' 
#' @param transformation_type A string variable that determines whether or not a log or VST transformation should be done on the original data. If this argument is used, it should be "vst" or "log" only. If no transformation is to be performed, the default value or any string other than "vst" or "log" can be used.
#' 
#' @param blind If a VST is to be done, this boolean value determines whether it is blind or not.
#' 
#' @return Returns the transformed data.
#'
#' @export
#' 

verify_and_transform_data <- function(data, clear_low_variance = FALSE, transformation_type = "", blind = TRUE){
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
  } else if (transformation_type == "log") {
    print("Applying a log transformation...")
    data <- log1p(data)
  }
  return(data)
}
