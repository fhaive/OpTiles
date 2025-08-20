#' Analyze and Plot Distribution of Coverage or Beta Values from Methylation Data
#'
#' This function computes and plots the distribution of either read coverage or beta values (percent methylation) 
#' from methylation data. It validates inputs and requires that, when investigating beta values, 
#' the input data frame columns match sample names provided in the metadata.
#'
#' @param methdata_df A data frame containing methylation data.
#'   - If `variable_to_investigate = "coverage"`, it should contain columns with "coverage" in their names.
#'   - If `variable_to_investigate = "betavalue"`, columns should correspond to sample names from metadata.
#' @param variable_to_investigate A character string specifying which variable to analyze: either `"coverage"` or `"betavalue"`.
#' @param position_text Numeric scalar indicating the X-axis position where statistics text will be placed in the plot. Default is 50.
#' @param metadata Either a data frame or character vector containing sample metadata.
#'   Required if `variable_to_investigate = "betavalue"`.
#' @param sample_name_variable Character string specifying the column name in `metadata` that contains sample names.
#'   Required if `variable_to_investigate = "betavalue"`.
#'
#' @return No return value. The function generates and displays a histogram plot of the selected variable distribution.
#'
#' @details
#' When `variable_to_investigate` is `"coverage"`, the function extracts coverage columns from `methdata_df` 
#' and plots their distribution.  
#' When `variable_to_investigate` is `"betavalue"`, it checks that `methdata_df` columns match sample names in `metadata` 
#' (using `sample_name_variable`), then plots the beta values distribution.
#'
#' @importFrom dplyr pull
#' @importFrom tidyr pivot_longer
#' @export
analyze_methylation_stats <- function(
    methdata_df, 
    variable_to_investigate, 
    metadata = NULL,
    sample_name_variable = NULL,
    position_text= 50) 
  {
  valid_options <- c("coverage", "betavalue")
  if (!variable_to_investigate %in% valid_options) {
    stop(sprintf("Invalid 'variable_to_investigate': '%s'. Choose either 'coverage' or 'betavalue'.", variable_to_investigate))
  }
  
  transform_in_vector <- function(df) {
    df %>%
      tidyr::pivot_longer(cols = everything(), values_to = "value") %>%
      dplyr::pull(value)
  }
  
  if (variable_to_investigate == "coverage") {
    meth_cov <- methdata_df[, grep("coverage", colnames(methdata_df)), drop = FALSE]
    if (ncol(meth_cov) == 0) {
      stop("No columns containing 'coverage' found in 'methdata_df'.")
    }
    meth_cov_vect <- transform_in_vector(meth_cov)
    
    do_histogram_with_statistics(
      vect = meth_cov_vect,
      title = "Coverage Distribution",
      x_axis = "Coverage",
      text_position = position_text,
      color_plot = "salmon"
    )
    
  } else if (variable_to_investigate == "betavalue") {
    if (!is.data.frame(methdata_df)) {
      stop("'methdata_df' must be a data.frame when 'variable_to_investigate' is 'betavalue'.")
    }
    if (is.null(metadata)) {
      stop("'metadata' must be provided when 'variable_to_investigate' is 'betavalue'.")
    }
    if (is.null(sample_name_variable)) {
      stop("'sample_name_variable' must be provided when 'variable_to_investigate' is 'betavalue'. It should be the column name in 'metadata' containing sample names.")
    }
    
    sample_names <- NULL
    if (is.data.frame(metadata)) {
      if (!(sample_name_variable %in% colnames(metadata))) {
        stop(sprintf("Column '%s' not found in 'metadata'.", sample_name_variable))
      }
      sample_names <- as.character(metadata[[sample_name_variable]])
    } else if (is.character(metadata)) {
      sample_names <- metadata
    } else {
      stop("'metadata' must be either a data.frame or a character vector.")
    }
    
    missing_samples <- dplyr::setdiff(colnames(methdata_df), sample_names)
    if (length(missing_samples) > 0) {
      stop(sprintf(
        "The following columns in 'methdata_df' are not present in 'metadata' sample names: %s",
        paste(missing_samples, collapse = ", ")
      ))
    }
    
    betavalue_vect <- transform_in_vector(methdata_df)
    
    do_histogram_with_statistics(
      vect = betavalue_vect,
      title = "Betavalues Distribution",
      x_axis = "Betavalues",
      text_position = position_text,
      color_plot = "tomato"
    )
  }
}
