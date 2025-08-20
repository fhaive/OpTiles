
#' Plot difference in methylation (diff.meth) Distribution from differential methylation
#'
#' This function generates histograms of diff.meth from differential 
#' methylation analysis results, allowing visualization of their distribution.
#'
#' @param DM_results A list containing differential methylation results. 
#' @param filtering A logical value indicating whether to use the filtered results (`TRUE`) 
#'                  or unfiltered results (`FALSE`).
#' @param threshold_pval A numeric value specifying the p-value threshold used for filtering (if applicable).
#'
#' @return A combined plot displaying histograms of diff.meth values for each comparison.
#'
#' @details
#' - If \code{filtering} is `TRUE`, the function plots significant results  
#' - If \code{filtering} is `FALSE`, the function plots all results.
#' - The histograms are arranged in a 3-column layout using \code{patchwork}.
#'
#' @examples
#' \dontrun{
#'   plot_FC_distribution(DM_results, filtering = TRUE, threshold_pval = 0.05)
#' }
#'
#' @import ggplot2
#' @import patchwork
#' @export

plot_FC_distribution <- function(DM_results, filtering, threshold_pval) {
  library(patchwork)  # Ensure patchwork is loaded for plot layout
  
  # Determine whether to use filtered or unfiltered Limma results
  if (filtering) {
    subtitle <- paste0("p-value < ", threshold_pval)
  } else {
    subtitle <- ""
  }
  plot_data = DM_results
  # Check if plot_data is empty
  if (length(plot_data) == 0) {
    stop("No data available for plotting. Check Limma results and filtering conditions.")
  }
  
  # Generate histograms for log fold change (logFC) distribution in different samples
  plot_list <- lapply(names(plot_data), function(name) {
    df <- plot_data[[name]]
    
    # Check if diff.meth column exists in the dataset
    if (!"diff.meth" %in% colnames(df)) {
      stop(paste("Column 'diff.meth' not found in dataset:", name))
    }
    
    ggplot(df, aes(x = diff.meth)) +
      geom_histogram(fill = "lightblue", bins = 200) +
      labs(title = paste("diff.meth", name),
           subtitle = subtitle,
           x = "diff.meth",
           y = "Frequency") +
      theme_minimal()
  })
  
  # Combine all histograms into a single layout (3 columns)
  final_plot <- Reduce(`+`, plot_list) & plot_layout(ncol = 3, widths = rep(1, 3))
  
  # Print the final combined plot
  print(final_plot)
  
  # Return the plot for further use if needed
  return(final_plot)
}


#' Combine Results into a Single Data Frame
#'
#' This function merges multiple differential methylation results into a single data frame.
#' The resulting data frame contains columns for diff.meth and adjusted p-values 
#' (qvalue) for each comparison.
#'
#' @param DM_results A list containing differential methylation results.
#' @param output_folder A character string specifying the directory where the output file should be saved.
#'
#' @return A combined data frame where each row represents a genomic location and each column 
#' contains diffmeth and adjusted p-value
#' for different comparisons.
#'
#' @details
#' - The function ensures that genomic locations are properly aligned across all comparisons.
#' - The resulting data frame is saved as an `.RData` file in the specified output folder.
#'
#' @examples
#' \dontrun{
#'   combine_into_df(DM_results, output_folder = "results/")
#' }
#'
#' @import dplyr
#' @export

combine_into_df <- function(DM_results, output_folder) {
  

    combine_data <- DM_results

  
  # Check if the Limma results are empty
  if (length(combine_data) == 0) {
    stop("No data available for merging. Check Limma results and filtering conditions.")
  }
  
  # Merge all data frames by "location"
  combined_df <- Reduce(function(x, y) merge(x, y, by = "location", all = TRUE), 
                        lapply(names(combine_data), function(name) {
                          df <- combine_data[[name]]
                          
                          # Ensure required columns exist
                          if (!all(c("location", "diff.meth", "qvalue") %in% colnames(df))) {
                            stop(paste("Missing required columns in dataset:", name))
                          }
                          
                          df <- df[, c("location", "diff.meth", "qvalue")]
                          colnames(df)[colnames(df) == "diff.meth"] <- paste0("diff.meth_", name)
                          colnames(df)[colnames(df) == "qvalue"] <- paste0("qvalue", name)
                          return(df)
                        }))
  
  # Save the combined data frame
  save(combined_df, file = paste0(output_folder, "/DM_results_df.RData"))
  
  # Return the combined data frame
  return(combined_df)
}
