#' Filter CpG Sites by Beta Value Range Across Groups
#'
#' This function filters a beta value matrix by computing the difference (max - min) 
#' in beta values within defined groups. CpGs with a difference greater than 
#' or equal to a specified threshold are excluded.
#'
#' @param betavalue_df A `data.frame` of beta values (percent methylation), where rows are CpG sites and columns are samples.
#' @param threshold_delta_beta A numeric threshold (in range min(betavalues)–max(betavalues)) above which a CpG site is considered too variable and is removed.
#' @param doplot Logical. If `TRUE`, a histogram of delta beta values across groups is plotted.
#' @param groups A character vector of group labels to search within column names of `betavalue_df` (e.g., control sample identifiers).
#'
#' @return A character vector of CpG site row names that passed the filtering (i.e., all groups have delta beta < threshold).
#'
#' @details
#' - For each group, all columns in `betavalue_df` whose names contain the group label are extracted.
#' - For each CpG (row), the range of beta values across replicates in each group is computed.
#' - CpGs where any group has a delta beta ≥ `threshold_delta_beta` are filtered out.
#'
#' @examples
#' \dontrun{
#' selected_cpgs <- filter_by_beta_in_controls(
#'   betavalue_df = beta_matrix,
#'   threshold_delta_beta = 25,
#'   doplot = TRUE,
#'   groups = c("control1", "control2")
#' )
#' }
#'
#' @import xfun readxl dplyr tidyr
#' @export

filter_by_beta_in_cpgs<- function(betavalue_df, threshold_delta_beta = 25, doplot = TRUE, groups) {
  
  # Input validation 
  if (!is.data.frame(betavalue_df)) {
    stop("'betavalue_df' must be a data.frame.")
  }
  
  if (!is.numeric(threshold_delta_beta) || length(threshold_delta_beta) != 1 || is.na(threshold_delta_beta)) {
    stop("'threshold_delta_beta' must be a single numeric value.")
  }
  
  if (!is.logical(doplot) || length(doplot) != 1 || is.na(doplot)) {
    stop("'doplot' must be a single logical value (TRUE or FALSE).")
  }
  
  if (missing(groups) || !is.character(groups) || length(groups) == 0) {
    stop("'groups' must be a non-empty character vector.")
  }
  
  # Initialize a results data frame for delta beta values
  betavalue_grouped_sample <- data.frame(matrix(nrow = nrow(betavalue_df), ncol = length(groups)))
  colnames(betavalue_grouped_sample) <- groups
  
  # Compute delta beta (max - min) for each group
  for (grp in groups) {
    cols <- grep(grp, colnames(betavalue_df), value = TRUE)
    if (length(cols) < 2) {
      warning(paste("Group", grp, "has fewer than 2 matching columns — skipping."))
      next
    }
    betavalue_grouped_sample[[grp]] <- apply(
      betavalue_df[, cols, drop = FALSE],
      1,
      function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
    )
  }
  
  # Optional plot: histogram of delta beta values
  if (doplot) {
    vector <- betavalue_grouped_sample %>%
      pivot_longer(cols = everything(), values_to = "value") %>%
      dplyr::pull(value)
    
    do_histogram_with_statistics(
      vect = na.omit(vector),
      title = "DeltaBeta across replicate groups",
      x_axis = "DeltaBeta",
      text_position = 20,
      color_plot = "lightblue"
    )
  }
  
  # Filter CpGs by threshold
  rownames(betavalue_grouped_sample) <- rownames(betavalue_df)
  
  # Keep CpGs where all group delta beta values are < threshold
  keep_cpgs <- apply(betavalue_grouped_sample, 1, function(x)
    all(x < threshold_delta_beta, na.rm = TRUE)
  )
  
  betavalue_filtered <- betavalue_df[keep_cpgs, ]
  location_CpGs_beta_filtered <- rownames(betavalue_filtered)
  
  return(location_CpGs_beta_filtered)
}

#' Compute the Standard Deviation of Beta Values Across Regions and Filter Based on Threshold
#'
#' This function computes the standard deviation (SD) of Beta values for CpG regions across different groups of samples.
#' It also visualizes the distribution of SD values across the groups.
#'
#' @param regions_cpgs A data frame or data.table containing methylation data for CpG regions.
#'                      The data should include columns for methylation values across different experimental groups, 
#'                      and a column `region_id` identifying the CpG region and a `cpg_id` column that map the CpGs inside the region.
#' @param meth A character string specifying the type of methylation data. Currently, only "betavalues" is supported.
#'             Default is "betavalues".
#' @param groups A character vector of group labels to search within column names of `regions_cpgs` (e.g., control sample identifiers).
#'
#' @return A dataframe containing the `region_id` and the SD of the experimental groups.
#'
#' @details
#' The function computes the standard deviation (SD) for Beta values across different groups (replicates) for each CpG region. 
#' It provides histograms to visualize the SD distributions average across samples.
#' 
#' The function assumes that `regions_cpgs` includes a column `region_id` for CpG region identifiers and columns 
#' for methylation values, which are grouped by experimental conditions (e.g., different time points or sample types).
#'
#' @import xfun
#' @import readxl
#' @import dplyr
#' @import tidyr
#' @import stringr
#' @import data.table
#'
#' @export
compute_beta_sd_regions <- function(regions_cpgs, groups, meth = "betavalues") {
  # Load necessary libraries
  library(xfun)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(data.table)  # Fast data handling
  
  # Remove cpg_id column
  betavalue_df <- regions_cpgs
  betavalue_df$cpg_id <- NULL
  data = betavalue_df
  
  setDT(data)  # Convert data to a data.table for efficiency
  
  # Function to calculate the standard deviation for each group in each region
  calculate_group_sd <- function(region_data) {
    group_sd <- lapply(groups, function(grp) {
      # Select relevant columns for the group
      cols <- grep(paste0("^", grp), names(region_data), value = TRUE)
      
      # Check if there are enough columns to calculate SD
      if (length(cols) > 0) {
        # Compute SD, ignoring NA values
        sd_values <- sd(unlist(region_data[, ..cols]), na.rm = TRUE)
        return(sd_values)
      } else {
        warning(paste("No columns found for group:", grp))
        return(NA)  # Return NA if no columns are found for the group
      }
    })
    
    # Convert list to named list for correct assignment in data.table
    names(group_sd) <- groups
    return(as.list(group_sd))
  }
  
  # Apply the function to compute SD for each region_id
  sd_results <- data[, calculate_group_sd(.SD), by = region_id]
  
  # Convert SD results into a vector for visualization
  sd_results_vector <- sd_results
  sd_results_vector$region_id <- NULL
  sd_results_vector <- sd_results_vector %>%
    pivot_longer(cols = everything(), values_to = "value") %>%
    dplyr::pull(value)
  
  # Plot histogram of SD values
  do_histogram_with_statistics(
    vect = na.omit(sd_results_vector),
    title = "SD betavalue within a region for each sample",
    x_axis = "SD betavalue",
    text_position = 20,
    color_plot = "lightblue"
  )
  
  # Filter regions: remove those with at least one SD greater than `significance_beta`
  # rownames(sd_results) <- sd_results$region_id
  # sd_results <- sd_results[,-1]
  # filtered_regions <- sd_results[rowSums(sd_results > significance_beta) == 0,]
  # 
  # # Convert to data frame and retain row names (CpG regions)
  # filtered_regions <- as.data.frame(filtered_regions, row.names = rownames(filtered_regions))

  # Return the filtered region IDs
  return(sd_results)
}
