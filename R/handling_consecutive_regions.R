#' Investigate Adjacent Regions in Methylation Data
#'
#' This function identifies and investigates adjacent regions in methylation data. It maps CpG sites to regions, 
#' identifies consecutive regions, and computes the distances between adjacent regions. The function returns a summary 
#' of adjacent regions with optional plotting for distance computation.
#'
#' @param meth A data frame or data.table containing filtered methylation data. It should include 
#'                              methylation values for CpG sites and associated metadata.
#' @param tiles A data frame or vector containing the region boundaries or "tiles" that correspond to CpG sites in 
#'              the `meth` data.
#' @param ... Additional parameters passed to the `compute_regions_distance` function, such as `doplot` for optional 
#'            plotting of the region distances. This allows for flexibility in passing arguments to the distance 
#'            computation function.
#'
#' @return A list with two elements:
#' \itemize{
#'   \item `regions_cpgs_consecutive`: A data frame or data.table of consecutive CpG regions identified in the 
#'                                     methylation data.
#'   \item `pairs_region_summary`: A summary data frame or list containing the computed distances between adjacent 
#'                                regions, based on the provided tiles.
#' }
#'
#' @details
#' This function first maps CpG sites in the `filtered_meth_control` data to regions using the provided `tiles`. 
#' Then, it identifies consecutive regions and computes distances between adjacent regions using the `compute_regions_distance` 
#' function. The results are returned as a list containing:
#' - A `regions_cpgs_consecutive` data frame with identified consecutive regions.
#' - A `pairs_region_summary` summary containing the distance metrics between adjacent regions.
#' 
#' You can provide additional parameters to `compute_regions_distance` using the `...` argument (e.g., for plotting).
#' 
#' @importFrom magrittr %>%
#' @import xfun
#' @import dplyr
#' @import data.table
#'
#' @export
adjacent_regions_investigation <- function(meth, tiles) {
  
  # Map CpG to regions using the provided tile boundaries
  regions_cpgs <- map_cpg_to_regions(methdf = meth,
                                     region = tiles)
  
  # Identify consecutive regions
  regions_cpgs_consecutive <- identify_consecutive_regions(regions_cpgs)
  
  # Compute distances between consecutive regions, allow additional parameters via "..."
  pairs_region_summary <- compute_regions_distance(regions_cpgs_consecutive, doplot = FALSE)
  
  not_consecutive_regions <- unique(regions_cpgs$region_id[-which(regions_cpgs$region_id %in% regions_cpgs_consecutive$region_id)])
  message(paste0("The number of regions that will not be merged are ", unique(length(not_consecutive_regions))))
  
  
  # Return the results as a list
  result_adj <- list(regions_cpgs_consecutive = regions_cpgs_consecutive,
                     pairs_region_summary = pairs_region_summary,
                     not_consecutive_regions = not_consecutive_regions)
  
  return(result_adj)
}

#' Identify Consecutive Regions in Methylation Data
#'
#' This function identifies consecutive regions in methylation data by comparing the region boundaries and chromosomes. 
#' It determines whether the end of a region is immediately followed by the start of another region on the same chromosome, 
#' thereby identifying consecutive regions. The function also identifies gaps between consecutive CpGs and marks the 
#' first CpG of the following region after a consecutive sequence.
#'
#' @param regions_cpgs A data frame or data.table containing CpG regions with `region_id`, and the location of the CpGs inside the regions, defined with `chr`, `start`, and `end` 
#'                     columns. It represents the methylation data with regions of interest, where each region is defined 
#'                     by its boundaries and associated chromosome.
#' @param ... Additional parameters passed to the function (currently not used explicitly).
#'
#' @return A data.table with the original `regions_cpgs` data plus additional columns:
#' \itemize{
#'   \item `consecutive`: A logical column indicating whether a region is consecutive to the previous one.
#'   \item `prev_end`: The `end` position of the previous region, shifted for each row.
#'   \item `first_true_after_false`: A logical column marking the first CpG of a region following a consecutive gap.
#' }
#'
#' @details
#' This function processes the input `regions_cpgs` data and identifies consecutive regions based on the following rules:
#' - A region is considered consecutive if the end of one region is immediately followed by the start of the next region 
#'   on the same chromosome.
#' - The function then identifies gaps between CpGs from the last region to the first CpG of the consecutive region.
#' - It also marks the first CpG of a consecutive region after a gap as `first_true_after_false`.
#' 
#' The result is a `data.table` with updated information regarding consecutive regions, gap identification, and a logical 
#' indicator for each region's consecutive status.
#'
#' @import xfun
#' @import dplyr
#' @import data.table
#'
#' @export
identify_consecutive_regions <- function(regions_cpgs, ...) {
  
  # Determine the region locations
  regions_df <- data.table(region_id = unique(regions_cpgs$region_id))
  regions_df[, c("chr", "start", "end") := tstrsplit(region_id, "[:-]", type.convert = TRUE)]
  
  # Identify consecutive regions based on chromosome and start-end positions
  regions_df$consecutive <- c(
    FALSE,  # First row is always FALSE (no previous row to compare)
    (regions_df$chr[-1] == regions_df$chr[-nrow(regions_df)]) &  # Same chromosome?
      (regions_df$start[-1] == regions_df$end[-nrow(regions_df)] + 1)  # Start of next == End of previous + 1?
  )
  
  regions_cpgs$consecutive <- regions_df$consecutive[match(regions_cpgs$region_id, regions_df$region_id)]
  
  # Define the gap between CpGs at the boundary of consecutive regions
  regions_cpgs$prev_end <- c(NA, regions_cpgs$end[-nrow(regions_cpgs)])
  regions_cpgs$first_true_after_false <- c(FALSE, diff(as.integer(regions_cpgs$consecutive)) == 1)  # Gap marker
  
  # Identify gap regions (function called inside)
  regions_cpgs <- identify_gap_regions(regions_cpgs)
  
  #Let's see the distribution of consecutive regions
  consecutive_regions <- number_of_consecutive_regions(regions_cpgs)
  
  #in this way we identify only the regions that can be merged; all the ones that are not in here are considered notconsecutive afterwards
  tomerge_regions <- select_couples_to_compare(consecutive_regions)
  
  message(paste0("The number of consecutive regions that will be merged are ", length(unique(tomerge_regions$region_id))))
    
  return(tomerge_regions)
}

#' Identify and Plot the Number of Consecutive Regions
#'
#' This function identifies stretches of consecutive regions (based on a logical \code{consecutive} column)
#' in the provided data frame and visualizes their frequency as a bar plot.
#' It assigns each group a count and returns the input data frame with an added column \code{groupN}.
#'
#' @param regions_cpgs A data frame containing at least the columns \code{region_id} and \code{consecutive}.
#'                     The \code{consecutive} column must be logical (\code{TRUE}/\code{FALSE}).
#'
#' @return A data frame identical to the input but with an additional column \code{groupN},
#'         indicating the size of the consecutive region group to which each row belongs.
#'         Also prints a bar plot showing the frequency of group sizes.
#'         
#' @import ggplot2
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' df <- data.frame(region_id = paste0("region_", 1:10), consecutive = c(FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE))
#' number_of_consecutive_regions(df)
#' }
number_of_consecutive_regions <- function(regions_cpgs) {
  
  # Input checks
  if (!is.data.frame(regions_cpgs)) stop("Input must be a data frame.")
  if (!all(c("region_id", "consecutive") %in% colnames(regions_cpgs))) {
    stop("Input must contain 'region_id' and 'consecutive' columns.")
  }
  if (!is.logical(regions_cpgs$consecutive)) {
    stop("'consecutive' column must be of type logical (TRUE/FALSE).")
  }
  
  # Keep only unique region_id-consecutive pairs to avoid duplicate group counting
  unique_regions <- unique(regions_cpgs[, c("region_id", "consecutive")])
  
  # Count how many TRUEs follow each FALSE in the 'consecutive' vector
  result <- count_true_following_false(unique_regions$consecutive)
  unique_regions$groupN <- result  # Assign group count to the unique regions
  
  # Build frequency table of group sizes, filtering only FALSE positions
  result_df <- data.frame(original = unique_regions$consecutive, count_true_after_false = result)
  result_df <- na.omit(result_df[result_df$original == FALSE, ])  # focus on FALSE positions
  freq_table <- as.data.frame(table(result_df$count_true_after_false))  # tabulate group sizes
  
  # Create the plot
  p <- ggplot(freq_table, aes(x = Var1, y = Freq)) +
    geom_col(fill = "steelblue") +
    labs(
      title = "Frequency of Consecutive Region Pairs",
      x = "Consecutive Pairs",
      y = "Counts"
    ) +
    theme_minimal()
  print(p)
  
  # Assign groupN values back to the full dataset based on region_id
  regions_cpgs$groupN <- unique_regions$groupN[match(regions_cpgs$region_id, unique_regions$region_id)]
  
  return(regions_cpgs)
}

#' Count Consecutive TRUE Values Following Each FALSE in a Logical Vector
#'
#' This function scans a logical vector and, for each `FALSE`, counts how many `TRUE` values
#' immediately follow it. The count is also assigned to the following `TRUE` values in that group.
#' This is useful for identifying stretches of consecutive `TRUE` values after a `FALSE`.
#'
#' @param vec A logical vector (containing only `TRUE` or `FALSE` values).
#'
#' @return An integer vector of the same length as `vec`, where `FALSE` positions are assigned
#'         the number of consecutive `TRUE`s that follow, and those `TRUE`s receive the same count.
#'         Positions that are not part of a counted group are filled with `NA`.
#'
#' @examples
#' count_true_following_false(c(FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE))
#'
#' @export
count_true_following_false <- function(vec) {
  
  # Input checks
  if (!is.logical(vec)) {
    stop("Input must be a logical vector (TRUE/FALSE).")
  }
  
  n <- length(vec)
  result <- rep(NA_integer_, n)  # Initialize output vector with NAs
  
  for (i in seq_len(n)) {
    if (!vec[i]) {  # If current value is FALSE
      count_true <- 0
      j <- i + 1
      # Count how many TRUEs follow immediately
      while (j <= n && vec[j]) {
        count_true <- count_true + 1
        j <- j + 1
      }
      # Assign count to FALSE and its following TRUEs (if any)
      if (count_true > 0) {
        result[i] <- count_true
        if (count_true == 1) {
          result[i + 1] <- count_true
        } else {
          result[(i + 1):(i + count_true)] <- count_true
        }
      } else {
        result[i] <- 0  # No TRUEs follow this FALSE
      }
    }
  }
  
  return(result)
}


#' Select Consecutive Region Pairs for Merging
#'
#' This function identifies region pairs to be merged from a data frame containing consecutive CpG regions.
#' It avoids reusing the same region more than once, and ensures consistent pairing across groups of consecutive regions.
#' Regions without a consecutive group (`groupN == 0`) are excluded. In groups with an odd number of regions,
#' the last unmatched region is skipped.
#'
#' @param consecutive_regions A data frame containing at least the columns \code{region_id}, \code{consecutive}, and \code{groupN}.
#'                            Each row corresponds to a genomic region with information about whether it is consecutive and its group.
#'
#' @return A filtered data frame of regions selected to be merged, based on consecutive group logic.
#'
#' @details
#' - Regions are grouped by \code{groupN}, with \code{groupN == 0} representing non-consecutive regions, which are removed.
#' - Within each group, regions are paired such that no region is reused.
#' - In groups with an even group number, every \code{(group index + 1)}\emph{th} region is skipped to avoid overlapping pairs.
#'
#' @examples
#' \dontrun{
#' selected <- select_couples_to_compare(regions_cpgs_with_groups)
#' }
#'
#' @export
select_couples_to_compare <- function(consecutive_regions) {
  
  # Input validation
  required_cols <- c("region_id", "consecutive", "groupN")
  missing_cols <- dplyr::setdiff(required_cols, colnames(consecutive_regions))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Split regions by group
  splitted_consecutive_groups <- split(consecutive_regions, consecutive_regions$groupN)
  
  # Remove non-consecutive group (groupN == 0)
  splitted_consecutive_groups[["0"]] <- NULL
  
  # Keep only unique region IDs per group
  splitted_consecutive_groups_reg <- lapply(splitted_consecutive_groups, function(df) {
    unique(df[, c("region_id", "consecutive")])
  })
  
  # Filter each group to avoid reusing regions
  processed_dfs <- lapply(seq_along(splitted_consecutive_groups_reg), function(nlist) {
    group_name <- names(splitted_consecutive_groups_reg)[[nlist]]
    df <- splitted_consecutive_groups_reg[[group_name]]
    n <- nrow(df)
    idx <- as.numeric(group_name)
    
    # For odd-numbered group indices, keep all regions
    if (idx %% 2 == 1) {
      return(df)
    }
    
    # For even-numbered group indices, skip every (idx + 1)th row
    skip_every <- idx + 1
    row_nums <- seq_len(n)
    keep_rows <- row_nums %% skip_every != 0
    
    df[keep_rows, , drop = FALSE]
  })
  
  # Merge all processed region subsets
  tomerge_regions <- do.call(rbind, processed_dfs)
  
  # Filter original input to include only selected regions
  consecutive_regions_tomerge <- consecutive_regions[consecutive_regions$region_id %in% tomerge_regions$region_id, ]
  
  return(consecutive_regions_tomerge)
}

#' Identify Gap Regions Between Consecutive CpGs
#'
#' This function identifies the gap between consecutive regions in methylation data. Specifically, it computes the 
#' distance between the end of one region and the start of the next region that is consecutively attached to it. 
#' The gap is calculated only for CpGs that mark the boundary between two consecutive regions.
#'
#' @param regions_cpgs A data frame or data.table containing the CpG regions data with columns `start`, `end`, and 
#'                     `first_true_after_false`. This data represents the methylation regions along with boundary 
#'                     information.
#' @param doplot A logical value. If `TRUE`, a histogram will be generated showing the distribution of the gap values 
#'               between consecutive regions. Default is `FALSE`.
#'
#' @return A data.table with an additional column:
#' \itemize{
#'   \item `gap`: A numeric column representing the distance (in CpGs) between the end of one region and the start 
#'               of the next consecutive region. `NA` is assigned for rows that do not mark a gap between consecutive regions.
#' }
#'
#' @details
#' This function calculates the gap between consecutive CpG regions, which is defined as the difference between the 
#' `start` position of the first CpG in the consecutive region and the `end` position of the last CpG in previous region. It applies the calculation only for CpGs 
#' where the `first_true_after_false` flag is `TRUE`, indicating that the region marks the beginning of a new consecutive region.
#' If the `doplot` argument is set to `TRUE`, the function will generate a histogram to visualize the distribution of gap values.
#' 
#' @import xfun
#' @import dplyr
#' @import data.table
#'
#' @export
identify_gap_regions <- function(regions_cpgs, doplot = FALSE) {
  
  # Compute the gap between consecutive regions
  regions_cpgs$gap <- ifelse(regions_cpgs$first_true_after_false, 
                             regions_cpgs$start - regions_cpgs$prev_end, 
                             NA)
  
  # Generate histogram plot if doplot is TRUE
  if (doplot) {
    do_histogram_with_statistics(
      vect = na.omit(regions_cpgs$gap),
      title = "Distance (CpGs) between two regions that are attached to each other",
      x_axis = "gap",
      text_position = 80,
      color_plot = "khaki3"
    )
  }
  
  return(regions_cpgs)
}

#' Compute the Distance Between Consecutive Regions
#'
#' This function computes the distance between two consecutive regions based on their boundary CpG locations.
#' Specifically, it calculates the distance from the last CpG of the first region to the first CpG of the second region.
#' Additionally, it checks for any inconsistencies, such as when the end of a region is before its start, and raises a warning.
#'
#' @param regions_cpgs_consecutive A data.table or data frame containing the CpG regions that are consecutive, with information 
#'                                about the region ( `region_id`, `consecutive` flag) and CpG location inside the region `start`, `end`.
#' @param doplot A logical value indicating whether a histogram of the computed distances should be plotted. Default is `TRUE`.
#'
#' @return A data.table containing:
#' \itemize{
#'   \item `region_pairs`: A string representing the pair of consecutive regions.
#'   \item `A`: The position of the first CpG in the first region.
#'   \item `B`: The position of the last CpG in the first region.
#'   \item `C`: The position of the first CpG in the second region.
#'   \item `D`: The position of the last CpG in the second region.
#'   \item `DA`: The distance (in base pairs) between the last CpG of the first region and the first CpG of the second region.
#'   \item `CB`: The distance (in base pairs) between the first CpG of the second region and the last CpG of the first region.
#' }
#'
#' @details
#' The function computes the distances between consecutive regions where the `consecutive` flag is `TRUE`. 
#' It first checks whether the end position of a region appears before its start, raising a warning if so (indicating the 
#' regions may come from different strands). It then computes the distances between consecutive regions using the boundary 
#' CpGs of each region (i.e., positions A, B, C, D). A histogram of the distance between consecutive regions (D-A) is 
#' optionally generated based on the `doplot` argument.
#'
#' @import data.table
#' @import dplyr
#'
#' @export
compute_regions_distance <- function(regions_cpgs_consecutive, doplot = TRUE) {
  
  # Input validation
  
  required_cols <- c("region_id", "start", "end", "consecutive")
  missing_cols <- dplyr::setdiff(required_cols, names(regions_cpgs_consecutive))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required column(s):", paste(missing_cols, collapse = ", ")))
  }
  
  # Check column types
  if (!is.numeric(regions_cpgs_consecutive$start)) {
    stop("'start' column must be numeric.")
  }
  if (!is.numeric(regions_cpgs_consecutive$end)) {
    stop("'end' column must be numeric.")
  }
  if (!is.logical(regions_cpgs_consecutive$consecutive)) {
    stop("'consecutive' column must be logical (TRUE/FALSE).")
  }
  
  # Find A, B points (first and last CpG locations for each region)
  region_summary <- regions_cpgs_consecutive %>%
    group_by(region_id) %>%
    summarise(
      A = min(end),   # First CpG location in the region
      B = max(end)    # Last CpG location in the region
    ) %>%
    ungroup()
  
  
  # Reorder region_summary to match the order in regions_cpgs_consecutive
  region_summary <- region_summary[match(regions_cpgs_consecutive$region_id, region_summary$region_id),]
  region_summary <- unique(region_summary) #since there is one row per each cpg and we want one raw per each region
  
  # Add consecutive information
  region_summary$consecutive <- regions_cpgs_consecutive$consecutive[match(region_summary$region_id, regions_cpgs_consecutive$region_id)]
  
  # Convert to data.table for efficiency
  library(data.table)
  region_summary <- data.table(region_summary)
  
  n <- nrow(region_summary)
  
  # Take odd and even indexes
  idx1 <- seq(1, n - 1, by = 2)  # 1, 3, 5, ...
  idx2 <- idx1 + 1              # 2, 4, 6, ...
  
  # Create pairs table
  pairs_region_summary <- data.table(
    region_pairs = paste(region_summary$region_id[idx1], region_summary$region_id[idx2], sep = ";"),
    A = region_summary$A[idx1],
    B = region_summary$B[idx1],
    C = region_summary$A[idx2],
    D = region_summary$B[idx2],
    DA = region_summary$B[idx2] - region_summary$A[idx1],
    CB = region_summary$A[idx2] - region_summary$B[idx1]
  )
  
  # Generate histogram if doplot is TRUE
  if (doplot) {
    do_histogram_with_statistics(
      vect = as.numeric(pairs_region_summary$DA),
      title = "Distance between two consecutive regions \nD-A : Distance LAST CpG SECOND region - FIRST CpG FIRST region",
      x_axis = "Delta base pair",
      text_position = 500,
      color_plot = "orange"
    )
  }
  
  return(pairs_region_summary)
}



#' Merge Consecutive and Non-Consecutive Regions Based on Distance Thresholds
#'
#' This function identifies and merges consecutive regions based on distance thresholds for the `D-A` and `C-B` values. 
#' It also handles non-consecutive regions by adding them to the final output, providing a combined data set of all regions, both merged and non-merged.
#' The merging process involves analyzing the distances between regions and applying the specified thresholds to determine which regions should be merged.
#'
#' @param filtered_meth_control A data frame or data table containing methylation data for control samples, which will be used to investigate adjacent regions.
#' @param tiles A data frame or data table containing tile information to map CpGs to regions.
#' @param da_thr A numeric value representing the threshold for the `D-A` distance. Regions with a `D-A` distance below this threshold will be merged.
#' @param cb_thr A numeric value representing the threshold for the `C-B` distance. Regions with a `C-B` distance below this threshold will be merged.
#' @param ... Additional arguments that can be passed to the internal functions.
#'
#' @return A data frame containing the combined information for merged regions and non-consecutive regions.
#'         This includes columns such as chromosome (`chr`), start position (`start`), end position (`end`), and region ID (`region_id`).
#'         The final dataset contains both merged regions and non-consecutive regions.
#'
#' @details
#' The function first investigates the adjacent regions in the data using the `adjacent_regions_investigation` function. 
#' It then splits the data based on the `D-A` and `C-B` distance thresholds, visualizes the results, and proceeds with the merging process. 
#' Non-consecutive regions that do not meet the merging criteria are kept in the final dataset.
#'
#' @import data.table
#' @import dplyr
#'
#' @export
merging_consecutive_regions <- function(filtered_meth_control, tiles, da_thr, cb_thr, ...) {
  # Investigate adjacent regions
  consecutive_regions <- adjacent_regions_investigation(filtered_meth_control, tiles)
  
  # Split data based on distance thresholds for D-A and C-B
  splitted <- split_data_regions_distances(df = consecutive_regions$pairs_region_summary, da_thr, cb_thr)
  
  # Plot counts based on the distance thresholds
  p <- plot_counts(splitted, da_thr, cb_thr)
  print(p)
  
  # Merge regions based on the distance thresholds
  merged_regions <- merging_regions(pairs_region_summary = consecutive_regions$pairs_region_summary, 
                                    regions_cpgs = consecutive_regions$regions_cpgs_consecutive, 
                                    da_thr, cb_thr)
  
  # Extract and prepare the non-consecutive regions for merging
  
  not_consecutive_regions <- data.frame(region_id = consecutive_regions$not_consecutive_regions)
  
  not_consecutive_regions <- data.frame(
    chr = sub("^chr([^:]+):.*", "\\1", not_consecutive_regions$region_id),
    start = as.integer(sub("^[^:]+:(\\d+)-.*", "\\1", not_consecutive_regions$region_id)),
    end = as.integer(sub("^.*-(\\d+)$", "\\1", not_consecutive_regions$region_id)),
    region_id = not_consecutive_regions$region_id
  )
  
  # Ensure uniqueness of non-consecutive regions
  not_consecutive_regions <- unique(not_consecutive_regions)
  
  # Combine merged and non-consecutive regions
  all_regions <- rbind(merged_regions, not_consecutive_regions)
  
  return(all_regions)
}

#' Merge Regions Based on Distance Thresholds and Optimal Density
#'
#' This function merges regions based on the distances between regions, applying the specified `D-A` and `C-B` distance thresholds. 
#' It handles different types of merging strategies based on the threshold values and processes regions to either merge them or keep them as separate regions.
#'
#' @param pairs_region_summary A data frame containing summaries of region pairs with distance information between consecutive regions.
#' @param regions_cpgs A data frame or data table containing CpG region information with corresponding coordinates and region IDs.
#' @param da_thr A numeric value representing the threshold for the `D-A` distance. Regions with a `D-A` distance below this threshold will be merged.
#' @param cb_thr A numeric value representing the threshold for the `C-B` distance. Regions with a `C-B` distance below this threshold will be merged.
#'
#' @return A data frame containing processed region information, including columns for chromosome (`chr`), start position (`start`), end position (`end`), and region ID (`region_id`). 
#'         The regions in the returned data frame are either merged or retained as separate based on the specified thresholds.
#'
#' @details
#' The function first splits the region pairs based on the `D-A` and `C-B` distance thresholds using the `split_data_regions_distances` function. 
#' It then merges regions according to different strategies:
#' - `D-A` less than the threshold: Regions are pulled together based on the `da_thr` threshold.
#' - `C-B` less than the threshold: Regions are merged using optimal density considerations.
#' - `C-B` greater than the threshold: These regions are not merged and are kept as separate regions.
#' The function returns a combined data frame containing the processed regions.
#'
#' @import data.table
#' @import dplyr
#'
#' @export
merging_regions <- function(pairs_region_summary, regions_cpgs, da_thr, cb_thr) {

  # Split data based on distance thresholds for D-A and C-B
  splitted <- split_data_regions_distances(df = pairs_region_summary, da_thr, cb_thr)
  
  # Merge regions with D-A less than threshold
  df_da_less <- rbind(splitted$da_less_cb_less, splitted$da_less_cb_more)
  merged_less_da <- pull_regions_together(df_da_less, dimension_region = da_thr)
  
  # Merge regions with C-B less than threshold
  df_cb_less <- splitted$da_more_cb_less
  merged_less_cb <- merge_regions_optimal_density(df_cb_less, regions_cpgs, dimension_region = da_thr)
  
  # Handle regions with C-B greater than threshold (not merged)
  df_cb_more <- splitted$da_more_cb_more
  split_regions <- strsplit(df_cb_more$region_pairs, ";")
  
  # Extract first and second elements for non-merged regions
  df_cb_more$first <- sapply(split_regions, `[`, 1)
  df_cb_more$second <- sapply(split_regions, `[`, 2)
  notmerged_more_cb <- c(df_cb_more$first, df_cb_more$second)
  
  # Prepare non-merged regions
  notmerged_regions <- data.frame(
    chr = sub("^chr([^:]+):.*", "\\1", notmerged_more_cb),
    start = sub("^[^:]+:(\\d+)-.*", "\\1", notmerged_more_cb),
    end = sub("^.*-(\\d+)$", "\\1", notmerged_more_cb),
    region_id = notmerged_more_cb
  )
  
  # Combine merged and non-merged regions
  processed_regions <- rbind(merged_less_da, merged_less_cb, notmerged_regions)
  
  # Ensure correct data types for start and end positions
  processed_regions$start <- as.integer(processed_regions$start)
  processed_regions$end <- as.integer(processed_regions$end)
  
  return(processed_regions)
}

#' Split Data Based on Distance Thresholds for D-A and C-B
#'
#' This function splits the input data frame into four subsets based on the `D-A` and `C-B` distance thresholds. 
#' The data is divided into categories where the `D-A` distance is either less than or greater than the `da_thr` threshold, 
#' and within those, the `C-B` distance is further split based on the `cb_thr` threshold.
#'
#' @param df A data frame containing region distance information, including columns for `DA` (D-A distance) and `CB` (C-B distance).
#' @param da_thr A numeric value representing the threshold for the `D-A` distance. Data with `DA` less than or greater than this threshold will be split into separate categories.
#' @param cb_thr A numeric value representing the threshold for the `C-B` distance. Data with `CB` less than or greater than this threshold will be further split.
#'
#' @return A list containing four data frames:
#' \describe{
#'   \item{da_less_cb_less}{Subset of the data where `DA` is less than or equal to `da_thr` and `CB` is less than or equal to `cb_thr`.}
#'   \item{da_less_cb_more}{Subset of the data where `DA` is less than or equal to `da_thr` and `CB` is greater than `cb_thr`.}
#'   \item{da_more_cb_less}{Subset of the data where `DA` is greater than `da_thr` and `CB` is less than or equal to `cb_thr`.}
#'   \item{da_more_cb_more}{Subset of the data where `DA` is greater than `da_thr` and `CB` is greater than `cb_thr`.}
#' }
#'
#' @details
#' This function is useful for categorizing regions based on their `D-A` and `C-B` distances into four distinct groups, which can then be processed separately for further analysis or merging.
#'
#' @import dplyr
#'
#' @export
split_data_regions_distances <- function(df, da_thr, cb_thr) {
  
  if (!all(c("DA", "CB") %in% names(df))) {
    stop("Input data frame must contain 'DA' and 'CB' columns.")
  }
  if (!is.numeric(df$DA) || !is.numeric(df$CB)) {
    stop("'DA' and 'CB' columns must be numeric.")
  }
  
  # Filter da < da_threshold
  da_less <- df[df$DA <= da_thr, ]
  
  # divide da_less in two lists according to the cb_threshold
  da_less_cb_less <- da_less[da_less$CB <= cb_thr, ]
  da_less_cb_more <- da_less[da_less$CB > cb_thr, ]
  
  # Filter da >= da_threshold
  da_more <- df[df$DA > da_thr, ]
  
  # divide da_more in two lists according to the cb_threshold
  da_more_cb_less <- da_more[da_more$CB <= cb_thr, ]
  da_more_cb_more <- da_more[da_more$CB > cb_thr, ]
  
  # Return a list with 4 data frames
  return(list(
    da_less_cb_less = da_less_cb_less,
    da_less_cb_more = da_less_cb_more,
    da_more_cb_less = da_more_cb_less,
    da_more_cb_more = da_more_cb_more
  ))
}

#' Plot the Counts of Region Pairs by DA and CB Thresholds
#'
#' This function generates a bar plot showing the count of region pairs based on their `D-A` (DA) and `C-B` (CB) distances relative to the specified thresholds. The data is grouped into four conditions based on the comparison of `DA` and `CB` with the given thresholds.
#'
#' @param split_data_result A list containing four data frames, as returned by the `split_data_regions_distances` function. Each data frame represents a subset of region pairs based on the `DA` and `CB` distance thresholds.
#' @param da_thr A numeric value representing the threshold for the `D-A` distance.
#' @param cb_thr A numeric value representing the threshold for the `C-B` distance.
#'
#' @return A `ggplot` object representing the bar plot.
#'
#' @details
#' This function is useful for visualizing the distribution of region pairs across the four conditions defined by the `D-A` and `C-B` distance thresholds:
#' \describe{
#'   \item{DA <= da_thr & CB <= cb_thr}{Condition where both the `D-A` and `C-B` distances are below or equal to the respective thresholds.}
#'   \item{DA <= da_thr & CB > cb_thr}{Condition where the `D-A` distance is below or equal to the `da_thr`, but the `C-B` distance is greater than the `cb_thr`.}
#'   \item{DA > da_thr & CB <= cb_thr}{Condition where the `D-A` distance is greater than the `da_thr`, but the `C-B` distance is below or equal to the `cb_thr`.}
#'   \item{DA > da_thr & CB > cb_thr}{Condition where both the `D-A` and `C-B` distances are greater than the respective thresholds.}
#' }
#'
#' The resulting plot provides a clear visualization of how the data is distributed based on these thresholds.
#'
#' @import ggplot2
#'
#' @export
plot_counts <- function(split_data_result, da_thr, cb_thr) {
  
  # Input validation
  required_names <- c("da_less_cb_less", "da_less_cb_more", "da_more_cb_less", "da_more_cb_more")
  
  if (!is.list(split_data_result) || !all(required_names %in% names(split_data_result))) {
    stop("must contain exactly four named elements")
  }
  
  if (!all(sapply(split_data_result, is.data.frame))) {
    stop("Each element of split_data_result must be a data frame.")
  }
  
  if (!is.numeric(da_thr) || !is.numeric(cb_thr)) {
    stop("Both da_thr and cb_thr must be numeric.")
  }
  
  count_data <- data.frame(
    condition = c(paste0("DA <=",da_thr," & CB <=",cb_thr),
                  paste0("DA <=",da_thr," & CB >",cb_thr),
                  paste0("DA >",da_thr," & CB <=",cb_thr),
                  paste0("DA >",da_thr," & CB >",cb_thr)),
    count = c(
      nrow(split_data_result$da_less_cb_less),
      nrow(split_data_result$da_less_cb_more),
      nrow(split_data_result$da_more_cb_less),
      nrow(split_data_result$da_more_cb_more)
    )
  )
  
  # Create the bar plot
  ggplot(count_data, aes(x = condition, y = count)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    labs(title = "Count of Regions Pairs by DA and CB Thresholds", 
         x = "Conditions", 
         y = "Count of Regions") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12))
}

#' Merge Regions Based on Distance
#'
#' This function merges two adjacent regions based on the distance specified by the `dimension_region` parameter. The function takes a data frame containing region pairs and merges them by adjusting the starting and ending positions to create new regions that span the distance from the first CpG in the pair to the end of the region defined by `dimension_region`.
#'
#' @param regions_pairs_ofi A data frame containing pairs of regions. Each pair is represented as a column `region_pairs` which contains a string in the format "region1;region2". Additional columns `A` and `B` represent the start and end positions of the first region.
#' @param dimension_region A numeric value specifying the distance to extend the new merged region, starting from the first CpG position (`A`) of the region.
#'
#' @return A data frame containing the merged regions with columns:
#' \describe{
#'   \item{chr}{Chromosome of the merged region.}
#'   \item{start}{Start position of the merged region.}
#'   \item{end}{End position of the merged region, calculated as `start + dimension_region`.}
#'   \item{region_id}{A unique identifier for the merged region, formatted as "chr:start-end".}
#' }
#'
#' @details
#' This function is used to combine two regions into one based on their proximity and a specified distance threshold. The merging occurs by starting at the first CpG position (`A`) of the first region and extending the region to a new end position defined by `dimension_region`. This is useful when adjacent regions are to be considered as a single larger region for further analysis.
#'
#' @examples
#' # Example data frame with region pairs
#' regions_pairs_ofi <- data.frame(region_pairs = c("chr1:100-200;chr1:201-300", "chr2:500-600;chr2:601-700"),
#'                                 A = c(100, 500),
#'                                 B = c(200, 600))
#' # Merge the regions with a specified dimension
#' merged_regions <- pull_regions_together(regions_pairs_ofi, dimension_region = 50)
#'
#' @export
pull_regions_together <- function(regions_pairs_ofi, dimension_region){
  # Input checks
  if (!is.data.frame(regions_pairs_ofi)) stop("Input must be a data frame")
  if (!"region_pairs" %in% colnames(regions_pairs_ofi)) stop("Missing 'region_pairs' column")
  if (!"A" %in% colnames(regions_pairs_ofi)) stop("Missing 'A' column")
  if (!is.numeric(regions_pairs_ofi$A)) stop("'A' column must be numeric")
  if (!is.numeric(dimension_region) || length(dimension_region) != 1) stop("'dimension_region' must be a single numeric value")
  
  # dimension_region is the distance from the first CpG until the end of the region
  
  # Split each row correctly
  split_regions <- strsplit(regions_pairs_ofi$region_pairs, ";")
  # Extract first and second elements
  regions_pairs_ofi$first <- sapply(split_regions, `[`, 1)
  regions_pairs_ofi$second <- sapply(split_regions, `[`, 2)
  
  # Now we need to merge the regions together
  regions_pairs_ofi$chr <- sub("^chr([^:]+):.*", "\\1", regions_pairs_ofi$region_pairs)
  regions_pairs_ofi$newregion <- paste0("chr", regions_pairs_ofi$chr, ":", regions_pairs_ofi$A, "-", regions_pairs_ofi$A + dimension_region-1)
  
  # Create the new data frame with merged regions
  merged_regions <- data.frame(
    chr = regions_pairs_ofi$chr,
    start = regions_pairs_ofi$A,
    end = regions_pairs_ofi$A + dimension_region-1,
    region_id = paste0("chr", regions_pairs_ofi$chr, ":", regions_pairs_ofi$A, "-", regions_pairs_ofi$A + dimension_region-1)
  )
  
  return(merged_regions)
}

#' Merge Regions Based on Optimal CpG Density
#'
#' This function merges adjacent regions by computing the CpG density within each region pair. The merge is performed in such a way that the resulting merged regions have the highest CpG density within a specified window (`dimension_region`).
#'
#' @param regions_pairs_ofi A data frame containing pairs of regions, represented as strings in the `region_pairs` column (e.g., "region1;region2"). Additional columns `first` and `second` represent the start positions of the first region in each pair.
#' @param regions_cpgs A data frame containing CpG information, with a column `region_id` corresponding to the region identifiers from `regions_pairs_ofi`.
#' @param dimension_region A numeric value specifying the window size (in base pairs) for computing CpG density within each region.
#'
#' @return A data frame containing the merged regions with the highest CpG density, with the following columns:
#' \describe{
#'   \item{chr}{Chromosome of the merged region.}
#'   \item{start}{Start position of the merged region.}
#'   \item{end}{End position of the merged region, calculated as the start position plus `dimension_region`.}
#'   \item{region_id}{A unique identifier for the merged region, formatted as "chr:start-end".}
#' }
#'
#' @details
#' This function merges two regions from a region pair based on their CpG density. First, it extracts the CpG data for both regions in each pair and computes the CpG density over a window defined by `dimension_region`. Then, it merges the regions with the highest CpG density and returns the resulting merged regions.
#'
#' The function works by splitting the region pairs and calculating the CpG density for each pair. It then selects the regions with the highest density and merges them into a single region.
#'
#' @examples
#' # Example data frame with region pairs
#' regions_pairs_ofi <- data.frame(region_pairs = c("chr1:100-200;chr1:201-300", "chr2:500-600;chr2:601-700"),
#'                                 first = c("chr1:100-200", "chr2:500-600"),
#'                                 second = c("chr1:201-300", "chr2:601-700"))
#' 
#' # Example CpG data
#' regions_cpgs <- data.frame(region_id = c("chr1:100-200", "chr1:201-300", "chr2:500-600", "chr2:601-700"),
#'                            cpg_count = c(50, 60, 30, 40), 
#'                            stringsAsFactors = FALSE)
#' 
#' # Merge regions with optimal CpG density
#' merged_regions <- merge_regions_optimal_density(regions_pairs_ofi, regions_cpgs, dimension_region = 300)
#'
#' @export
merge_regions_optimal_density <- function(regions_pairs_ofi, regions_cpgs, dimension_region){
  # Split each row correctly
  split_regions <- strsplit(regions_pairs_ofi$region_pairs, ";")
  # Extract first and second elements
  regions_pairs_ofi$first <- sapply(split_regions, `[`, 1)
  regions_pairs_ofi$second <- sapply(split_regions, `[`, 2)
  
  roi <- c(regions_pairs_ofi$first, regions_pairs_ofi$second)
  regions_cpgs_roi <- regions_cpgs[regions_cpgs$region_id %in% roi, ]
  
  # Combine two consecutive regions in a list to create a region pair
  if (length(regions_pairs_ofi$first) == length(regions_pairs_ofi$second)) {
    group_regions <- mapply(c, regions_pairs_ofi$first, regions_pairs_ofi$second, SIMPLIFY = FALSE)
  }
  
  # Get all the CpGs belonging to that region pair
  group_regions_cpgs <- lapply(group_regions, function(x){
    df <- regions_cpgs_roi[regions_cpgs_roi$region_id %in% unlist(x),]
  })
  
  # Compute CpG density for each region pair
  cpg_density_list <- compute_cpg_density(df_list = group_regions_cpgs, window_size = dimension_region)
  
  # Compute the merged regions based on the highest CpG density
  merged_regions <- compute_optimal_density(cpg_density_list)
  
  return(merged_regions)
}

#' Compute CpG Density for Regions
#'
#' This function calculates the CpG density for regions within a specified window size. The CpG density is computed for each region by counting the number of CpGs within the defined window.
#'
#' @param df_list A data frame or a list of data frames, where each data frame contains CpG information. Each data frame must have columns `chr`, `start`, and `end`, representing the chromosome and the start and end positions of the CpGs.
#' @param window_size A numeric value indicating the window size (in base pairs) for computing CpG density. Default is 300.
#'
#' @return A data frame (or a list of data frames) containing the CpG density for each window, with the following columns:
#' \describe{
#'   \item{chr}{Chromosome of the CpG region.}
#'   \item{start}{Start position of the window.}
#'   \item{end}{End position of the window, calculated as `start + window_size`.}
#'   \item{num_cpgs}{The number of CpGs within the window.}
#' }
#'
#' @details
#' This function calculates the CpG density within a specified window size for each CpG in a given data frame. The function can either be applied to a single data frame or to a list of data frames. The results are sorted in descending order by the number of CpGs in each window.
#' 
#' The function works by iterating through each CpG in the input data and defining a window around each CpG. It then counts how many CpGs fall within the defined window and stores the results.
#'
#' @examples
#' # Example data frame with CpG data
#' cpg_data <- data.frame(chr = rep(1, 5),
#'                        start = c(100, 200, 300, 400, 500),
#'                        end = c(100, 200, 300, 400, 500),
#'                        stringsAsFactors = FALSE)
#' 
#' # Compute CpG density with a window size of 300
#' cpg_density <- compute_cpg_density(cpg_data, window_size = 300)
#'
#' @export
compute_cpg_density <- function(df_list, window_size = 300) {
  # Helper function to compute CpG density for a single dataframe
  compute_cpg_density_single <- function(df, window_size) {
    compute_cpgs_density <- function(df, start_idx, window_size) {
      # Extract chromosome and start position of the CpG
      chr <- df$chr[start_idx]
      start_cpg <- df$start[start_idx]
      
      # Define the window
      window_start <- start_cpg
      window_end <- start_cpg + window_size-1
      
      # Find CpGs within the window
      cpgs_in_window <- df[df$start >= window_start & df$end <= window_end, ]
      
      # Return dataframe with CpG density details
      return(data.frame(chr = chr, 
                        start = window_start, 
                        end = window_end, 
                        num_cpgs = nrow(cpgs_in_window)))
    }
    
    # Apply function to all CpGs in the dataframe
    cpgs_density <- lapply(1:nrow(df), function(i) compute_cpgs_density(df, i, window_size))
    
    # Combine results into a single dataframe and sort by number of CpGs
    cpgs_density_df <- do.call(rbind, cpgs_density)
    cpgs_density_df <- cpgs_density_df[order(cpgs_density_df$num_cpgs, decreasing = TRUE), ]
    
    return(cpgs_density_df)
  }
  
  # If input is a single dataframe, apply function directly
  if (is.data.frame(df_list)) {
    return(compute_cpg_density_single(df_list, window_size))
  }
  
  # If input is a list, apply function to each dataframe in the list
  if (is.list(df_list)) {
    return(lapply(df_list, function(df) compute_cpg_density_single(df, window_size)))
  }
  
  stop("Input must be a dataframe or a list of dataframes")
}

#' Compute Optimal CpG Density for Region Pairs
#'
#' This function processes a list of data frames containing CpG densities, identifies the regions with the highest CpG density, and returns the optimal regions based on those densities.
#'
#' The function distinguishes between regions with a unique maximum CpG density and those with multiple regions sharing the same maximum density. The user can modify the approach for handling multiple maxima in the future, but currently, the first maximum is chosen.
#'
#' @param cpg_density_list A list of data frames, where each data frame contains CpG density information with columns `chr`, `start`, `end`, and `num_cpgs`. Each data frame represents the density information for one region pair.
#'
#' @return A data frame with the optimal regions based on the maximum CpG density. The data frame contains the following columns:
#' \describe{
#'   \item{chr}{Chromosome of the region.}
#'   \item{start}{Start position of the region.}
#'   \item{end}{End position of the region.}
#'   \item{region_id}{A unique identifier for the region, formatted as `chr:start-end`.}
#' }
#'
#' @details
#' The function separates regions with a unique maximum CpG density from those where multiple regions share the same maximum density. For regions with a unique maximum, the region with the highest CpG count is returned. For regions with multiple maxima, the first region with the maximum CpG count is returned (this approach can be modified in the future).
#'
#' After processing, the function combines the results from both types of maxima and returns a data frame containing the optimal regions with the highest CpG density.
#'
#' @examples
#' # Example of CpG density data
#' cpg_density_example <- list(
#'   data.frame(chr = rep("chr1", 3), start = c(100, 200, 300), end = c(100, 200, 300), num_cpgs = c(5, 10, 15)),
#'   data.frame(chr = rep("chr1", 2), start = c(400, 500), end = c(400, 500), num_cpgs = c(20, 10))
#' )
#' 
#' # Compute the optimal regions based on CpG density
#' optimal_regions <- compute_optimal_density(cpg_density_example)
#' 
#' @export
compute_optimal_density <- function(cpg_density_list){
  # There is the possibility that some regions pair have same density for some regions. 
  # Thus we can split them and decide how to take the max;
  # I splitted the thing in two because in the future we can decide how to take the max differently
  
  # Lists to store results
  unique_max_density_list <- list()
  multiple_max_density_list <- list()
  
  # Process each dataframe in cpg_density_list
  for (i in seq_along(cpg_density_list)) {
    df <- as.data.frame(cpg_density_list[[i]])
    
    if (nrow(df) == 0) next  # Skip empty dataframes
    
    max_cpgs <- max(df$num_cpgs)  # Find max CpG count
    df_max <- df[which(df$num_cpgs == max_cpgs),]  # Filter rows with max CpG count
    
    # If there's only one unique max, store in unique_cpg_density_list
    if (nrow(df_max) == 1) {
      unique_max_density_list[[i]] <- as.data.frame(cpg_density_list[[i]])
    } else {
      # If multiple windows share the max, store in multiple_max_density_list
      multiple_max_density_list[[i]] <- as.data.frame(cpg_density_list[[i]])
    }
  }
  
  # Remove NULL elements from lists
  unique_max_density_list <- unique_max_density_list[lengths(unique_max_density_list) > 0]
  multiple_max_density_list <- multiple_max_density_list[lengths(multiple_max_density_list) > 0]
  
  # Process unique max density list
  if(length(unique_max_density_list) > 0){
    
    # For the ones that have a unique max
    unique_max_density_list_optimal <- list()
    for (i in seq_along(unique_max_density_list)) {
      df <- as.data.frame(unique_max_density_list[[i]])
      
      if (nrow(df) == 0) next  # Skip empty dataframes
      
      max_cpgs <- max(df$num_cpgs)  # Find max CpG count
      df_max <- df[which(df$num_cpgs == max_cpgs),]  # Filter rows with max CpG count
      
      unique_max_density_list_optimal[[i]] <- df_max
    }
    unique_max_density <- do.call(rbind, unique_max_density_list_optimal)
    unique_max_density$region_id <- paste0("chr",unique_max_density$chr, ":", unique_max_density$start, "-", unique_max_density$end)
    
  } else { unique_max_density <- NULL }
  
  # Process multiple max density list
  if(length(multiple_max_density_list) > 0){
    
    # First (location) more highly dense, this can be changed according to the plan
    multiple_max_density_list_optimal <- list()
    for (i in seq_along(multiple_max_density_list)) {
      df <- as.data.frame(multiple_max_density_list[[i]])
      
      if (nrow(df) == 0) next  # Skip empty dataframes
      
      max_cpgs <- max(df$num_cpgs)  # Find max CpG count
      df_max <- df[which(df$num_cpgs == max_cpgs),]  # Filter rows with max CpG count
      
      multiple_max_density_list_optimal[[i]] <- df_max[1,]  # Take the first row if multiple maxima
    }
    
    multiple_max_density <- do.call(rbind, multiple_max_density_list_optimal)
    multiple_max_density$region_id <- paste0("chr",multiple_max_density$chr, ":", multiple_max_density$start, "-", multiple_max_density$end)
    
  } else { multiple_max_density <- NULL }
  
  # Combine both unique and multiple max density results
  regions_optimal_density <- rbind(unique_max_density, multiple_max_density)
  regions_optimal_density$num_cpgs <- NULL  # Remove the `num_cpgs` column as it's not needed
  
  return(regions_optimal_density)
}



#' Assign Beta Values to Merged Regions
#'
#' This function assigns beta values to merged regions based on CpG locations. It associates each CpG in the provided beta value data with its corresponding region and returns the regions with their associated beta values.
#'
#' @param regions A data frame containing merged region information, including columns for `chr`, `start`, `end`, and `region_id`.
#' @param cpgs A data frame containing CpG data with columns for `chr`, `start`, `end`, and `cpg_id`, which are used to map CpGs to their respective regions.
#' @param betavalue A data frame or matrix containing beta values for CpGs. The row names should correspond to the `cpg_id` of the CpGs.
#'
#' @return A list containing two data frames:
#' \describe{
#'   \item{regions_cpgs}{A data frame containing merged regions with associated CpG information, including a `cpg_id` and `region_id`.}
#'   \item{regions_merged_cpg_beta}{A data frame containing the merged regions with their associated beta values. The columns include `region_id`, `cpg_id`, and beta values for each CpG.}
#' }
#'
#' @details
#' The function first maps CpGs to their respective regions by using the provided `regions` and `cpgs` data. It then merges the CpG data with the beta values based on `cpg_id` and associates each CpG with the appropriate region.
#'
#' The function outputs two data frames:
#' \itemize{
#'   \item The `regions_cpgs` data frame contains the merged regions with CpG information, including a `region_id`.
#'   \item The `regions_merged_cpg_beta` data frame includes the merged regions with the corresponding beta values for each CpG.
#' }
#'
#' @export
assign_beta_merged_regions <- function(regions, cpgs, betavalue){
  
  # Map CpGs to regions
  regions_cpgs <- map_cpg_to_regions(methdf = cpgs, region = regions)
  
  # Create a unique cpg_id column for regions_cpgs
  regions_cpgs$cpg_id <- paste0("chr", regions_cpgs$chr, ":", regions_cpgs$start, "-", regions_cpgs$end)
  
  # Assign cpg_id to betavalue
  betavalue$cpg_id <- rownames(betavalue)
  
  regions_merged_cpg_beta <- data.frame(region_id = regions_cpgs$region_id,
                                        cpg_id = regions_cpgs$cpg_id)
  
  # Merge betavalue with regions_cpgs by cpg_id
  regions_merged_cpg_beta2 <- betavalue[match(regions_merged_cpg_beta$cpg_id,betavalue$cpg_id),]
  regions_merged_cpg_beta2$cpg_id <- NULL
  regions_merged_cpg_beta <- cbind(regions_merged_cpg_beta,regions_merged_cpg_beta2)
  
  # Identify overlapping CpGs (where cpg_id is duplicated)
  #overlapping_cpg <- regions_cpgs[duplicated(regions_cpgs$cpg_id),]
  #regions_cpgs_overlaps <- regions_cpgs[which(regions_cpgs$cpg_id %in% overlapping_cpg$cpg_id),]
  
  # Update rownames to sequential numbers
  rownames(regions_merged_cpg_beta) <- 1:nrow(regions_merged_cpg_beta)
  
  # Return the results as a list
  results <- list(regions_cpgs = regions_cpgs, regions_merged_cpg_beta = regions_merged_cpg_beta)
  
  return(results)
}

#' Compute Mean or Median Beta Values for Merged CpG Regions
#'
#' This function summarizes beta values across CpG sites within each region by computing either the mean or median.
#' Replicate-level beta values are grouped using metadata, and summaries are computed per region for each defined group.
#'
#' @param regions_cpgs A \code{data.frame} or \code{data.table} containing beta values for CpG sites. Must include
#'   a \code{region_id} column, a \code{cpg_id} column, and one or more columns with beta values per replicate.
#' @param method A character string specifying the summary method. Either \code{"mean"} or \code{"median"}.
#' @param metadata A \code{data.frame} containing sample annotations. Must include a column matching the sample names
#'   used in \code{regions_cpgs}, and another column with group identifiers.
#' @param sample_name A character string indicating the name of the column in \code{metadata} that contains sample labels 
#'
#' @return A \code{data.table} summarizing the beta values for each \code{region_id} and each group. Each cell contains
#'   the mean or median beta value for that group in that region.
#'
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   region_id = c("region1", "region1", "region2"),
#'   cpg_id = c("cg001", "cg002", "cg003"),
#'   T0_24_1 = c(0.5, 0.6, 0.7),
#'   T0_24_2 = c(0.4, 0.7, 0.8),
#'   T0_48_1 = c(0.6, 0.5, 0.6)
#' )
#' compute_beta_merged_regions(df, method = "mean")
#' }
#'
#' @importFrom stringr str_extract
#' @importFrom data.table setDT .SD :=
#'
#' @export

compute_beta_merged_regions <- function(regions_cpgs, method, metadata, sample_name) {
  
  # Validate meth
  if (!is.character(method) || length(method) != 1 ||
      !method %in% c("mean", "median")) {
    stop("Error: method must be one of 'mean', 'median'")
  }
  
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
  setDT(betavalue_df)  # Convert betavalue_df to a data.table for efficiency
  
  # Extract unique groups (e.g., "T0_24_1", "T0_24_2")
  colnames_no_region <- dplyr::setdiff(names(betavalue_df), c("region_id"))
  
  if (!(sample_name %in% colnames(metadata))) {
    stop("Error: sample_name needs to be character string indicating the name of the column in metadata that contains sample labels ")
  }
  
  groups <- unique(metadata[[sample_name]])
  # Function to calculate the median for each group in each region
  calculate_group_means <- function(region_data) {
    group_means <- lapply(groups, function(grp) {
      # Select relevant columns for the group
      cols <- grep(paste0("^", grp), names(region_data), value = TRUE)
      
      # Check if there are enough columns to calculate the median
      if (length(cols) > 0) {
        # Calculate the median for each group, ignoring NA values
        if(method == "mean"){
          mean_values <- mean(unlist(region_data[, ..cols]), na.rm = TRUE)
        } else if(method == "median"){
          mean_values <- median(unlist(region_data[, ..cols]), na.rm = TRUE)
        }
        return(mean_values)
      } else {
        return(NA)  # Return NA if no columns are found for the group
      }
    })
    
    # Convert the list to a named list for correct assignment in data.table
    names(group_means) <- groups
    return(as.list(group_means))
  }
  
  # Apply the function to calculate means for each region_id
  mean_results <- betavalue_df[, calculate_group_means(.SD), by = region_id]
  
  # Return the result with the mean values for each group in each region
  return(mean_results)
}

#' Check for Overlaps Between CpG Regions Based on Shared CpGs
#'
#' This function identifies CpG sites (`cpg_id`) that appear in more than one CpG region (`region_id`), 
#' quantifies the extent of this overlap, and visualizes the number of overlapping CpGs per region.
#'
#' @param regions A data frame or data.table with at least two columns: `region_id` (identifier of each region) 
#'                and `cpg_id` (identifier of each CpG site). Each row corresponds to a CpG-region mapping.
#'
#' @return No value is returned. The function outputs messages to the console indicating:
#'   - How many CpGs are present in more than one region
#'   - How many regions contain overlapping CpGs
#' 
#' It also generates a histogram of the number of overlapping CpGs per region, to help assess whether overlaps 
#' need to be addressed.
#'
#' @details
#' Overlapping CpG regions can bias downstream analyses if not handled appropriately. This function provides a 
#' quick diagnostic by identifying and quantifying overlaps. Users can inspect the histogram to decide whether 
#' overlaps are frequent enough to require filtering or correction.
#'
#' The histogram is created using the `do_histogram_with_statistics()` function, which must be available in the environment.
#'
#' @import data.table
#'
#' @export

check_if_any_overlap_between_regions <- function(regions){
  
  if (!all(c("cpg_id", "region_id") %in% colnames(regions))) {
    stop("regions must contain columns: 'cpg_id' and 'region_id'")
  }
  
  library(data.table)
  duplicated_cpg <- unique(regions$cpg_id[duplicated(regions$cpg_id)])
  message(paste0("There are ",length(unique(duplicated_cpg)), " CpGs that are present in multiple regions"))
  regions_dupl <- regions[which(regions$cpg_id %in% duplicated_cpg),]
  length(unique(regions_dupl$region_id))
  message(paste0("There are ",length(unique(regions_dupl$region_id)), " overlapping regions"))
  
  distribution_regions_overlap <- as.data.frame(table(regions_dupl$region_id))
  
  if (length(unique(duplicated_cpg)) != 0) {
    #message("Give a look to the plot 'Number of overlapping CpGs per region' because if you have many CpGs overlapping you can decide to handle this, otherwise to ignore")
    
    # Load and execute the histogram function
    do_histogram_with_statistics(
      vect = as.numeric(distribution_regions_overlap$Freq),
      title = "Number of overlapping CpGs per region",
      x_axis = "Number regions",
      text_position = 25,
      color_plot = "orange"
    )
  }
  
  
  
}

