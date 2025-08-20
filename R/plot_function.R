#' Plot a Histogram with Summary Statistics Overlay
#'
#' This function plots a histogram of a numeric vector and overlays text labels with summary
#' statistics (Min, Q1, Mean, Median, Q3, Max). It can optionally save the plot to a file
#' in PNG, JPEG, or PDF format.
#'
#' @param vect A numeric vector to be plotted.
#' @param title A string for the plot title.
#' @param x_axis A string for the x-axis label.
#' @param text_position Numeric, the x-coordinate where summary statistics will be displayed.
#' @param color_plot Color for the bars and statistics text.
#' @param savefile Optional. A file path (ending in `.png`, `.jpg`, or `.pdf`) to save the plot.
#'
#' @return No return value. Displays (and optionally saves) a histogram with annotated statistics.
#' @export
#'
#' @examples
#' \dontrun{
#' do_histogram_with_statistics(
#'   vect = rnorm(100),
#'   title = "Random Values",
#'   x_axis = "Values",
#'   text_position = 0,
#'   color_plot = "steelblue",
#'   savefile = "histogram_output.png"
#' )
#' }
do_histogram_with_statistics <- function(vect, title, x_axis, text_position, color_plot, savefile = NULL) {
  
  data <- vect
  media <- mean(data, na.rm = TRUE)
  quantili <- quantile(data, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  massimo <- max(data, na.rm = TRUE)
  minimo <- min(data, na.rm = TRUE)
  
  # Open graphics device if saving the plot
  if (!is.null(savefile) && nzchar(savefile)) {
    ext <- tools::file_ext(savefile)
    if (ext == "png") {
      png(savefile, width = 800, height = 600)
    } else if (ext %in% c("jpeg", "jpg")) {
      jpeg(savefile, width = 800, height = 600)
    } else if (ext == "pdf") {
      pdf(savefile)
    } else {
      stop("Unsupported file format. Please use .png, .jpeg, or .pdf.")
    }
  }
  
  # Create the histogram
  hist_data <- hist(data, main = title, xlab = x_axis, ylab = "Frequency", col = color_plot)
  
  # Add text annotations with statistics
  max_y <- max(hist_data$counts)
  y_positions <- seq(from = max_y, length.out = 6, by = -max_y * 0.1)
  labels <- c(
    paste("Min:", round(minimo, 2)),
    paste("Q1:", round(quantili[1], 2)),
    paste("Mean:", round(media, 2)),
    paste("Median:", round(quantili[2], 2)),
    paste("Q3:", round(quantili[3], 2)),
    paste("Max:", round(massimo, 2))
  )
  
  for (i in seq_along(labels)) {
    text(x = text_position, y = y_positions[i], labels = labels[i], pos = 4, col = color_plot, cex = 1.5)
  }
  
  # Close graphics device if saving
  if (!is.null(savefile) && nzchar(savefile)) {
    dev.off()
    message("Plot saved to ", savefile)
  }
}


#' Generate Pie Charts Summarizing Genomic Annotations
#'
#' This function creates pie charts showing the proportion of regions that are annotated
#' vs not annotated for each sample. It takes two lists of data frames (annotated and not annotated),
#' and plots one pie chart per sample.
#'
#' @param not_annotated A named list of data frames, where each data frame contains all the genomic regions (annotated and not annotated).
#'        Each data frame must contain a column named \code{region_id}, and the list names must correspond to sample names.
#' @param annotated A named list of data frames, where each data frame contains the annotated genomic regions.
#'        Each data frame must also contain a column named \code{region_id}, and match the names in \code{not_annotated}.
#' @param ncol_plot Integer, the number of columns to use when arranging the pie charts.
#'
#' @return No return value. Generates and displays a grid of pie charts using \code{ggplot2.multiplot()}.
#' @export
#'
#' @import ggplot2
#' @importFrom easyGgplot2 ggplot2.multiplot
#'
#' @examples
#' \dontrun{
#' annotation_piechart(not_annotated_list, annotated_list, ncol_plot = 3)
#' }
annotation_piechart <- function(not_annotated, annotated, ncol_plot) {
  
  # Checks
  if (!is.list(not_annotated) || !is.list(annotated)) {
    stop("'not_annotated' and 'annotated' must be lists")
  }
  
  if (is.null(names(not_annotated)) || is.null(names(annotated))) {
    stop("Both 'not_annotated' and 'annotated' must be named lists with sample names as names")
  }
  
  if (!all(names(not_annotated) %in% names(annotated)) ||
      !all(names(annotated) %in% names(not_annotated))) {
    stop("'not_annotated' and 'annotated' must have the same sample names")
  }
  
  for (nm in names(not_annotated)) {
    if (!is.data.frame(not_annotated[[nm]]) || !"region_id" %in% colnames(not_annotated[[nm]])) {
      stop(paste0("Element '", nm, "' in 'not_annotated' must be a data frame with a 'region_id' column"))
    }
    if (!is.data.frame(annotated[[nm]]) || !"region_id" %in% colnames(annotated[[nm]])) {
      stop(paste0("Element '", nm, "' in 'annotated' must be a data frame with a 'region_id' column"))
    }
  }
  
  if (!is.numeric(ncol_plot) || length(ncol_plot) != 1 || ncol_plot < 1) {
    stop("'ncol_plot' must be a single positive integer")
  }
  
  # Each plot corresponds to one sample
  plot_vect_annotation <- lapply(names(annotated), function(sample_name) {
    
    # Count regions that are annotated vs not annotated
    annotated_count <- sum(not_annotated[[sample_name]]$region_id %in% annotated[[sample_name]]$region_id)
    not_annotated_count <- sum(!(not_annotated[[sample_name]]$region_id %in% annotated[[sample_name]]$region_id))
    
    values <- c(annotated_count, not_annotated_count)
    labels <- c("Annotated", "Not Annotated")
    
    # Create a data frame for plotting
    df <- data.frame(
      category = labels,
      count = values
    )
    
    # Generate the pie chart
    plotti <- ggplot(df, aes(x = "", y = count, fill = category)) +
      geom_bar(width = 1, stat = "identity") +
      coord_polar("y") +
      labs(title = paste0("DMR Annotation Summary in ", sample_name)) +
      theme_void() +
      scale_fill_manual(values = c("skyblue", "salmon"))
    
    return(plotti)
  })
  
  # Arrange all plots in a grid using easyGgplot2
  do.call(easyGgplot2::ggplot2.multiplot, c(plot_vect_annotation, cols = ncol_plot))
}


#' Generate Pie Charts Showing Annotation Frequency Per Region
#'
#' This function creates pie charts summarizing how frequently each region is annotated across 
#' samples. Regions are grouped into bins based on the number of times they are annotated.
#'
#' @param annotated A named list of data frames. Each data frame should contain a column called 
#'        \code{region_id}, and the list names must correspond to sample names.
#' @param ncol_plot Integer, the number of columns to use when arranging the pie charts.
#'
#' @return No return value. Displays a grid of pie charts using \code{ggplot2.multiplot()}.
#' @export
#'
#' @import ggplot2
#' @importFrom easyGgplot2 ggplot2.multiplot
#'
#' @examples
#' \dontrun{
#' frequency_annotation_piechart(annotated_list, ncol_plot = 3)
#' }
frequency_annotation_piechart <- function(annotated, ncol_plot) {
  
  # Input validation
  if (!is.list(annotated)) {
    stop("'annotated' must be a list")
  }
  if (is.null(names(annotated))) {
    stop("'annotated' must be a named list with sample names")
  }
  for (nm in names(annotated)) {
    if (!is.data.frame(annotated[[nm]]) || !"region_id" %in% colnames(annotated[[nm]])) {
      stop(paste0("Element '", nm, "' in 'annotated' must be a data frame with a 'region_id' column"))
    }
  }
  if (!is.numeric(ncol_plot) || length(ncol_plot) != 1 || ncol_plot < 1) {
    stop("'ncol_plot' must be a single positive integer")
  }
  
  # Each plot corresponds to one sample
  plot_vect_freq_annotation <- lapply(names(annotated), function(sample_name) {
    
    df <- annotated[[sample_name]]
    
    # Count how many times each region_id appears
    region_counts <- as.data.frame(table(df$region_id))
    
    # Group frequencies into bins: 1, 2, 3, 4, 5+
    region_counts$group <- cut(
      region_counts$Freq,
      breaks = c(0, 1, 2, 3, 4, Inf),
      labels = c("1", "2", "3", "4", "5+"),
      right = FALSE
    )
    
    # Count how many regions fall into each group
    group_counts <- as.data.frame(table(region_counts$group))
    
    # Generate the pie chart
    ggplot(group_counts, aes(x = "", y = Freq, fill = Var1)) +
      geom_bar(width = 1, stat = "identity") +
      coord_polar("y") +
      labs(title = paste0("Annotation frequency per region in ", sample_name)) +
      theme_void()
  })
  
  # Arrange all plots in a grid using easyGgplot2
  do.call(easyGgplot2::ggplot2.multiplot, c(plot_vect_freq_annotation, cols = ncol_plot))
}

#' Plot Pie Charts of Gene Loci Distribution per Sample
#'
#' This function creates pie charts showing the relative distribution of annotated gene loci
#' (from the \code{gene_loci} column) for each sample in the input list.
#'
#' @param annotated A named list of data frames. Each data frame must contain a column named 
#'        \code{gene_loci}, and the list names should correspond to sample identifiers.
#' @param ncol_plot Integer, specifying the number of columns in the final multiplot layout.
#'
#' @return No return value. Displays pie charts using \code{ggplot2.multiplot()}.
#' @export
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom easyGgplot2 ggplot2.multiplot
#'
#' @examples
#' \dontrun{
#' plot_geneloci_annotation(annotated_list, ncol_plot = 2)
#' }
plot_geneloci_annotation <- function(annotated, ncol_plot) {
  # Input validation
  if (!is.list(annotated)) {
    stop("'annotated' must be a list")
  }
  if (is.null(names(annotated))) {
    stop("'annotated' must be a named list with sample names")
  }
  for (nm in names(annotated)) {
    if (!is.data.frame(annotated[[nm]]) || !"gene_loci" %in% colnames(annotated[[nm]])) {
      stop(paste0("Element '", nm, "' in 'annotated' must be a data frame with a 'gene_loci' column"))
    }
  }
  if (!is.numeric(ncol_plot) || length(ncol_plot) != 1 || ncol_plot < 1) {
    stop("'ncol_plot' must be a single positive integer")
  }
  # Create one pie chart per sample
  # Palette coerente
  all_loci <- unique(unlist(lapply(annotated, function(df) unique(df$gene_loci))))
  n_colors <- length(all_loci)
  if (n_colors > 12) {
    palette <- setNames(scales::hue_pal()(n_colors), all_loci)
  } else {
    palette <- setNames(RColorBrewer::brewer.pal(n_colors, "Set3"), all_loci)
  }
  
  # Pie plots
  plot_geneloci_annotation_pie <- lapply(names(annotated), function(sample_name) {
    df <- annotated[[sample_name]]
    
    df_pie <- df %>%
      count(gene_loci) %>%
      mutate(
        prop = n / sum(n),
        percent_label = ifelse(prop < 0.02, "", paste0(round(prop * 100), "%"))
      )
    
    ggplot(df_pie, aes(x = "", y = prop, fill = gene_loci)) +
      geom_col(width = 1, color = "white") +
      geom_text(aes(label = percent_label), position = position_stack(vjust = 0.5), color = "black", size = 4) +
      coord_polar(theta = "y") +
      scale_fill_manual(values = palette) +
      labs(
        title = paste0("Distribution of gene_loci in ", sample_name),
        fill = "gene_loci"
      ) +
      theme_void() +
      theme(
        strip.text = element_text(size = 12, face = "bold"),
        legend.position = "right"
      )
  })
  
  # Multiplot
  do.call(easyGgplot2::ggplot2.multiplot, c(plot_geneloci_annotation_pie, cols = ncol_plot))
  
}

#' Plot Multidimensional Scaling (MDS) of Methylation Data
#'
#' This function creates an MDS plot (using `plotMDS`) from methylation data (`df`) 
#' colored by levels of a grouping variable in the metadata. Optionally saves the plot to a file.
#'
#' @param df A data frame or matrix containing methylation values or similar numeric data.
#' @param metadata A data frame containing sample metadata with columns matching sample IDs.
#' @param grouping_variable A character vector of one or more column names in `metadata` to group/color samples.
#' @param top_plot Integer, number of top features (rows) to use for MDS calculation (default 1000).
#' @param save_file Character string file path to save the plot (optional). If missing or NULL, plot is displayed.
#' @param width Numeric, width of saved plot in inches (default 20).
#' @param height Numeric, height of saved plot in inches (default 10).
#'
#' @return NULL. Displays the MDS plot, and optionally saves it to `save_file`.
#' 
#' @importFrom limma plotMDS
#' @export
#'
#' @examples
#' \dontrun{
#' plot_MSD(df = methylation_matrix, metadata = sample_metadata, grouping_variable = "Treatment", top_plot = 500)
#' plot_MSD(df = methylation_matrix, metadata = sample_metadata, grouping_variable = c("Treatment", "Batch"), save_file = "MDSplot.png")
#' }

plot_MSD <- function(df, metadata, grouping_variable,top_plot = 1000, save_file, width = 20, height = 10){
 
  pal <- rainbow(length(unique(metadata[[grouping_variable]]))) # Example palette
  
  # If saving file, open graphics device
  if (!missing(save_file) && !is.null(save_file)) {
    png(filename = save_file, width = width, height = height, units = "in", res = 300)
  }
  
  # Set 1x1 plotting area
  par(mfrow = c(1, 1))
  
  # Plot MDS
  plotMDS(df, top = top_plot, gene.selection = "common",
          col = pal[factor(metadata[[grouping_variable]])])
  
  # Add legend
  legend("topright", legend = levels(factor(metadata[[grouping_variable]])),
         text.col = pal, bg = "white", cex = 1)
  
  # Close graphics device if saving
  if (!missing(save_file) && !is.null(save_file)) {
    dev.off()
  }
  
  invisible(NULL)
  
}

#' Plot Number of Differentially Methylated CpGs or Regions
#'
#' This function generates a bar plot displaying the number of differentially methylated CpGs or regions 
#' across different samples based on Limma analysis results.
#'
#' @param DM_results A list containing differential methylation results from Limma.
#' @param filtering A logical value indicating whether to use filtered results (`TRUE`) or unfiltered results (`FALSE`).
#' @param threshold_pval A numeric value representing the p-value threshold for filtering (if `filtering = TRUE`).
#'
#' @return A ggplot2 bar plot showing the number of differentially methylated CpGs or regions per sample.
#'
#' @details
#' - If `filtering = TRUE`, the plot is based on `DM_results$DM_limma_filtered` with a subtitle indicating the p-value threshold.
#' - If `filtering = FALSE`, the plot is based on `DM_results$DM_limma` with no subtitle.
#' - The x-axis represents samples, and the y-axis represents the number of differentially methylated sites.
#'
#' @examples
#' \dontrun{
#'   plot_differential_methylation(DM_results, filtering = TRUE, threshold_pval = 0.05)
#' }
#'
#' @import ggplot2
#' @export

plot_differential_methylation <- function(DM_results, filtering, threshold_pval) {
  
  # Determine whether to use filtered or unfiltered results
  if (filtering) {
    plot_data <- lapply(DM_results, function(df){
      df[which(df$qvalue < threshold_pval),]
    })
    subtitle <- paste0("p-value < ", threshold_pval)
  } else {
    plot_data <- DM_results
    subtitle <- ""
  }
  
  # Compute the number of differentially methylated CpGs or regions per sample
  lengths_list <- sapply(plot_data, nrow)
  
  # Convert the result into a data frame for plotting
  lengths_df <- data.frame(lengths = lengths_list, sample = names(lengths_list))
  
  # Load ggplot2 for visualization
  library(ggplot2)
  
  # Create the bar plot
  p <- ggplot(lengths_df, aes(x = sample, y = lengths)) +
    geom_col(color = "black", fill = "blue", alpha = 0.7) +  # Create bar chart with blue fill
    labs(
      title = "Number of Differentially Methylated CpGs or Regions",
      subtitle = subtitle,
      x = "Sample",
      y = "Number of CpGs or Regions"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability
  
  # Print the plot
  print(p)
}