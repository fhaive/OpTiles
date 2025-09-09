#' Load and Process Methylation Data from Bismark Coverage Files
#'
#' Reads Bismark coverage files (`*.bismark.cov.gz`) from a directory and processes methylation data
#' using metadata with sample names and treatment information.
#'
#' @param metadata A data frame containing metadata for the samples. It must include:
#'   - A column (name given by `sample_name_variable`) whose values are reported also in
#'       the Bismark coverage file names (excluding the `.bismark.cov.gz` extension).
#'   - A numeric vector column (name given by `treatment_variable`) encoding the treatment groups,
#'       with length equal to the number of samples/files.
#' @param repository Path to directory containing Bismark coverage files where the samples are ordered as in the metadata(`*.bismark.cov.gz`).
#' @param sample_name_variable Character string specifying the metadata column with sample names.
#' @param treatment_variable Character string specifying the metadata column with numeric treatment encoding.
#' @param mincov Minimum coverage threshold per CpG site (integer >= 1).
#' @param individual_sample_plot Logical, whether to generate coverage and methylation plots per sample.
#' @param filter_high_percentage Logical, whether to filter out sites with extremely high coverage (top 0.1%).
#' @param normalization_coverage Logical, whether to normalize coverage across samples before uniting.
#' @param method_normalization Character string specifying the normalization method to use if \code{normalization_coverage = TRUE}. Must be one of: `"median"` (default) or `"mean"`.
#' @param output_folder Path to output folder for saving intermediate files (currently unused).
#' @param assembly Genome assembly name (default: `"hg38"`).
#'
#' @return A \code{methylBase} object representing united methylation data across all samples.
#'
#' @details
#' The metadata must meet the following requirements:
#'   The sample names in `metadata[[sample_name_variable]]` must exactly match
#'         the coverage file names (without extensions) in `repository`.
#'   The treatment column specified by `treatment_variable` must be numeric and have the same length as the sample names.
#'
#' The function uses \code{methylKit} to read, normalize, filter, and unite methylation data.
#' If \code{individual_sample_plot} is TRUE, plots for methylation and coverage statistics per sample are generated.
#' If \code{normalization_coverage} is TRUE, coverage will be normalized using the specified method (\code{"median"} or \code{"mean"}) before uniting data.
#' @importFrom methylKit methRead unite normalizeCoverage filterByCoverage
#' @export
load_methylation_data <- function(metadata, repository, sample_name_variable, 
                                  treatment_variable,
                                  mincov, 
                                  individual_sample_plot = TRUE, 
                                  filter_high_percentage = TRUE,
                                  normalization_coverage = TRUE,
                                  method_normalization = "median",
                                  output_folder, assembly = "hg38") {
  
  # Validate inputs
  if (!is.data.frame(metadata)) stop("metadata must be a data frame")
  if (!(sample_name_variable %in% colnames(metadata))) stop("sample_name_variable column not found in metadata")
  if (!(treatment_variable %in% colnames(metadata))) stop("treatment_variable column not found in metadata")
  if (length(metadata[[sample_name_variable]]) != length(metadata[[treatment_variable]])) {
    stop("Length of sample_name_variable and treatment_variable columns must be equal")
  }
  if (!is.numeric(metadata[[treatment_variable]])) stop("treatment_variable column must be numeric")
  if (!dir.exists(repository)) stop("repository folder does not exist")
  if (!is.numeric(mincov) || mincov < 1) stop("mincov must be numeric and >= 1")
  
  # List Bismark coverage files and extract base names
  file_list <- list.files(path = repository, pattern = "bismark.cov.gz$", full.names = TRUE)
  if (length(file_list) != nrow(metadata)) {
    stop("Number of Bismark files in 'repository' does not match the number of samples in metadata.")
  }
  if (length(file_list) != length(metadata[[sample_name_variable]])) {
    stop("Mismatch between number of Bismark files and sample names in metadata.")
  }
  file_names <- basename(file_list)
  file_names <- sub("\\.bismark\\.cov\\.gz$", "", file_names)
  
  # Prepare lists for methylKit input
  lst_files <- lapply(file_list, function(x) x)
  ord2 <- sapply(metadata[[sample_name_variable]], function(x) grep(paste0(x,"_"), lst_files))
  metadata <- metadata[order(ord2),]
  lst_file_names <- lapply(metadata[[sample_name_variable]], identity)
  
  treatment_vector <- metadata[[treatment_variable]]
  
  # Read methylation data
  myobj <- methRead(location = lst_files,
                    sample.id = lst_file_names,
                    assembly = assembly,
                    treatment = treatment_vector,
                    context = "CpG",
                    mincov = mincov,
                    pipeline = "bismarkCoverage",
                    header = FALSE)
  
  # Generate QC plots per sample if requested
  if (individual_sample_plot) {
    for (i in seq_along(myobj)) {
      cat("Generating plot for sample:", myobj[[i]]@sample.id, "\n")
      getMethylationStats(myobj[[i]], plot = FALSE, both.strands = FALSE)
      getCoverageStats(myobj[[i]], plot = TRUE, both.strands = FALSE)
    }
    cat("Plots can be viewed in the 'Plots' section.\n")
  }
  
  # Filter out extremely high coverage CpGs (top 0.1%) if enabled
  if (filter_high_percentage) {
    myobj <- filterByCoverage(myobj, hi.perc = 99.9)
  }
  
  # Normalize coverage across samples using median method
  if (normalization_coverage) {
    myobj <- normalizeCoverage(myobj, method = method_normalization)
  }
  
  # Unite data, keeping CpGs present in all samples
  meth <- methylKit::unite(myobj, destrand = FALSE)
  
  return(meth)
}
