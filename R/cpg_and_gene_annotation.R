
#' Download Annotations from BioMart
#'
#' This function retrieves gene annotations from the BioMart database based on user-specified parameters.
#'
#' @param biomart_database A string specifying the BioMart database to use (e.g., "ensembl").
#' @param biomart_dataset A string specifying the dataset to query (e.g., "hsapiens_gene_ensembl").
#' @param chromosome_vector A character vector containing the chromosomes to filter the results.
#' @param selected_attributes A character vector of attribute names to be retrieved from BioMart.
#' @param attributes_page A character vector specifying the attribute pages to filter relevant attributes.
#' @param gene_biotype A string specifying the filter name for gene biotypes (e.g., "biotype").
#' @param filter_gene_biotype A character vector containing the gene biotypes to filter for retrieval.
#'
#' @return A named list where each element corresponds to a different attribute page and contains a data frame with the requested annotations.
#'
#' @import biomaRt dplyr tidyr gridExtra grid ggplot2
#' @export
#'
#' @examples
#' \dontrun{
#' biomart_info <- download_biomart_annotation(
#'   biomart_database = "ensembl",
#'   biomart_dataset = "hsapiens_gene_ensembl",
#'   chromosome_vector = c("1", "2", "3"),
#'   selected_attributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name"),
#'   attributes_page = c("feature_page"),
#'   gene_biotype = "biotype",
#'   filter_gene_biotype = c("protein_coding")
#' )
#' }
#' 

download_biomart_annotation <- function(biomart_database, biomart_dataset,chromosome_vector, selected_attributes, attributes_page, gene_biotype, filter_gene_biotype){
  
  library(biomaRt)
  library(dplyr)
  library(tidyr)
  library(gridExtra)
  library(grid)
  library(ggplot2)
  
  #Downloading the data from biomaRt
  ensembl = useEnsembl(biomart = biomart_database, dataset = biomart_dataset)
  normal_chrs = chromosome_vector
  
  attributes <- biomaRt::listAttributes(ensembl)
  attributes <- attributes[attributes$page %in% attributes_page,]
  
  pages_retrive_attributed <- split(attributes, attributes$page )
  biomart_info <- list()
  for(page in names(pages_retrive_attributed)){
    
    aa <- selected_attributes[selected_attributes %in% pages_retrive_attributed[[page]][["name"]]]
    values <- list()
    biomart_info[[page]] <- biomaRt::getBM(attributes = aa,
                                  filters = gene_biotype,
                                  values = filter_gene_biotype,
                                  mart = ensembl)
  }
  
  return(biomart_info)
  
}

#' Extract Gene Regions from a DataFrame
#'
#' This function extracts different gene regions (entire gene, UTRs, exons, and promoter) 
#' from a given dataframe containing genomic annotations.
#'
#' @param df A dataframe containing genomic annotations. Required columns:
#'   \itemize{
#'     \item \code{hgnc_symbol} - Gene symbol
#'     \item \code{chromosome_name} - Chromosome identifier
#'     \item \code{start_position}, \code{end_position} - Gene start and end positions
#'     \item \code{strand} - Strand information (+1 or -1)
#'     \item \code{5_utr_start}, \code{5_utr_end} - 5' UTR start and end positions
#'     \item \code{3_utr_start}, \code{3_utr_end} - 3' UTR start and end positions
#'     \item \code{exon_chrom_start}, \code{exon_chrom_end} - Exon start and end positions
#'     \item \code{rank} - Exon rank
#'     \item \code{transcription_start_site} - Transcription start site
#'   }
#' @param promoter_distance A numeric vector of length 2 specifying the upstream 
#'   and downstream distances to define the promoter region.
#'
#' @return A dataframe with the extracted gene regions, including:
#'   \itemize{
#'     \item \code{gene} - Gene symbol
#'     \item \code{location} - Genomic location (formatted as chr:start-end)
#'     \item \code{gene_loci} - Type of region (e.g., entire_gene, 5UTR, 3UTR, exon#, promoter)
#'     \item \code{strand} - Strand information
#'   }
#'
#' @examples
#' df <- data.frame(
#'   hgnc_symbol = c("GENE1", "GENE2"),
#'   chromosome_name = c("1", "X"),
#'   start_position = c(1000, 2000),
#'   end_position = c(5000, 7000),
#'   strand = c(1, -1),
#'   transcription_start_site = c(1200, 2100)
#' )
#' promoter_distance <- c(-2000, 500)
#' extract_gene_regions(df, promoter_distance)
#'
#' @export
extract_gene_regions <- function(df, promoter_distance) {
  
  required_columns <- c("chromosome_name","start_position","end_position","strand",
                        "5_utr_start","5_utr_end",
                        "3_utr_start","3_utr_end",
                        "exon_chrom_start","exon_chrom_end","rank",
                        "transcription_start_site")
  
  missing_columns <- dplyr::setdiff(required_columns, colnames(df))
  extra_columns <- dplyr::setdiff(colnames(df), required_columns)
  
  if (length(missing_columns) > 0) {
    stop("Missing columns: ", paste(missing_columns, collapse = ", "))
  }
  if (length(extra_columns) > 0) {
    message("Extra columns found in df: ", paste(extra_columns, collapse = ", "))
  }
  
  gene_name_col <- grep("symbol",colnames(df))
  gene <- unique(df[,gene_name_col])
  gene_name_col <- colnames(df)[gene_name_col]
  
  df[df == ""] <- NA
  df <- df[!is.na(df[[gene_name_col]]),]
  
  # Entire gene
  eg <- unique(df[, c(gene_name_col, "chromosome_name", "start_position", "end_position", "strand","transcription_start_site","3_utr_start","3_utr_end")])

  entire_gene <- data.frame(
    gene = eg[[gene_name_col]],
    location = paste0("chr", eg$chromosome_name, ":", eg$start_position, "-", eg$end_position),
    gene_loci = "entire_gene",
    strand = eg$strand
  )
  
  entire_gene <- unique(entire_gene)
  
  tss <- unique(df[, c(gene_name_col, "chromosome_name","transcription_start_site","strand")])
  tss_df <- data.frame(
    gene = tss[[gene_name_col]],
    location = paste0("chr", tss$chromosome_name, ":", tss[["transcription_start_site"]], "-", tss[["transcription_start_site"]]),
    gene_loci = "TSS",
    strand = tss$strand
  )
  
  
  # 5' UTR
  utr5 <- na.omit(unique(df[, c(gene_name_col, "chromosome_name", "5_utr_start", "5_utr_end", "strand")]))
  utr5_df <- if (nrow(utr5) > 0) {
    data.frame(
      gene = utr5[[gene_name_col]],
      location = paste0("chr", utr5$chromosome_name, ":", utr5[["5_utr_start"]], "-", utr5[["5_utr_end"]]),
      gene_loci = "5UTR",
      strand = utr5$strand
    )
  } else {
    data.frame(gene = character(), location = character(), gene_loci = character())
  }
  
  # 3' UTR
  utr3 <- na.omit(unique(df[, c(gene_name_col, "chromosome_name", "3_utr_start", "3_utr_end", "strand")]))
  utr3_df <- if (nrow(utr3) > 0) {
    data.frame(
      gene = utr3[[gene_name_col]],
      location = paste0("chr", utr3$chromosome_name, ":", utr3[["3_utr_start"]], "-", utr3[["3_utr_end"]]),
      gene_loci = "3UTR",
      strand = utr3$strand
    )
  } else {
    data.frame(gene = character(), location = character(), gene_loci = character())
  }
  
  # Exons
  exons <- na.omit(unique(df[, c(gene_name_col, "chromosome_name", "exon_chrom_start", "exon_chrom_end", "rank", "strand")]))
  exon_df <- if (nrow(exons) > 0) {
    data.frame(
      gene = exons[[gene_name_col]],
      location = paste0("chr", exons$chromosome_name, ":", exons$exon_chrom_start, "-", exons$`exon_chrom_end`),
      gene_loci = paste0("exon", exons$rank),
      strand = exons$strand
    )
  } else {
    data.frame(gene = character(), location = character(), gene_loci = character())
  }

  if(!any(is.na(promoter_distance)) && length(promoter_distance) == 2){
    # Promoter
    promoter <- unique(df[, c(gene_name_col, "chromosome_name", "transcription_start_site", "strand")])
    promoter_df <- data.frame(
      gene = promoter[[gene_name_col]],
      # location =  paste0("chr", promoter$chromosome_name, ":", 
      #                    promoter$transcription_start_site - promoter_distance[1], "-",
      #                    promoter$transcription_start_site + promoter_distance[2]),
      location = ifelse(promoter$strand == 1,
                        paste0("chr", promoter$chromosome_name, ":",
                               promoter$transcription_start_site - promoter_distance[1], "-",
                               promoter$transcription_start_site + promoter_distance[2]),
                        paste0("chr", promoter$chromosome_name, ":",
                               promoter$transcription_start_site - promoter_distance[2], "-",
                               promoter$transcription_start_site + promoter_distance[1])),
      gene_loci = "promoter",
      strand = promoter$strand
    )
    final_df <- rbind(entire_gene, utr5_df, utr3_df, exon_df, promoter_df,tss_df)
  } else {
    final_df <- rbind(entire_gene, utr5_df, utr3_df, exon_df, tss_df)
  }
  
  # Combine all dataframes
  final_df <- final_df[order(final_df$gene), ]
  
  return(final_df)
}


#' Retrieve and Annotate Gene Regions Using BioMart
#'
#' This function retrieves gene annotation data from the Ensembl BioMart database, merges the results, 
#' and extracts specific gene regions such as promoters, UTRs, exons, and entire genes.
#'
#' @param biomart_database A string specifying the BioMart database to use (e.g., \code{"ensembl"}).
#' @param biomart_dataset A string specifying the dataset within the BioMart database (e.g., \code{"hsapiens_gene_ensembl"}).
#' @param chromosome_vector A character vector specifying the chromosomes to filter (e.g., \code{c("1", "2", "X")}).
#' @param selected_attributes A character vector of attributes to retrieve from BioMart.
#' @param attributes_page A character string specifying the attribute page from which to select attributes.
#' @param gene_biotype A string specifying the gene biotype filter (e.g., \code{"protein_coding"}).
#' @param filter_gene_biotype A character vector of gene biotypes to filter.
#' @param promoter_distance A numeric vector of length 2, respectively the upstream and downstream distances 
#'   to define the promoter region.
#'
#' @return A dataframe with annotated gene regions, including:
#'   \itemize{
#'     \item \code{gene} - Gene symbol
#'     \item \code{location} - Genomic location (formatted as chr:start-end)
#'     \item \code{gene_loci} - Type of region (e.g., entire_gene, 5UTR, 3UTR, exon#, promoter)
#'     \item \code{strand} - Strand information
#'   }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Downloads gene annotation data from the specified BioMart database using \code{download_biomart_annotation()}.
#'   \item Identifies common attributes across retrieved datasets and merges them into a single dataframe.
#'   \item Extracts specific gene regions using \code{extract_gene_regions()}.
#' }
#'
#' @examples
#' \dontrun{
#' biomart_annotation(
#'   biomart_database = "ensembl",
#'   biomart_dataset = "hsapiens_gene_ensembl",
#'   chromosome_vector = c("1", "2", "X"),
#'   selected_attributes = c("hgnc_symbol", "chromosome_name", "start_position"),
#'   attributes_page = "feature_page",
#'   gene_biotype = "protein_coding",
#'   filter_gene_biotype = c("protein_coding"),
#'   promoter_distance = c(2000, 500)
#' )
#' }
#'
#' @export
biomart_annotation <- function(biomart_database, biomart_dataset, chromosome_vector, 
                               selected_attributes, attributes_page, gene_biotype, 
                               filter_gene_biotype, promoter_distance = NA) {
  
  biomart_df <- download_biomart_annotation(biomart_database, biomart_dataset, 
                                            chromosome_vector, selected_attributes, 
                                            attributes_page, gene_biotype, filter_gene_biotype)
  
  common_attributes <- Reduce(dplyr::intersect, lapply(biomart_df, colnames))
  
  merge_df <- Reduce(function(x, y) merge(x, y, by = common_attributes, all = TRUE), biomart_df)
  
  annotated_gene_regions <- extract_gene_regions(merge_df, promoter_distance)
  
  return(annotated_gene_regions)
}

#' Map CpG Sites to Genomic Regions
#'
#' This function maps CpG sites to user-defined genomic regions and returns all overlapping CpG-to-region mappings.
#'
#' @param methdf A `data.frame` or `MethylBase` object with methylation data. Must include the columns \code{chr}, \code{start}, and \code{end}
#'   indicating the chromosome and coordinates of each CpG site.
#' @param region A `data.frame` or `MethylBase` object with genomic regions. Must include the columns \code{chr}, \code{start}, and \code{end}.
#'
#' @return A `data.frame` containing the CpGs that overlap each region, including:
#' \describe{
#'   \item{`region_id`}{Region identifier in the format \code{"chr:start-end"}.}
#'   \item{`chr`, `start`, `end`}{Genomic coordinates of the CpG site.}
#' }
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Converts CpG and region data into `GRanges` objects.
#'   \item Identifies overlaps using \code{\link[GenomicRanges]{findOverlaps}}.
#'   \item Returns a unique list of CpG-to-region mappings.
#' }
#'
#' @import GenomicRanges
#' @import data.table
#' @importFrom IRanges IRanges
#'
#' @examples
#' \dontrun{
#' methdf <- data.frame(chr = "1", start = c(100, 150, 250), end = c(100, 150, 250))
#' region <- data.frame(chr = "1", start = c(90, 200), end = c(160, 300))
#' map_cpg_to_regions(methdf, region)
#' }
#'
#' @export
map_cpg_to_regions <- function(methdf, region) {
  
  if (!is.numeric(methdf$start) || !is.numeric(methdf$end)) {
    stop("The 'start' and 'end' columns in 'methdf' must be numeric.")
  }
  
  if (!is.numeric(region$start) || !is.numeric(region$end)) {
    stop("The 'start' and 'end' columns in 'region' must be numeric.")
  }
  # Convert CpG sites to GRanges object
  sites2 <- GRanges(
    seqnames = methdf$chr,
    ranges = IRanges(start = methdf$start, end = methdf$end)
  )
  
  # Convert genomic regions to GRanges with a unique region ID
  regions2 <- GRanges(
    seqnames = region$chr,
    ranges = IRanges(start = region$start, end = region$end),
    region_id = paste0("chr", region$chr, ":", region$start, "-", region$end)
  )
  
  # Identify overlaps between CpG sites and regions
  overlaps <- GenomicRanges::findOverlaps(regions2, sites2)
  
  # Extract mapping of CpG sites to regions
  results <- data.frame(
    region_id = regions2$region_id[queryHits(overlaps)],
    chr       = as.character(seqnames(sites2[subjectHits(overlaps)])),
    start     = start(sites2[subjectHits(overlaps)]),
    end       = start(sites2[subjectHits(overlaps)])
  )
  
  # Convert to data.table for performance and remove duplicates
  setDT(results)
  regions_cpgs <- as.data.frame(unique(results))
  
  # Return final region-CpG mapping
  return(regions_cpgs)
}



#' Map Genomic Regions to Annotated Regions and Calculate Overlap Statistics
#'
#' This function maps a set of genomic regions (e.g., CpG tiles or peaks) to annotated genomic regions 
#' (e.g., genes, enhancers) based on coordinate overlaps, and computes overlap lengths and percentages.
#'
#' @param regions_to_map A `data.frame` containing genomic regions to be mapped. Must include a `region_id` 
#'   column in the format `"chrN:start-end"` (e.g., `"chr1:100-200"`).
#' @param annotation A `data.frame` containing annotated regions. Must include a `location` column in 
#'   the same `"chrN:start-end"` format.
#'
#' @return A `list` with two elements:
#' \describe{
#'   \item{`annotated_df`}{A `data.table` with overlaps between input regions and annotations, including: 
#'   \code{annotation_regions}, \code{merged_region}, \code{overlap_bp}, \code{pct_tile}, 
#'   \code{pct_merged}, and \code{overlap_window}.}
#'   \item{`overlap_regions`}{A `data.frame` containing the genomic coordinates of the overlapping regions.}
#' }
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Parses region coordinates from `chr:start-end` format.
#'   \item Converts input data to `GRanges` objects.
#'   \item Identifies overlapping regions using \code{\link[GenomicRanges]{findOverlaps}}.
#'   \item Computes the width of overlap and its percentage with respect to both input and annotation regions.
#'   \item Returns a table of overlaps and summarized information.
#' }
#'
#' @import data.table
#' @import GenomicRanges
#' @importFrom IRanges IRanges
#'
#' @examples
#' \dontrun{
#' regions_to_map <- data.frame(region_id = c("chr1:100-200", "chr1:150-250"))
#' annotation <- data.frame(location = c("chr1:180-300", "chr1:50-120"))
#' result <- map_regions_to_regions(regions_to_map, annotation)
#' print(result$annotated_df)
#' }
#'
#' @export
map_regions_to_regions <- function(regions_to_map, annotation) {
  
  # Load required package
  library(data.table)
  # Input validation
  if (!is.data.frame(regions_to_map)) {
    stop("'regions_to_map' must be a data.frame")
  }
  if (!"region_id" %in% colnames(regions_to_map)) {
    stop("'regions_to_map' must contain a 'region_id' column")
  }
  if (!is.character(regions_to_map$region_id)) {
    stop("'region_id' column in 'regions_to_map' must be of type character")
  }
  
  if (!is.data.frame(annotation)) {
    stop("'annotation' must be a data.frame")
  }
  if (!"location" %in% colnames(annotation)) {
    stop("'annotation' must contain a 'location' column")
  }
  if (!is.character(annotation$location)) {
    stop("'location' column in 'annotation' must be of type character")
  }
  
  # Validate format of 'region_id' and 'location' columns (basic check)
  validate_region_format <- function(x) {
    all(grepl("^chr?[0-9XYM]+:[0-9]+-[0-9]+$", x))
  }
  
  if (!validate_region_format(regions_to_map$region_id)) {
    message("All 'region_id' entries in 'regions_to_map' must be in the format 'chrN:start-end'")
  }
  if (!validate_region_format(annotation$location)) {
    message("All 'location' entries in 'annotation' that are mapped are in the format 'chrN:start-end'")
  }
  # -- Step 1: Parse 'region_id' strings into chrom, start, end (remove 'chr' prefix first)
  split_sample <- tstrsplit(gsub("^chr", "", regions_to_map$region_id), "[:-]", type.convert = TRUE)
  merged_regions <- data.frame(
    chr = split_sample[[1]],
    start = split_sample[[2]],
    end = split_sample[[3]],
    region_id = regions_to_map$region_id
  )
  
  # -- Step 2: Parse 'location' strings from annotation
  split_annot <- tstrsplit(gsub("^chr", "", annotation$location), "[:-]", type.convert = TRUE)
  annotation_regions <- data.table(
    chr = split_annot[[1]],
    start = split_annot[[2]],
    end = split_annot[[3]],
    region_id = annotation$location
  )
  
  annotation_regions <- na.omit(annotation_regions)
  
  # -- Step 3: Convert parsed regions into GRanges objects
  gr_merged <- GRanges(
    seqnames = merged_regions$chr,
    ranges = IRanges(start = merged_regions$start, end = merged_regions$end),
    region_id = merged_regions$region_id
  )
  
  gr_annotation <- GRanges(
    seqnames = annotation_regions$chr,
    ranges = IRanges(start = annotation_regions$start, end = annotation_regions$end),
    region_id = annotation_regions$region_id
  )
  
  # -- Step 4: Find overlaps between annotation and merged regions
  hits <- findOverlaps(gr_annotation, gr_merged)
  
  # Handle no overlaps case:
  if (length(hits) == 0) {
    # Create empty data.frame for annotated_df with the correct columns
    annotated_df <- data.table(
      annotation_regions = character(),
      merged_region      = character(),
      overlap_bp         = numeric(),
      pct_tile           = numeric(),
      pct_merged         = numeric(),
      overlap_window     = character()
    )
    
    overlap_regions <- data.frame(
      chr = character(),
      start = integer(),
      end = integer(),
      overlap_location = character()
    )
    
    return(list(
      annotated_df = annotated_df,
      overlap_regions = overlap_regions
    ))
  }
  
  # -- Step 5: Extract overlapping regions from GRanges hits
  annotation_hit <- gr_annotation[queryHits(hits)]
  merged_hit     <- gr_merged[subjectHits(hits)]
  
  # -- Step 6: Calculate intersecting (overlapping) regions
  overlap <- pintersect(annotation_hit, merged_hit)
  
  # -- Step 7: Format overlapping regions into chr:start-end format
  overlap_regions <- data.frame(
    chr   = as.character(seqnames(overlap)),
    start = start(overlap),
    end   = end(overlap)
  )
  
  overlap_regions$overlap_location <- paste0(
    "chr", overlap_regions$chr, ":", overlap_regions$start, "-", overlap_regions$end
  )
  
  # -- Step 8: Compute widths and overlap percentages
  overlap_width <- width(overlap)
  width_tiles   <- width(annotation_hit)
  width_merged  <- width(merged_hit)
  
  pct_tiles  <- overlap_width / width_tiles * 100
  pct_merged <- overlap_width / width_merged * 100
  
  # -- Step 9: Create final annotated data frame
  annotated_df <- data.frame(
    annotation_regions = annotation_hit$region_id,
    merged_region      = merged_hit$region_id,
    overlap_bp         = overlap_width,
    pct_tile           = pct_tiles,
    pct_merged         = pct_merged,
    overlap_window     = overlap_regions$overlap_location
  )
  
  # -- Step 10: Ensure unique rows
  setDT(annotated_df)
  annotated_df <- unique(annotated_df)
  
  # -- Step 11: Return results as a named list
  result <- list(
    annotated_df    = annotated_df,
    overlap_regions = overlap_regions
  )
  
  return(result)
}

#' Map Genomic Regions and CpG Sites to Annotated Genomic Windows
#'
#' This function maps sample regions and CpG sites to annotated genomic regions,
#' computes overlaps, and summarizes how many CpG sites overlap each annotation window.
#'
#' @param sample_to_map A `data.frame` containing the sample regions and CpG coordinates to be mapped.
#'   Must contain columns: \code{region_id}, \code{chr}, \code{start}, \code{end}.
#' @param annotation A `data.frame` of annotated genomic regions. Must contain at least the column:
#'   \code{location} (in "chr:start-end" format). If \code{origin_annotation = "gene"}, it must also contain
#'   columns \code{gene} and \code{gene_loci}.
#' @param origin_annotation A `character` indicating the type of annotation:
#'   \itemize{
#'     \item \code{"gene"} — gene-based annotations. Requires \code{gene} and \code{gene_loci} columns in \code{annotation}.
#'     \item \code{"enhancer"} — adds "enhancer" as value to \code{gene_loci}.
#'     \item \code{"CpGI"} or \code{"CpGIslands"} — adds "CpGI" as value to \code{gene_loci}.
#'     \item \code{NULL} or {NA} — gene and gene_loci columns will be NA.
#'   }
#'
#' @return A `data.frame` with:
#' \describe{
#'   \item{`annotation_regions`}{Original annotation region ID.}
#'   \item{`region_id`}{ID of the region from \code{sample_to_map} overlapping the annotation.}
#'   \item{`overlap_window`}{Genomic coordinates of the overlapping window.}
#'   \item{`overlap_bp`}{Number of overlapping base pairs.}
#'   \item{`pct_tile`}{Percent overlap relative to the annotation region.}
#'   \item{`pct_merged`}{Percent overlap relative to the sample region.}
#'   \item{`overlapping_NCpGs`}{Number of CpG sites overlapping the region.}
#'   \item{`gene`, `gene_loci`}{(Optional) Gene and gene loci associated with the annotation region, if applicable.}
#' }
#'
#' @details
#' This function performs a two-step mapping:
#' \enumerate{
#'   \item Maps sample regions to annotations via overlap (using \code{\link{map_regions_to_regions}}).
#'   \item Maps CpG sites to the overlapping windows (using \code{\link{map_cpg_to_regions}}).
#' }
#' It then merges the results and optionally adds gene or enhancer information.
#'
#' @seealso \code{\link{map_regions_to_regions}}, \code{\link{map_cpg_to_regions}}
#'
#' @import data.table
#' @export
map_regions <- function(sample_to_map, annotation, origin_annotation = NA) {
  
  # Input checks
  required_sample_cols <- c("chr", "start", "end")
  missing_sample_cols <- dplyr::setdiff(required_sample_cols, colnames(sample_to_map))
  if (length(missing_sample_cols) > 0) {
    stop("sample_to_map is missing required columns: ", paste(missing_sample_cols, collapse = ", "))
  }
  
  if(!any(colnames(sample_to_map) %in% "region_id")){
    sample_to_map$region_id <- paste0("chr",sample_to_map$chr,":",sample_to_map$start,"-",sample_to_map$end)
  }
  
  if (!("location" %in% colnames(annotation))) {
    stop("annotation must contain a 'location' column")
  }
  
  valid_origins <- c(NULL, "gene", "enhancer", "CpGI", "CpGIslands", NA)
  if (!(origin_annotation %in% valid_origins)) {
    stop("Invalid 'origin_annotation' argument. Allowed values: NULL, NA, 'gene', 'enhancer', 'CpGI', 'CpGIslands'")
  }
  
  if (origin_annotation == "gene") {
    required_gene_cols <- c("gene", "gene_loci")
    missing_gene_cols <- dplyr::setdiff(required_gene_cols, colnames(annotation))
    if (length(missing_gene_cols) > 0) {
      stop("For gene-based annotation, 'annotation' must contain columns: ", paste(missing_gene_cols, collapse = ", "))
    }
  }
  # Step 1: Map sample regions to annotation regions based on overlaps
  annotation_reg_to_reg <- map_regions_to_regions(
    regions_to_map = sample_to_map,
    annotation = annotation
  )
  
  # Step 2: Map CpG sites from sample to overlap windows from step 1
  annotation_cpgs_to_overlap_regions <- map_cpg_to_regions(
    methdf = sample_to_map,
    region = annotation_reg_to_reg$overlap_regions
  )
  
  # Create CpG identifiers
  annotation_cpgs_to_overlap_regions <- data.frame(
    overlap_window = annotation_cpgs_to_overlap_regions$region_id,
    cpg_id = paste0("chr",
                    annotation_cpgs_to_overlap_regions$chr, ":",
                    annotation_cpgs_to_overlap_regions$start, "-",
                    annotation_cpgs_to_overlap_regions$end)
  )
  
  # Merge region-region annotation with CpG mappings
  merged_annotation <- merge(
    as.data.frame(annotation_reg_to_reg$annotated_df),
    annotation_cpgs_to_overlap_regions,
    by = "overlap_window",
    all.x = TRUE
  )
  
  # Count number of CpG overlaps per annotation + window
  setDT(merged_annotation)
  ov_freq <- merged_annotation[, .(overlapping_NCpGs = .N),
                               by = .(annotation_regions, overlap_window)]
  
  # Merge counts into the full mapping table
  merged_df2 <- merge(
    merged_annotation,
    as.data.frame(ov_freq),
    by = c("annotation_regions", "overlap_window"),
    all.x = TRUE
  )
  
  # Handle regions with no CpG overlap
  merged_df2$overlapping_NCpGs[is.na(merged_df2$cpg_id)] <- 0
  
  # Select relevant columns and ensure uniqueness
  final_df <- merged_df2[, c("annotation_regions", "merged_region",
                             "overlap_window", "overlap_bp", "pct_tile",
                             "pct_merged", "overlapping_NCpGs")]
  final_df <- unique(final_df)
  
  # Optional annotation columns based on origin_annotation
  if (!is.null(origin_annotation)) {
    if (origin_annotation == "gene") {
      if (all(c("gene", "gene_loci") %in% colnames(annotation))) {
        final_df[, c("gene", "gene_loci")] <- annotation[
          match(final_df$annotation_regions, annotation$location),
          c("gene", "gene_loci")
        ]
      } else {
        stop("Missing required columns in 'annotation'. To perform gene-based annotation, the data frame must contain both 'gene' and 'gene_loci' columns.")
      }
    } else if (origin_annotation == "enhancer") {
      final_df$gene_loci <- "enhancer"
      final_df$gene <- NA
    } else if (origin_annotation %in% c("CpGI", "CpGIslands")) {
      final_df$gene_loci <- "CpGI"
      final_df$gene <- NA
    }
  } else {
    final_df$gene_loci <- NA
    final_df$gene <- NA
  }
  
  # Rename column for consistency
  colnames(final_df)[2] <- "region_id"
  
  return(final_df)
}

