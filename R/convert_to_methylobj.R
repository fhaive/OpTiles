#' Convert Genomic Region Strings to a Methylation Object
#'
#' This function converts a character vector of genomic region strings into a format compatible with 
#' methylation analysis using the \code{methylKit} package. It constructs genomic ranges from the input,
#' matches them to CpG methylation data, and returns a region-level methylation object.
#'
#' @param regions_location A character vector of genomic regions in the format \code{"chr:start-end"} 
#'                         (e.g., \code{"chr1:10000-10100"}). Each element represents a genomic region.
#' 
#' @param cpgs_methyl_obj A \code{methylRawList} or \code{methylRawListDB} object containing single-base resolution 
#'                        CpG methylation data from the \code{methylKit} package. This object is used to extract methylation 
#'                        counts over the specified regions.
#'
#' @return A \code{methylRegion} object containing summarized methylation information over the specified regions.
#'
#' @details
#' The function parses each region string into chromosome, start, and end coordinates, constructs a 
#' \code{GRanges} object, and uses \code{regionCounts} from the \code{methylKit} package to summarize methylation
#' values for each region. It is useful for summarizing methylation across custom genomic regions.
#'
#' @importFrom dplyr select
#' @importFrom tidyr separate
#' @import GenomicRanges
#' @importFrom IRanges IRanges
#' @importFrom methylKit regionCounts
#'
#' @examples
#' \dontrun{
#' regions <- c("chr1:1000-1100", "chr2:2000-2100")
#' cpg_meth <- unite(my_methylRawList)
#' region_meth <- convert_to_methylobj(regions, cpg_meth)
#' }
#'
#' @export

convert_to_methylobj <- function(regions_location, cpgs_methyl_obj){
  
  individual_regiond_location <- data.frame(region = regions_location) %>%
    separate(region, into = c("chr", "positions"), sep = ":", remove = FALSE) %>%
    separate(positions, into = c("start", "end"), sep = "-") %>%
    dplyr::select(chr,start,end)
  
  individual_regiond_location$start <- as.integer(individual_regiond_location$start)
  individual_regiond_location$end <- as.integer(individual_regiond_location$end)
  individual_regiond_location$chr <- (sub("^chr([^:]+).*", "\\1", individual_regiond_location$chr))
  individual_regiond_location <- unique(individual_regiond_location)
  
  # Convert CpG sites to GRanges object
  regions <- GRanges(
    seqnames = individual_regiond_location$chr,
    ranges = IRanges(start = as.integer(individual_regiond_location$start), 
                     end = as.integer(individual_regiond_location$end)),
    cpg_id = individual_regiond_location$region
  )
  
  region_obj <- regionCounts(cpgs_methyl_obj, regions)
  
  return(region_obj)
}