################################################################################
# Methylation Analysis using OpTiles
# This script demonstrates how to perform methylation data analysis, 
# region optimization, filtering, and annotation using OpTiles and related packages.
################################################################################

#------------------------------------------------------------------------------
# Load necessary libraries ----
#------------------------------------------------------------------------------
# These packages provide tools for methylation analysis, data manipulation, 
# visualization, annotation, and input/output operations.

# Load necessary libraries ----

library(OpTiles)
library(dplyr)
library(methylKit)
library(ggplot2)
library(tidyr)
library(genomation)
library(readxl)
library(writexl)
library(xfun)
library(limma)
library(data.table)
library(openxlsx)
library(stringr)

#------------------------------------------------------------------------------
# Load metadata ----
#------------------------------------------------------------------------------
# The metadata file should contain information about samples, 
# conditions, and other experimental variables.

metadata <- read_excel("how_to_use//metadata.xlsx")

#------------------------------------------------------------------------------
# CpG analysis 
# Load preprocessed methylation files ----
#------------------------------------------------------------------------------
# Load the CpG-level methylation data (downloaded from Zenodo repository https://doi.org/10.5281/zenodo.16961293).
# Each sample is processed, filtered, and normalized.

meth <- load_methylation_data(
  repository = "/downloaded_from_zenodo/",  # Folder containing methylation data files
  sample_name_variable = "sample",           # Column name for sample identifiers
  treatment_variable = "treatment_vector",   # Column for experimental groups
  mincov = 10,                               # Minimum coverage threshold for CpGs
  individual_sample_plot = FALSE,            # Skip individual QC plots
  filter_high_percentage = TRUE,             # Filter CpGs with extreme methylation
  normalization_coverage = TRUE,             # Apply coverage normalization
  method_normalization = "median",           # Median normalization method
  metadata = metadata                        # Link metadata to samples
)

# Convert methylKit object to data frame for inspection
methdf <- getData(meth)

#------------------------------------------------------------------------------
# Visualize data distribution ----
#------------------------------------------------------------------------------
# Generate quality control plots (coverage distribution, etc.)
# This helps identify biases or outliers across samples.

analyze_methylation_stats(
  meth = methdf,
  variable_to_investigate = "coverage", 
  position_text = 100)

#------------------------------------------------------------------------------
# Extract percent methylation (beta) values ----
#------------------------------------------------------------------------------
# Convert methylation counts to percentage methylation (beta values)
# for each CpG and each sample.

betavalue <- percMethylation(meth)
betavalue_df <- as.data.frame(betavalue)
# Assign genomic coordinates as row names (chr:start-end)
location_all_cpgs <- paste0("chr", meth$chr, ":", meth$start, "-", meth$end)
rownames(betavalue_df) <- location_all_cpgs

# Visualize beta value distributions across samples
analyze_methylation_stats(
  meth = betavalue_df,
  variable_to_investigate = "betavalue", 
  position_text = 60,
  metadata = metadata,
  sample_name_variable = "sample")


#------------------------------------------------------------------------------
# Filtering CpGs ----
#------------------------------------------------------------------------------
# Filter CpGs based on differences in beta values among control groups.

groups_control <- unique(metadata$condition[metadata$sample_type == "control"])
print(groups_control)
filtered_locations <- filter_by_beta_in_cpgs(
  betavalue_df = betavalue_df,
  threshold_delta_beta = 25, # Minimum methylation difference threshold (%)
  groups = groups_control
)

# Subset CpG matrix and methylation object
idx_filtering <- which(rownames(betavalue_df)%in% filtered_locations)
betavalue_filtered <- betavalue_df[idx_filtering,]
meth_filtered <- meth[idx_filtering,]

#------------------------------------------------------------------------------
# Build genomic regions (tiles) ----
#------------------------------------------------------------------------------
# Aggregate CpGs into fixed-length genomic tiles.
# Each tile represents a genomic region with summarized CpG methylation data.

length_region = 300
tiles = tileMethylCounts(meth_filtered,win.size=length_region,step.size=length_region,cov.bases = 0)

tilesdf <- getData(tiles)

# Visualize tile coverage distribution
analyze_methylation_stats(
  meth = tilesdf,
  variable_to_investigate = "coverage", 
  position_text = 5000)

#------------------------------------------------------------------------------
# Extract beta values for regions ----
#------------------------------------------------------------------------------
# Compute average methylation percentage per tile/region.
betavalue_tiles <- percMethylation(tiles)
betavalue_tilesdf <- as.data.frame(betavalue_tiles)
location_all_regions <- paste0("chr", tilesdf$chr, ":", tilesdf$start, "-", tilesdf$end)
rownames(betavalue_tilesdf) <- location_all_regions

# Visualize regional methylation distribution
analyze_methylation_stats(
  meth = betavalue_tilesdf,
  variable_to_investigate = "betavalue", 
  position_text = 60,
  metadata = metadata,
  sample_name_variable = "sample")

#------------------------------------------------------------------------------
# Assign CpGs to regions ----
#------------------------------------------------------------------------------
# Map each CpG to its corresponding region/tile.
regions_cpgs <- map_cpg_to_regions(methdf = meth_filtered,
                                   region = tiles)
# Count number of CpGs per region before merging
N_cpg_regions <- as.data.frame(table(regions_cpgs$region_id))
N_cpg_regions_before_merging <- N_cpg_regions

# Plot distribution of CpG counts per region
do_histogram_with_statistics(
  vect = as.numeric(N_cpg_regions$Freq),
  title = "Number of Sites in Regions before merging of consecutive regions",
  x_axis = "Number of Sites",
  text_position = 80,
  color_plot = "red"
)

#------------------------------------------------------------------------------
# Optimize regions ----
#------------------------------------------------------------------------------
# Merge consecutive regions based on distance thresholds.
# This reduces fragmentation and identifies larger functional methylation blocks.
# This step can be computationally intensive for large datasets.

## Default full-data approach ----

regions_plus_merged_regions <- merging_consecutive_regions(meth_filtered, 
                                                           tiles,
                                                           da_thr = length_region,
                                                           cb_thr = length_region)

## Alternative 1: Run only on chromosome 1 for quick testing ----

idx_chr1 <- grep("chr1:", rownames(betavalue_filtered))
meth_filtered_chr1 <- meth_filtered[idx_chr1,]

idx_chr1_tiles <- grep("chr1:",rownames(betavalue_tilesdf))
tiles_chr1 <- tiles[idx_chr1_tiles,]


regions_plus_merged_regions_chr1 <- merging_consecutive_regions(meth_filtered_chr1, 
                                                                tiles_chr1,
                                                                da_thr = length_region,
                                                                cb_thr = length_region)

## Alternative 2: Run in parallel by chromosome ----
# Uses multiple CPU cores to parallelize the merging process.
# Adjust 'mc.cores' according to available computational resources.

idx_chr_cpgs <- list()
meth_subset <- list()
idx_chr_tiles <- list()
tiles_subset <- list()

chr_vector <- unique(meth_filtered$chr)

for(chrnum in chr_vector){
  print(chrnum)
  idx_chr_cpgs[[chrnum]] <- grep(paste0("chr",chrnum,":"),rownames(betavalue_filtered))
  meth_subset[[chrnum]] <- meth_filtered[idx_chr_cpgs[[chrnum]],]
  
  idx_chr_tiles[[chrnum]] <- grep(paste0("chr",chrnum,":"),rownames(betavalue_tilesdf))
  tiles_subset[[chrnum]] <- tiles[idx_chr_tiles[[chrnum]],]
}

results <- mclapply(chr_vector, function(chr) {
  meth_chunks<- meth_subset[[chr]]
  tiles_chunks <- tiles_subset[[chr]]
  merging_consecutive_regions(meth_chunks, tiles_chunks,
                              da_thr = length_region, cb_thr = length_region)
}, mc.cores = 25)
names(results) <- chr_vector

regions_plus_merged_regions_parallel <- do.call(rbind, results)

## Alternative 3: Load precomputed results ----

load("downloaded_from_zenodo/results_example_script.RData")

#------------------------------------------------------------------------------
# Assign CpGs to merged regions ----
#------------------------------------------------------------------------------

# Recompute beta values of optimized regions

betavalue_filtered <- percMethylation(meth_filtered)
location_all_cpgs_filtered <- paste0("chr", meth_filtered$chr, ":", meth_filtered$start, "-", meth_filtered$end)
betavalue_filtered <- as.data.frame(betavalue_filtered)
rownames(betavalue_filtered) <- location_all_cpgs_filtered

regions_cpgs_merged <- assign_beta_merged_regions(regions = regions_plus_merged_regions,
                                                  cpgs = meth_filtered,
                                                  betavalue = betavalue_filtered)

#------------------------------------------------------------------------------
# Compute intra-region variability (standard deviation) ----
#------------------------------------------------------------------------------
# Calculate variability across samples for each region.

groups_sd<- unique(metadata$condition)
sd_regions <- compute_beta_sd_regions(regions_cpgs = regions_cpgs_merged$regions_merged_cpg_beta,
                                      groups = groups_sd) 

#------------------------------------------------------------------------------
# CpG count per merged region ----
#------------------------------------------------------------------------------

N_cpg_regions <- as.data.frame(table(regions_cpgs_merged$regions_cpgs$region_id))
N_cpg_regions_after_merging <- N_cpg_regions

# Plot CpG count histogram (after merging)
do_histogram_with_statistics(
  vect = as.numeric(N_cpg_regions$Freq),
  title = "Number of Sites in Regions after merging of consecutive regions",
  x_axis = "Number of Sites",
  text_position = 80,
  color_plot = "red"
)


#------------------------------------------------------------------------------
# Filtering regions by CpG count and variability ----
#------------------------------------------------------------------------------

min_NCpGs = 5  # Minimum number of CpGs per region
sd_threshold = 25  # SD threshold for variability

regions_filtered_NCpGs<- as.vector(N_cpg_regions$Var1[N_cpg_regions$Freq>min_NCpGs])

high_sd_regions <- sd_regions %>%
  filter(if_any(-region_id, ~ . > sd_threshold)) %>%
  pull(region_id)

low_sd_regions <- sd_regions %>%
  filter(if_all(-region_id, ~ . <= sd_threshold)) %>%
  pull(region_id)

#------------------------------------------------------------------------------
# Visualize intra-region methylation patterns ----
#------------------------------------------------------------------------------

# Example region of interest
roi <- "chr1:852082-852381"
roi <- "chr16:2026484-2026783"
roi <- "chr10:100185886-100186185"

roidf <- regions_cpgs_merged$regions_merged_cpg_beta[regions_cpgs_merged$regions_merged_cpg_beta$region_id == roi,]

roidf$pos <-  sub(".*:(\\d+)-.*", "\\1", roidf$cpg_id)
roidf$pos <- as.numeric(roidf$pos)

df_long <- roidf %>%
  pivot_longer(
    cols = starts_with("72_C_"),
    names_to = "sample",
    values_to = "value"
  )

df_long$sample <- gsub("72_C_", "", df_long$sample)

# Plot
ggplot(df_long, aes(x = pos, y = value/100, color = sample)) +
  geom_point(size = 3) +
  geom_line(aes(group = sample), alpha = 0.5) +
  scale_color_manual(values = c("red", "blue", "green", "purple")) +
  theme_minimal() +
  labs(
    x = paste0(roi),
    y = "Beta Value",
    color = "Control Sample Replicates"
    #,title = "Example of lowSD region"
  )+
  scale_x_continuous(
    breaks = seq(min(df_long$pos), max(df_long$pos), by = 30)
  )+
  scale_y_continuous(limits = c(0, 1))

#------------------------------------------------------------------------------
# Apply region filters ----
#------------------------------------------------------------------------------
filtered_regions <- unique(intersect(regions_filtered_NCpGs, low_sd_regions))

# Convert selected regions to methylKit object if interested
filtered_regions_obj <- convert_to_methylobj(regions_location = filtered_regions, 
                                             cpgs_methyl_obj = meth_filtered)

#------------------------------------------------------------------------------
# Annotation of filtered regions ----
#------------------------------------------------------------------------------
# Retrieve and map ENSEMBL gene annotations to methylation regions.

## Download gene annotation from Ensembl BioMart ----

annotation_df <- biomart_annotation(biomart_database = "ENSEMBL_MART_ENSEMBL",
                                    biomart_dataset = "hsapiens_gene_ensembl",
                                    chromosome_vector = c(as.character(1:22), "X", "Y", "MT"),
                                    selected_attributes = c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol', 
                                                            'strand', 'transcription_start_site', 
                                                            '5_utr_start', '5_utr_end', '3_utr_start', '3_utr_end', 'exon_chrom_start', 'exon_chrom_end', 
                                                            'rank'),
                                    attributes_page = c("feature_page","structure"),
                                    gene_biotype = "chromosome_name",
                                    filter_gene_biotype = c(as.character(1:22), "X", "Y", "MT"), 
                                    promoter_distance = c(500,100))
## Map methylation regions to annotated genes ----
annotated <- map_regions(sample_to_map =  regions_plus_merged_regions, annotation = annotation_df, origin_annotation = "gene")

#------------------------------------------------------------------------------
# Visualization of annotation overlap ----
#------------------------------------------------------------------------------
# Count overlap between all methylation regions and annotated genes.

library(ggplot2)
all_regions <- unique(regions_plus_merged_regions$region_id)
annotated_regions <- unique(annotated$region_id)
overlap <- length(intersect(all_regions, annotated_regions))
non_overlap <- length(setdiff(all_regions, annotated_regions))

df <- data.frame(
  category = c("Overlap", "No Overlap"),
  count = c(overlap, non_overlap)
)

# Bar plot of overlap counts

ggplot(df, aes(x = category, y = count, fill = category)) +
  geom_col() +
  scale_fill_manual(values = c("Overlap" = "#2F3E9E", "No Overlap" = "#F25C54")) +
  geom_text(aes(label = count), vjust = -0.5) +
  labs(title = "Overlap of regions with annotated data",
       y = "Number of regions",
       x = "") +
  theme_minimal()

# Visualization of overlap percentage

ggplot(annotated, aes(x = pct_merged)) +
  geom_histogram(fill = "#2F3E9E", bins = 30) +
  labs(title = "Distribution %overlap",
       x = "% overlap",
       y = "Freq") +
  theme_minimal()

#------------------------------------------------------------------------------
# Generate heatmap of gene × region relationships ----
#------------------------------------------------------------------------------
# Simplify gene location categories

annotated[which(!(annotated$gene_loci %in% c("3UTR", "5UTR", "promoter","entire_gene"))), "gene_loci"] <- "gene_body"
annotated <- as.data.frame(annotated)
# Count occurrences per gene-region combination
count_long <- annotated %>%
  dplyr::count(gene, gene_loci)%>%
  complete(gene, gene_loci, fill = list(n = 0))

count_long <- count_long[order(count_long$n, decreasing = T),]
sample_genes <- sample(count_long$gene,size = 10)

# Heatmap visualization
p<-ggplot(count_long[count_long$gene %in% sample_genes,], aes(x = gene_loci, y = gene, fill = n)) +
  geom_tile(color = "grey90") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 10, face = "bold")
  ) +
  labs(
    title = "Gene × Region Heatmap",
    x = "Genomic Region",
    y = "Gene",
    fill = "Count"
  ) 
plot(p)

save.image(file = "results_example_script.RData")

