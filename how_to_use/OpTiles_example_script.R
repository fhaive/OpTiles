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

# Load metadata ----
metadata <- read_excel("/project/qkgimi/optiles/data/preprocessed_files//metadata.xlsx")

# CpG analysis 
## Load preprocessed files ----
meth <- load_methylation_data(
  repository = "/project/qkgimi/optiles/data/preprocessed_files//",
  sample_name_variable = "sample",
  treatment_variable = "treatment_vector",
  mincov = 10,
  individual_sample_plot = FALSE, #suggested to give a look to the data
  filter_high_percentage = TRUE,
  normalization_coverage = TRUE,
  method_normalization = "median",
  metadata = metadata
)

methdf <- getData(meth)

## Visualize data distribution ----
analyze_methylation_stats(
  meth = methdf,
  variable_to_investigate = "coverage", 
  position_text = 100)

## Extract percent methylation values (beta values)
betavalue <- percMethylation(meth)
betavalue_df <- as.data.frame(betavalue)
location_all_cpgs <- paste0("chr", meth$chr, ":", meth$start, "-", meth$end)
rownames(betavalue_df) <- location_all_cpgs

analyze_methylation_stats(
  meth = betavalue_df,
  variable_to_investigate = "betavalue", 
  position_text = 60,
  metadata = metadata,
  sample_name_variable = "sample")


## Filtering CpGs ----
groups_control <- unique(metadata$condition[metadata$sample_type == "control"])
print(groups_control)
filtered_locations <- filter_by_beta_in_cpgs(
  betavalue_df = betavalue_df,
  threshold_delta_beta = 25,
  groups = groups_control
)

idx_filtering <- which(rownames(betavalue_df)%in% filtered_locations)

betavalue_filtered <- betavalue_df[idx_filtering,]
meth_filtered <- meth[idx_filtering,]

# Build regions ----
length_region = 300
tiles = tileMethylCounts(meth_filtered,win.size=length_region,step.size=length_region,cov.bases = 0)

tilesdf <- getData(tiles)

## Visualize data distribution ----
analyze_methylation_stats(
  meth = tilesdf,
  variable_to_investigate = "coverage", 
  position_text = 5000)

# Extract percent methylation values (beta values)
betavalue_tiles <- percMethylation(tiles)
betavalue_tilesdf <- as.data.frame(betavalue_tiles)
location_all_regions <- paste0("chr", tilesdf$chr, ":", tilesdf$start, "-", tilesdf$end)
rownames(betavalue_tilesdf) <- location_all_regions

analyze_methylation_stats(
  meth = betavalue_tilesdf,
  variable_to_investigate = "betavalue", 
  position_text = 60,
  metadata = metadata,
  sample_name_variable = "sample")

## Assign to each region the CpGs
regions_cpgs <- map_cpg_to_regions(methdf = meth_filtered,
                                   region = tiles)
## Compute the number of CpG sites per region
N_cpg_regions <- as.data.frame(table(regions_cpgs$region_id))
N_cpg_regions_before_merging <- N_cpg_regions

# Load and execute the histogram function
do_histogram_with_statistics(
  vect = as.numeric(N_cpg_regions$Freq),
  title = "Number of Sites in Regions before merging of consecutive regions",
  x_axis = "Number of Sites",
  text_position = 80,
  color_plot = "red"
)

#Compute sd regions
# tiles2 = regions_cpgs
# tiles_location <- data.frame(
#   chr = sub("^chr([^:]+):.*", "\\1", tiles2$region_id),
#   start = as.integer(sub("^[^:]+:(\\d+)-.*", "\\1", tiles2$region_id)),
#   end = as.integer(sub("^.*-(\\d+)$", "\\1", tiles2$region_id)),
#   region_id = tiles2$region_id
# )
# betavalue <- percMethylation(meth_filtered)
# betavalue_df <- as.data.frame(betavalue)
# location_all_cpgs <- paste0("chr", meth_filtered$chr, ":", meth_filtered$start, "-", meth_filtered$end)
# rownames(betavalue_df) <- location_all_cpgs
# 
# tiles_cpgs_beta <- assign_beta_merged_regions(regions = tiles_location,
#                                               cpgs = meth_filtered,
#                                               betavalue = betavalue_df)
# groups_sd<- unique(metadata$condition)
# sd_tiles <- compute_beta_sd_regions(regions_cpgs = tiles_cpgs_beta$regions_merged_cpg_beta, groups = groups_sd)

# Optimize regions ----
regions_plus_merged_regions <- merging_consecutive_regions(meth_filtered, 
                                                           tiles,
                                                           da_thr = length_region,
                                                           cb_thr = length_region)

betavalue_filtered <- percMethylation(meth_filtered)
location_all_cpgs_filtered <- paste0("chr", meth_filtered$chr, ":", meth_filtered$start, "-", meth_filtered$end)
betavalue_filtered <- as.data.frame(betavalue_filtered)
rownames(betavalue_filtered) <- location_all_cpgs_filtered

regions_cpgs_merged <- assign_beta_merged_regions(regions = regions_plus_merged_regions,
                                                  cpgs = meth_filtered,
                                                  betavalue = betavalue_filtered)

#if you want to investigate the standard deviation within region across the groups
groups_sd<- unique(metadata$condition)
sd_regions <- compute_beta_sd_regions(regions_cpgs = regions_cpgs_merged$regions_merged_cpg_beta,
                                      groups = groups_sd) 

## Compute the number of CpG sites per region
N_cpg_regions <- as.data.frame(table(regions_cpgs_merged$regions_cpgs$region_id))
N_cpg_regions_after_merging <- N_cpg_regions

# Load and execute the histogram function
do_histogram_with_statistics(
  vect = as.numeric(N_cpg_regions$Freq),
  title = "Number of Sites in Regions after merging of consecutive regions",
  x_axis = "Number of Sites",
  text_position = 80,
  color_plot = "red"
)

#filter for minimum of CpGs
min_NCpGs = 5
regions_filtered_NCpGs<- as.vector(N_cpg_regions$Var1[N_cpg_regions$Freq>min_NCpGs])

#filter lowSD
sd_threshold = 25

high_sd_regions <- sd_regions %>%
  filter(if_any(-region_id, ~ . > sd_threshold)) %>%
  pull(region_id)

low_sd_regions <- sd_regions %>%
  filter(if_all(-region_id, ~ . <= sd_threshold)) %>%
  pull(region_id)

## Regions intra-variability

roi = "chr1:852082-852381"
roi <- "chr16:2026484-2026783"
roi <- "chr10:100185886-100186185"

roidf <- regions_cpgs_merged$regions_merged_cpg_beta[regions_cpgs_merged$regions_merged_cpg_beta$region_id == roi,]

roidf$pos <-  sub(".*:(\\d+)-.*", "\\1", roidf$cpg_id)

library(ggplot2)
library(tidyr)

# Esempio: supponiamo che il tuo dataframe si chiami df
roidf$pos <- as.numeric(roidf$pos)

df_long <- roidf %>%
  pivot_longer(
    cols = starts_with("72_C_"),
    names_to = "sample",
    values_to = "value"
  )

# Rimuovo "72_C_" e lascio solo 1, 2, 3, 4 come fattore
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

## apply filters----
filtered_regions <- unique(intersect(regions_filtered_NCpGs, low_sd_regions))

#convert to methylkitobject
filtered_regions_obj <- convert_to_methylobj(regions_location = filtered_regions, 
                                             cpgs_methyl_obj = meth_filtered)

# Annotate ----

## Download annotation of interest ----
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
## Map annotation
annotated <- map_regions(sample_to_map =  regions_plus_merged_regions, annotation = annotation_df, origin_annotation = "gene")

library(ggplot2)
all_regions <- unique(regions_plus_merged_regions$region_id)
annotated_regions <- unique(annotated$region_id)
overlap <- length(intersect(all_regions, annotated_regions))
non_overlap <- length(setdiff(all_regions, annotated_regions))

df <- data.frame(
  category = c("Overlap", "No Overlap"),
  count = c(overlap, non_overlap)
)

ggplot(df, aes(x = category, y = count, fill = category)) +
  geom_col() +
  scale_fill_manual(values = c("Overlap" = "#2F3E9E", "No Overlap" = "#F25C54")) +
  geom_text(aes(label = count), vjust = -0.5) +
  labs(title = "Overlap of regions with annotated data",
       y = "Number of regions",
       x = "") +
  theme_minimal()


ggplot(annotated, aes(x = pct_merged)) +
  geom_histogram(fill = "#2F3E9E", bins = 30) +
  labs(title = "Distribution %overlap",
       x = "% overlap",
       y = "Freq") +
  theme_minimal()
#Example of plot

annotated[which(!(annotated$gene_loci %in% c("3UTR", "5UTR", "promoter","entire_gene"))), "gene_loci"] <- "gene_body"
annotated <- as.data.frame(annotated)
# Count occurrencies
count_long <- annotated %>%
  dplyr::count(gene, gene_loci)%>%
  complete(gene, gene_loci, fill = list(n = 0))

count_long <- count_long[order(count_long$n, decreasing = T),]
sample_genes <- sample(count_long$gene,size = 10)

# Heatmap
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
save.image(file ="OpTiles/results_example_script.RData")
