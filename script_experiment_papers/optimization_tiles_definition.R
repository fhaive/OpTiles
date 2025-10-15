# Optimization tiles definition

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

# Load data ----
## Load metadata ----
metadata <- read_excel("how_to_use/metadata.xlsx")
# CpG analysis 

## Load preprocessed files ----
meth <- load_methylation_data(
  repository = "/downloaded_from_zenodo/", #folder where there are the downloaded file from zenodo https://doi.org/10.5281/zenodo.16961293
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

## Extract percent methylation values (beta values)
betavalue <- percMethylation(meth)
betavalue_df <- as.data.frame(betavalue)
location_all_cpgs <- paste0("chr", meth$chr, ":", meth$start, "-", meth$end)
rownames(betavalue_df) <- location_all_cpgs

# Filtering CpGs ----
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

# Extract percent methylation values (beta values)
betavalue_tiles <- percMethylation(tiles)
betavalue_tilesdf <- as.data.frame(betavalue_tiles)
location_all_regions <- paste0("chr", tilesdf$chr, ":", tilesdf$start, "-", tilesdf$end)
rownames(betavalue_tilesdf) <- location_all_regions

## Assign to each region the CpGs
regions_cpgs <- map_cpg_to_regions(methdf = meth_filtered,
                                   region = tiles)
## Compute the number of CpG sites per region
N_cpg_regions <- as.data.frame(table(regions_cpgs$region_id))
N_cpg_regions_before_merging <- N_cpg_regions

#Compute sd tiles
tiles2 = regions_cpgs
tiles_location <- data.frame(
  chr = sub("^chr([^:]+):.*", "\\1", tiles2$region_id),
  start = as.integer(sub("^[^:]+:(\\d+)-.*", "\\1", tiles2$region_id)),
  end = as.integer(sub("^.*-(\\d+)$", "\\1", tiles2$region_id)),
  region_id = tiles2$region_id
)
betavalue <- percMethylation(meth_filtered)
betavalue_df <- as.data.frame(betavalue)
location_all_cpgs <- paste0("chr", meth_filtered$chr, ":", meth_filtered$start, "-", meth_filtered$end)
rownames(betavalue_df) <- location_all_cpgs

tiles_cpgs_beta <- assign_beta_merged_regions(regions = tiles_location,
                                              cpgs = meth_filtered,
                                              betavalue = betavalue_df)
groups_sd<- unique(metadata$condition)
sd_tiles <- compute_beta_sd_regions(regions_cpgs = tiles_cpgs_beta$regions_merged_cpg_beta, groups = groups_sd)

# Optimize regions ----
library(tictoc)
# Consider to parallelize this step as shown in the example_script.R
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


# Comparison tiles with merged regions ----
N_cpg_regions_tiles <- N_cpg_regions_before_merging
N_cpg_regions_tiles$Var1 <- as.vector(N_cpg_regions_tiles$Var1)
NCpGs_tiles <- data.frame(region_id = as.vector(N_cpg_regions_tiles$Var1),
                          Freq = as.vector(N_cpg_regions_tiles$Freq))
colnames(NCpGs_tiles) <- colnames(sd_tiles)
NCpGs_tiles <- NCpGs_tiles[which(NCpGs_tiles$region_id %in% sd_tiles$region_id),]

N_cpg_regions_merged <- N_cpg_regions_after_merging
N_cpg_regions_merged$Var1 <- as.vector(N_cpg_regions_after_merging$Var1)
NCpGs_merged <- data.frame(region_id = as.vector(N_cpg_regions_merged$Var1),
                           Freq = as.vector(N_cpg_regions_merged$Freq))
colnames(NCpGs_merged) <- colnames(sd_regions)
NCpGs_merged <- NCpGs_merged[which(NCpGs_merged$region_id %in% sd_regions$region_id),]

transform_in_vector <- function(df){
  
  vector <- df
  if ("region_id" %in% colnames(vector)){
    vector$region_id <- NULL
  }
  vector <- vector %>%
    pivot_longer(cols = everything(), values_to = "value") %>%
    pull(value)
  
  return(vector)
}

sd_regions <- as.data.frame(sd_regions)
sd_tiles <- as.data.frame(sd_tiles)

sd_vect_merged <- transform_in_vector(sd_regions)
NCpGs_vect_merged <- transform_in_vector(NCpGs_merged)

sd_vect_tile <- transform_in_vector(sd_tiles)
NCpGs_vect_tiles <- transform_in_vector(NCpGs_tiles)

plol_hex <- function(sd,ncpg,title){
  
  df_long <- pivot_longer(sd,
                          cols = -region_id,
                          names_to = "sample",
                          values_to = "sd") 
  
  
  df_long$NCpG <- ncpg$Freq[match(df_long$region_id, ncpg$Var1)]
  

  ggplot(df_long, aes(x = NCpG, y = sd)) +
    geom_density_2d_filled(alpha = 0.8) +
    labs(title = title,
         x = "NCpGs/region",
         y = "std/region", fill = "Density") +
    theme_minimal()+
    xlim(0,30) +
    ylim(0,20) +
    scale_fill_viridis_d(name = "Estimated Density\n(kernel 2D)")+
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = 6),axis.text.y = element_text(size = 6))

  
}
plol_hex(sd = sd_tiles[,c("72_C","region_id")],N_cpg_regions_tiles,"Tiling-genome")
plol_hex(sd = sd_regions[,c("72_C","region_id")],N_cpg_regions_merged,"OpTiles")


save.image(file = "script_experiment_papers/results_optimization_tiles_definition.RData")
