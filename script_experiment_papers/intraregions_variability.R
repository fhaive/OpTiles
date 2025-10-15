#Intra-regions methylation variability

#This script continue from the results obtained from "optimization_tiles_definition.R" script.

load("script_experiment_papers/results_optimization_tiles_definition.RData")

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

## Regions intra-variability, do a plot for each region of interest.

roi = "chr1:852082-852381"
roi <- "chr16:2026484-2026783"
roi <- "chr10:100185886-100186185"

roidf <- regions_cpgs_merged$regions_merged_cpg_beta[regions_cpgs_merged$regions_merged_cpg_beta$region_id == roi,]

roidf$pos <-  sub(".*:(\\d+)-.*", "\\1", roidf$cpg_id)

library(ggplot2)
library(tidyr)

#
roidf$pos <- as.numeric(roidf$pos)

df_long <- roidf %>%
  pivot_longer(
    cols = starts_with("72_C_"),
    names_to = "sample",
    values_to = "value"
  )

# remove the replicates from the name
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
