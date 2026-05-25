library(tidyverse)
library(readxl)
library(ggplot2)
library(dplyr)
library(vegan)
library(patchwork)

meta_per <- read_excel("meta_2025_v7.xlsx", sheet = "R") %>%
  dplyr::select(
    SN, rawseqID, sample_code, subject_code, household, reads, timepoint, carrier, st131_detect, age,
    sex, abx_6months, index_pt, st131_qpcr_trpa_pabb, st131_mnth, st131_mnth_isolate_count, st131_wgs,
    st131pos_density_wgs, raw_counts, old_SN
  ) %>%
  filter(carrier %in% c("non_carrier", "per", "int"), !is.na(carrier)) %>%
  filter(st131_detect %in% c("yes", "no"), !is.na(st131_detect)) %>%
  mutate(
    st131_detect = factor(st131_detect, levels = c("no", "yes")),
    index_pt = factor(index_pt, levels = c("yes", "no")),
    carrier = factor(carrier, levels = c("non_carrier", "int", "per")),
    abx_6months = factor(abx_6months, levels = c("yes", "no")),
    timepoint = factor(timepoint, levels = c("Baseline", "Repeat")),
    sex = factor(sex, levels = c("M", "F")),
    subject_code = factor(subject_code),
    household = factor(household),
    Kneadcount_mil = reads / 1000000
  ) # Convert reads to millions

# Load EC data
ec_raw <- read_excel("meta_2025_v7.xlsx", sheet = "ec") %>%
  gather(SN, relab, -EC) %>%
  spread(EC, relab) %>%
  column_to_rownames(var = "SN")

ec_ready_to_join <- ec_raw %>%
  rownames_to_column(var = "SN")

ec_left_joined_df <- left_join(ec_ready_to_join, meta_per, by = "SN") %>%
  filter(carrier %in% c("non_carrier", "per", "int"), !is.na(carrier))

# Identify EC columns and exclude UNMAPPED/UNGROUPED from sum calculation
ec_cols <- grep("^EC_", colnames(ec_left_joined_df), value = TRUE)
# Exclude UNMAPPED and UNGROUPED if they exist
ec_cols_for_sum <- ec_cols[!ec_cols %in% c("UNMAPPED", "UNGROUPED", "EC_UNMAPPED", "EC_UNGROUPED")]

# Calculate sum of EC values per sample (excluding UNMAPPED and UNGROUPED)
ec_left_joined_df <- ec_left_joined_df %>%
  rowwise() %>%
  mutate(ec_sum = sum(c_across(all_of(ec_cols_for_sum)), na.rm = TRUE)) %>%
  ungroup()

# Normalize EC values: (value / sum) * 1,000,000, then add 1 and log2 transform
normalized_ec_df <- ec_left_joined_df %>%
  mutate(
    across(
      all_of(ec_cols_for_sum),
      ~ {
        # Avoid division by zero
        if_else(ec_sum > 0,
          log2((.x / ec_sum * 1000000) + 1),
          0
        )
      }
    )
  ) %>%
  # Remove the temporary sum column
  select(-ec_sum)

# normalized_ec_df already contains all metadata columns from the initial join
# No need to join again with meta_per
normalized_ec_df2 <- normalized_ec_df

# Prepare data for plotting
# For ST131 analysis: use both baseline and repeat samples
df_ec_boxplot_st131 <- normalized_ec_df2 %>%
  pivot_longer(
    cols = starts_with("EC_"),
    names_to = "EC",
    values_to = "Count"
  ) %>%
  na.omit()

# For carrier analysis: use only baseline samples
df_ec_boxplot_carrier <- normalized_ec_df2 %>%
  filter(timepoint == "Baseline") %>%
  pivot_longer(
    cols = starts_with("EC_"),
    names_to = "EC",
    values_to = "Count"
  ) %>%
  na.omit()

selected_ec_columns <- c("EC__1.14.13.59", "EC__2.3.1.102", "EC__6.3.2.38", "EC__6.3.2.39")

# For ST131 detection analysis (both baseline and repeat)
filtered_long_df_st131 <- df_ec_boxplot_st131 %>%
  filter(EC %in% selected_ec_columns) %>%
  mutate(
    # Create display names for EC genes
    EC_display = case_when(
      EC == "EC__1.14.13.59" ~ "iucD, EC__1.14.13.59",
      EC == "EC__2.3.1.102" ~ "iucB, EC__2.3.1.102",
      EC == "EC__6.3.2.38" ~ "iucA, EC__6.3.2.38",
      EC == "EC__6.3.2.39" ~ "iucC, EC__6.3.2.39",
      TRUE ~ EC
    ),
    # Set factor levels to control order
    EC_display = factor(EC_display, levels = c(
      "iucD, EC__1.14.13.59",
      "iucB, EC__2.3.1.102",
      "iucA, EC__6.3.2.38",
      "iucC, EC__6.3.2.39"
    ))
  )

# For carrier status analysis (baseline only)
filtered_long_df_carrier <- df_ec_boxplot_carrier %>%
  filter(EC %in% selected_ec_columns) %>%
  mutate(
    # Create display names for EC genes
    EC_display = case_when(
      EC == "EC__1.14.13.59" ~ "iucD, EC__1.14.13.59",
      EC == "EC__2.3.1.102" ~ "iucB, EC__2.3.1.102",
      EC == "EC__6.3.2.38" ~ "iucA, EC__6.3.2.38",
      EC == "EC__6.3.2.39" ~ "iucC, EC__6.3.2.39",
      TRUE ~ EC
    ),
    # Set factor levels to control order
    EC_display = factor(EC_display, levels = c(
      "iucD, EC__1.14.13.59",
      "iucB, EC__2.3.1.102",
      "iucA, EC__6.3.2.38",
      "iucC, EC__6.3.2.39"
    ))
  )

# ==========================================
# SPEARMAN CORRELATION: EC vs E. coli abundance
# ==========================================

# Load taxa data using the same method as linear_regression_univariate.R
taxa_model <- read_excel("meta_2025_v7.xlsx", sheet = "taxa") %>%
  select(-Kingdom, -Phyla, -Class, -Order, -Family, -Genus) %>%
  gather(rawseqID, relab, -Species) %>%
  mutate(relab = relab / 100) %>%
  spread(Species, relab) %>%
  filter(rawseqID %in% c(meta_per$rawseqID))

# Replace NAs with 0
taxa_model[is.na(taxa_model)] <- 0

# Convert to matrix and ensure numeric
species_mat_model <- taxa_model %>% column_to_rownames(var = "rawseqID")
species_mat_model <- as.data.frame(species_mat_model)
species_mat_model[] <- lapply(species_mat_model, function(x) as.numeric(x))
species_mat_model <- species_mat_model[rowSums(species_mat_model, na.rm = TRUE) > 0, , drop = FALSE]

# Calculate relative abundance using decostand
species_relab_model <- decostand(species_mat_model, MARGIN = 1, method = "total")

# Filter species present in >5% of samples
bugs_f_model <- species_relab_model[, which(apply(species_relab_model, 2, function(x) length(which(x > 0)) / length(x)) > 0.05), drop = FALSE]

# Replace zeros with half of minimum non-zero value
bugs_ff_model <- apply(bugs_f_model, 2, function(x) {
  x[x == 0] <- min(x[x > 0]) / 2
  x
})

# Log10 transformation
bugs_log10_model <- log10(bugs_ff_model)

# Convert back to tibble with rawseqID
mgx.log10.model <- bugs_log10_model %>% as_tibble(rownames = "rawseqID")

# Get E. coli abundance (log10 transformed) and merge with metadata
ecoli_abundance <- mgx.log10.model %>%
  select(rawseqID, s__Escherichia_coli) %>%
  inner_join(meta_per %>% select(rawseqID, SN, carrier), by = "rawseqID")

# Merge E. coli abundance with EC data
ec_ecoli_data <- normalized_ec_df2 %>%
  select(SN, carrier, all_of(selected_ec_columns)) %>%
  inner_join(ecoli_abundance %>% select(SN, s__Escherichia_coli, carrier), by = c("SN", "carrier")) %>%
  pivot_longer(
    cols = all_of(selected_ec_columns),
    names_to = "EC",
    values_to = "EC_TPM"
  ) %>%
  mutate(
    EC_display = case_when(
      EC == "EC__1.14.13.59" ~ "iucD, EC__1.14.13.59",
      EC == "EC__2.3.1.102" ~ "iucB, EC__2.3.1.102",
      EC == "EC__6.3.2.38" ~ "iucA, EC__6.3.2.38",
      EC == "EC__6.3.2.39" ~ "iucC, EC__6.3.2.39",
      TRUE ~ EC
    ),
    EC_display = factor(EC_display, levels = c(
      "iucD, EC__1.14.13.59",
      "iucB, EC__2.3.1.102",
      "iucA, EC__6.3.2.38",
      "iucC, EC__6.3.2.39"
    )),
    carrier = factor(carrier, levels = c("non_carrier", "int", "per"))
  ) %>%
  na.omit()

# Calculate Spearman correlations for each EC and carrier group
spearman_results <- ec_ecoli_data %>%
  group_by(EC_display, carrier) %>%
  summarise(
    n = n(),
    rho = cor(s__Escherichia_coli, EC_TPM, method = "spearman"),
    p_value = cor.test(s__Escherichia_coli, EC_TPM, method = "spearman")$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    label = sprintf("ρ = %.3f %s", rho, significance)
  )

cat("\nSpearman correlation results:\n")
print(spearman_results)

# Save correlation results
write.csv(spearman_results, "ec_ecoli_spearman_correlation.csv", row.names = FALSE)

# Create scatter plots with fitted lines for each EC
color_palette <- c("non_carrier" = "#2a9df4", "int" = "#accc87", "per" = "#b03535")

# Function to create individual plots (reversed axes: x = EC, y = E. coli)
create_ec_plot <- function(ec_name) {
  plot_data <- ec_ecoli_data %>% filter(EC_display == ec_name)
  corr_data <- spearman_results %>% filter(EC_display == ec_name)

  p <- ggplot(plot_data, aes(x = EC_TPM, y = s__Escherichia_coli, color = carrier)) +
    geom_point(size = 1, alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.2, linewidth = 1) +
    scale_color_manual(
      values = color_palette,
      name = "Carrier Status",
      labels = c("non_carrier" = "Non-carrier", "int" = "Intermittent", "per" = "Persistent")
    ) +
    labs(
      title = ec_name,
      x = "E.C. Normalized Abundance (Log2 TPM)",
      y = "E.coli abundance (log10)"
    ) +
    # Zoom: EC on x ≥ 0, E. coli y ≥ -3.6
    #    coord_cartesian(xlim = c(0, NA), ylim = c(-3.6, NA)) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )

  # Add correlation annotations (anchor to the top-right *corner* of the panel)
  for (i in seq_len(nrow(corr_data))) {
    carrier_label <- corr_data$carrier[i]
    label_text <- corr_data$label[i]

    p <- p + annotate(
      "text",
      x = Inf,
      y = Inf,
      label = paste0(carrier_label, ": ", label_text),
      color = color_palette[as.character(carrier_label)],
      hjust = 1.1,
      vjust = 1.5 + (i - 1) * 2.0, # start a bit lower + keep larger spacing
      size = 3
    )
  }

  return(p)
}

# Create plots for each EC
ec_plots <- list()
for (ec in levels(ec_ecoli_data$EC_display)) {
  ec_plots[[ec]] <- create_ec_plot(ec)
}

# Combine plots using patchwork
combined_plot <- wrap_plots(ec_plots, ncol = 2)

print(combined_plot)

ggsave("ec_ecoli_correlation_combined.svg", combined_plot, width = 14, height = 12, dpi = 300)

# Save individual plots
# for (ec in names(ec_plots)) {
#  ec_clean <- gsub("[^[:alnum:]_]", "_", ec)
#  ggsave(paste0("ec_ecoli_correlation_", ec_clean, ".png"),
#       ec_plots[[ec]], width = 8, height = 6, dpi = 300)
#  ggsave(paste0("ec_ecoli_correlation_", ec_clean, ".svg"),
#         ec_plots[[ec]], width = 8, height = 6, dpi = 300)
# }
