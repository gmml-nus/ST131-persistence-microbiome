library(tidyverse)
library(readxl)
library(ggplot2)
library(dplyr)
library(vegan)
library(patchwork)
library(rstatix)

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
  ) 

# load EC data
ec_raw <- read_excel("meta_2025_v7.xlsx", sheet = "ec") %>%
  gather(SN, relab, -EC) %>%
  spread(EC, relab) %>%
  column_to_rownames(var = "SN")

ec_ready_to_join <- ec_raw %>%
  rownames_to_column(var = "SN")
ec_left_joined_df <- left_join(ec_ready_to_join, meta_per, by = "SN") %>%
  filter(carrier %in% c("non_carrier", "per", "int"), !is.na(carrier))
ec_cols <- grep("^EC_", colnames(ec_left_joined_df), value = TRUE)
ec_cols_for_sum <- ec_cols[!ec_cols %in% c("UNMAPPED", "UNGROUPED", "EC_UNMAPPED", "EC_UNGROUPED")]

ec_left_joined_df <- ec_left_joined_df %>%
  rowwise() %>%
  mutate(ec_sum = sum(c_across(all_of(ec_cols_for_sum)), na.rm = TRUE)) %>%
  ungroup()

normalized_ec_df <- ec_left_joined_df %>%
  mutate(
    across(
      all_of(ec_cols_for_sum),
      ~ {
        if_else(ec_sum > 0,
          log2((.x / ec_sum * 1000000) + 1),
          0
        )
      }
    )
  ) %>%
  select(-ec_sum)

# normalized_ec_df already contains all metadata columns from the initial join
normalized_ec_df2 <- normalized_ec_df

df_ec_boxplot_st131 <- normalized_ec_df2 %>%
  pivot_longer(
    cols = starts_with("EC_"),
    names_to = "EC",
    values_to = "Count"
  ) %>%
  na.omit()

df_ec_boxplot_carrier <- normalized_ec_df2 %>%
  filter(timepoint == "Baseline") %>%
  pivot_longer(
    cols = starts_with("EC_"),
    names_to = "EC",
    values_to = "Count"
  ) %>%
  na.omit()

selected_ec_columns <- c("EC__1.14.13.59", "EC__2.3.1.102", "EC__6.3.2.38", "EC__6.3.2.39")

filtered_long_df_st131 <- df_ec_boxplot_st131 %>%
  filter(EC %in% selected_ec_columns) %>%
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
    ))
  )

filtered_long_df_carrier <- df_ec_boxplot_carrier %>%
  filter(EC %in% selected_ec_columns) %>%
  mutate(
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

# Load taxa data using the same method as linear_regression_univariate.R
taxa_model <- read_excel("meta_2025_v7.xlsx", sheet = "taxa") %>%
  select(-Kingdom, -Phyla, -Class, -Order, -Family, -Genus) %>%
  gather(rawseqID, relab, -Species) %>%
  mutate(relab = relab / 100) %>%
  spread(Species, relab) %>%
  filter(rawseqID %in% c(meta_per$rawseqID))

taxa_model[is.na(taxa_model)] <- 0
species_mat_model <- taxa_model %>% column_to_rownames(var = "rawseqID")
species_mat_model <- as.data.frame(species_mat_model)
species_mat_model[] <- lapply(species_mat_model, function(x) as.numeric(x))
species_mat_model <- species_mat_model[rowSums(species_mat_model, na.rm = TRUE) > 0, , drop = FALSE]
species_relab_model <- decostand(species_mat_model, MARGIN = 1, method = "total")

bugs_f_model <- species_relab_model[, which(apply(species_relab_model, 2, function(x) length(which(x > 0)) / length(x)) > 0.05), drop = FALSE]

bugs_ff_model <- apply(bugs_f_model, 2, function(x) {
  x[x == 0] <- min(x[x > 0]) / 2
  x
})
bugs_log10_model <- log10(bugs_ff_model)
mgx.log10.model <- bugs_log10_model %>% as_tibble(rownames = "rawseqID")
ecoli_abundance <- mgx.log10.model %>%
  select(rawseqID, s__Escherichia_coli) %>%
  inner_join(meta_per %>% select(rawseqID, SN, carrier), by = "rawseqID")

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

color_palette <- c("non_carrier" = "#2a9df4", "int" = "#accc87", "per" = "#b03535")

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
      vjust = 1.5 + (i - 1) * 2.0, 
      size = 3
    )
  }

  return(p)
}

ec_plots <- list()
for (ec in levels(ec_ecoli_data$EC_display)) {
  ec_plots[[ec]] <- create_ec_plot(ec)
}

combined_plot <- wrap_plots(ec_plots, ncol = 2)

ggsave("ec_ecoli_correlation_combined.svg", combined_plot, width = 14, height = 12, dpi = 300)

# Export p/q values
spearman_results %>%
  mutate(q_value_BH = p.adjust(p_value, method = "BH")) %>%
  write.csv("ec_ecoli_spearman_correlation.csv", row.names = FALSE)

st131_boxplot_stats <- filtered_long_df_st131 %>%
  group_by(EC, EC_display) %>%
  wilcox_test(Count ~ st131_detect) %>%
  add_significance("p") %>%
  ungroup() %>%
  transmute(
    EC, EC_display,
    analysis = "ST131 detection", test = "Wilcoxon",
    Group1 = group1, Group2 = group2,
    p_value = p, q_value_BH = NA_real_, Significance = p.signif
  )

carrier_kw_stats <- filtered_long_df_carrier %>%
  group_by(EC, EC_display) %>%
  kruskal_test(Count ~ carrier) %>%
  add_significance("p") %>%
  ungroup() %>%
  transmute(
    EC, EC_display,
    analysis = "Carrier status", test = "Kruskal-Wallis",
    Group1 = "All", Group2 = "All",
    p_value = p, q_value_BH = NA_real_, Significance = p.signif
  )

carrier_dunn_stats <- filtered_long_df_carrier %>%
  group_by(EC, EC_display) %>%
  dunn_test(Count ~ carrier, p.adjust.method = "BH") %>%
  add_significance("p.adj") %>%
  ungroup() %>%
  transmute(
    EC, EC_display,
    analysis = "Carrier status", test = "Dunn (post-hoc)",
    Group1 = group1, Group2 = group2,
    p_value = p, q_value_BH = p.adj, Significance = p.adj.signif
  )

bind_rows(st131_boxplot_stats, carrier_kw_stats, carrier_dunn_stats) %>%
  write.csv("ec_boxplot_pvalues.csv", row.names = FALSE)
