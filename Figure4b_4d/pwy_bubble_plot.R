library(tidyverse)
library(readxl)
library(ggplot2)
library(viridis)
library(scales)
library(readr)

# Load metadata
meta <- read_excel("metadata.xlsx", sheet = "R") %>%
  dplyr::select(
    SN, rawseqID, sample_code, subject_code, household, reads, timepoint, carrier, st131_detect, age,
    sex, abx_6months, index_pt, st131_qpcr_trpa_pabb, st131_mnth, st131_mnth_isolate_count, st131_wgs,
    st131pos_density_wgs, raw_counts, old_SN
  ) %>%
  mutate(
    st131_detect = factor(st131_detect, levels = c("yes", "no", "NA")),
    index_pt = factor(index_pt, levels = c("yes", "no")),
    carrier = factor(carrier, levels = c("non_carrier", "per", "int")),
    abx_6months = factor(abx_6months, levels = c("yes", "no")),
    timepoint = factor(timepoint, levels = c("Baseline", "Repeat")),
    sex = factor(sex, levels = c("M", "F")),
    subject_code = factor(subject_code),
    household = factor(household)
  )

meta_filtered <- meta %>%
  filter(carrier %in% c("non_carrier", "per", "int"), !is.na(carrier)) %>%
  filter(st131_detect %in% c("yes", "no"))

pwy_raw <- read_excel("metadata.xlsx", sheet = "pwy")
pwy_abundance <- pwy_raw %>%
  column_to_rownames("PWY") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("SN")

pwy_cols <- setdiff(colnames(pwy_abundance), c("SN", "UNMAPPED", "UNINTEGRATED"))
pwy_abundance <- pwy_abundance %>% select(SN, all_of(pwy_cols))

pwy_map <- read_csv("pwy.csv", show_col_types = FALSE) %>%
  distinct(PWY, pwy)
pwy_list <- pwy_map %>% pull(PWY)

available_pwys <- intersect(pwy_list, colnames(pwy_abundance))
pwy_data <- pwy_abundance %>%
  select(SN, all_of(available_pwys)) %>%
  inner_join(meta_filtered %>% select(SN, st131_detect, carrier), by = "SN")

calculate_group_stats <- function(data, group_var) {
  group_data <- data %>%
    select(SN, all_of(available_pwys), {{ group_var }}) %>%
    rename(group = {{ group_var }}) %>%
    filter(!is.na(group))

  species_cols <- available_pwys
  min_nonzero <- sapply(group_data[species_cols], function(x) {
    x_pos <- x[x > 0 & !is.na(x)]
    if (length(x_pos) == 0) {
      return(1e-10)
    }
    return(min(x_pos) / 2)
  })

  log10_mat <- as.data.frame(group_data[species_cols])
  for (col_name in species_cols) {
    x <- log10(ifelse(group_data[[col_name]] > 0, group_data[[col_name]], min_nonzero[[col_name]]))
    log10_mat[[col_name]] <- x
  }

  prevalence_df <- group_data %>%
    group_by(group) %>%
    summarise(across(all_of(species_cols), ~ mean(.x > 0, na.rm = TRUE) * 100), .groups = "drop") %>%
    pivot_longer(-group, names_to = "pwy_id", values_to = "prevalence")

  log10_abun_df <- bind_cols(group = group_data$group, log10_mat) %>%
    group_by(group) %>%
    summarise(across(all_of(species_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    pivot_longer(-group, names_to = "pwy_id", values_to = "log10_abundance")

  stats <- prevalence_df %>%
    left_join(log10_abun_df, by = c("group", "pwy_id")) %>%
    mutate(pwy_id = factor(pwy_id, levels = available_pwys))

  return(stats)
}

# Calculate stats for ST131 detection
st131_stats <- calculate_group_stats(pwy_data, st131_detect) %>%
  mutate(
    group_type = "ST131 Detection",
    group_label = case_when(
      group == "yes" ~ "ST131+",
      group == "no" ~ "ST131-"
    )
  )

# Calculate stats for carrier status
carrier_stats <- calculate_group_stats(pwy_data, carrier) %>%
  mutate(
    group_type = "Carrier Status",
    group_label = case_when(
      group == "non_carrier" ~ "Non-carrier",
      group == "int" ~ "Intermittent",
      group == "per" ~ "Persistent"
    )
  )

# Map to display names from pwy.csv for y-axis labels
display_levels <- pwy_map %>%
  filter(PWY %in% available_pwys) %>%
  mutate(display = ifelse(is.na(pwy) | pwy == "", PWY, pwy)) %>%
  pull(display)

# Prepare data for carrier-only plot
carrier_stats_only <- carrier_stats %>%
  left_join(pwy_map, by = c("pwy_id" = "PWY")) %>%
  mutate(
    display_name = ifelse(is.na(pwy) | pwy == "", as.character(pwy_id), pwy),
    group_label = factor(group_label, levels = c("Non-carrier", "Intermittent", "Persistent")),
    display_name = factor(display_name, levels = display_levels)
  )

ab_min_carrier <- suppressWarnings(min(carrier_stats_only$log10_abundance, na.rm = TRUE))
ab_max_carrier <- suppressWarnings(max(carrier_stats_only$log10_abundance, na.rm = TRUE))
mid_val_carrier <- if (is.finite(ab_min_carrier) && is.finite(ab_max_carrier) &&
  ab_min_carrier < 0 && ab_max_carrier > 0) {
  0
} else {
  (ab_min_carrier + ab_max_carrier) / 2
}

carrier_bubble_plot <- ggplot(carrier_stats_only, aes(x = group_label, y = display_name)) +
  geom_point(
    aes(size = prevalence, fill = log10_abundance),
    alpha = 0.9,
    shape = 21,
    color = "grey30",
    stroke = 0.4
  ) +
  scale_size_continuous(
    name = "Prevalence (%)",
    range = c(4, 18),
    breaks = c(0, 25, 50, 75, 100),
    limits = c(0, 100)
  ) +
  scale_fill_gradient2(
    name = "log10(Abundance)",
    low = "#8C7BD9",
    mid = "#FFFFFF",
    high = "#F58E7E",
    midpoint = mid_val_carrier,
    limits = c(ab_min_carrier, ab_max_carrier),
    oob = scales::squish
  ) +
  scale_x_discrete(position = "top", drop = FALSE) +
  coord_cartesian(clip = "off") +
  labs(
    x = element_blank(),
    y = "Pathways"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank(),
    strip.text = element_text(size = 12, face = "bold", color = "black"),
    axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0, size = 10, margin = margin(b = 6)),
    axis.text.x.bottom = element_blank(),
    axis.ticks.x.bottom = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.margin = margin(t = 34, r = 12, b = 12, l = 12),
    legend.position = "right",
    legend.box = "vertical"
  ) +
  guides(
    size = guide_legend(order = 1, override.aes = list(shape = 21, fill = "white", color = "grey30")),
    fill = guide_colorbar(order = 2, barwidth = 0.8, barheight = 4)
  )

group_levels <- c("ST131+", "ST131-", " ", "Non-carrier", "Intermittent", "Persistent")
vline_positions <- which(group_levels != " ")

combined_stats <- bind_rows(st131_stats, carrier_stats) %>%
  left_join(pwy_map, by = c("pwy_id" = "PWY")) %>%
  mutate(display_name = ifelse(is.na(pwy) | pwy == "", as.character(pwy_id), pwy)) %>%
  mutate(
    group_type = factor(group_type, levels = c("ST131 Detection", "Carrier Status")),
    group_label = factor(group_label, levels = group_levels),
    display_name = factor(display_name, levels = display_levels)
  )

ab_min <- suppressWarnings(min(combined_stats$log10_abundance, na.rm = TRUE))
ab_max <- suppressWarnings(max(combined_stats$log10_abundance, na.rm = TRUE))
mid_val <- if (is.finite(ab_min) && is.finite(ab_max) && ab_min < 0 && ab_max > 0) 0 else (ab_min + ab_max) / 2

bubble_plot <- ggplot(combined_stats, aes(x = group_label, y = display_name)) +
  geom_vline(xintercept = vline_positions, color = "#CFCFCF", linewidth = 0.3, inherit.aes = FALSE) +
  geom_point(
    aes(size = prevalence, fill = log10_abundance),
    alpha = 0.9,
    shape = 21,
    color = "grey30",
    stroke = 0.4
  ) +
  scale_size_continuous(
    name = "Prevalence (%)",
    range = c(4, 18),
    breaks = c(0, 25, 50, 75, 100),
    limits = c(0, 100)
  ) +
  scale_fill_gradient2(
    name = "log10(Abundance)",
    low = "#8C7BD9",
    mid = "#FFFFFF",
    high = "#F58E7E",
    midpoint = mid_val,
    limits = c(ab_min, ab_max),
    oob = scales::squish
  ) +
  scale_x_discrete(position = "top", limits = group_levels, drop = FALSE) +
  coord_cartesian(clip = "off") +
  labs(
    # title = "Pathway Prevalence and Abundance by ST131 Detection and Carrier Status",
    x = element_blank(),
    y = "Pathways"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank(),
    strip.text = element_text(size = 12, face = "bold", color = "black"),
    axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0, size = 10, margin = margin(b = 6)),
    axis.text.x.bottom = element_blank(),
    axis.ticks.x.bottom = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.margin = margin(t = 34, r = 12, b = 12, l = 12),
    legend.position = "right",
    legend.box = "vertical"
  ) +
  guides(
    size = guide_legend(order = 1, override.aes = list(shape = 21, fill = "white", color = "grey30")),
    fill = guide_colorbar(order = 2, barwidth = 0.8, barheight = 4)
  )

ggsave("pwy_bubble_plot.svg", bubble_plot, width = 12, height = 8, dpi = 300)
ggsave("pwy_bubble_plot_carrier_only.svg", carrier_bubble_plot, width = 8, height = 8, dpi = 300)


# Print summary statistics
cat("Summary of prevalence and abundance by group:\n")
print(combined_stats %>%
  group_by(group_type, group_label) %>%
  summarise(
    mean_prevalence = mean(prevalence, na.rm = TRUE),
    mean_log10_abundance = mean(log10_abundance, na.rm = TRUE),
    .groups = "drop"
  ))
cat("Plots saved as pwy_bubble_plot.svg and pwy_bubble_plot_carrier_only.svg\n")
