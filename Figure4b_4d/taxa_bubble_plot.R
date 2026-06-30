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

# Filter data - remove NA in st131_detect and use only baseline samples
meta_filtered <- meta %>%
  filter(carrier %in% c("non_carrier", "per", "int"), !is.na(carrier)) %>%
  filter(st131_detect %in% c("yes", "no"))

# Load taxa data
taxa_raw <- read_excel("metadata.xlsx", sheet = "taxa") %>%
  select(-Kingdom, -Phyla, -Class, -Order, -Family, -Genus)

# Get the species from taxa.csv (only the species listed there)
taxa_list <- read_csv("taxa.csv", show_col_types = FALSE) %>%
  distinct(Species) %>%
  pull(Species)

# Filter taxa data for selected species and samples
taxa_data <- taxa_raw %>%
  filter(Species %in% taxa_list) %>%
  gather(rawseqID, relab, -Species) %>%
  mutate(relab = relab / 100) %>%
  filter(rawseqID %in% meta_filtered$rawseqID) %>%
  spread(Species, relab) %>%
  left_join(meta_filtered %>% select(rawseqID, st131_detect, carrier), by = "rawseqID")

# Calculate prevalence and log10 abundance for each group
calculate_group_stats <- function(data, group_var) {
  group_data <- data %>%
    select(rawseqID, all_of(taxa_list), {{ group_var }}) %>%
    rename(group = {{ group_var }})

  # Compute per-species pseudo-counts (half of the minimum non-zero across all samples)
  species_cols <- taxa_list
  min_nonzero <- sapply(group_data[species_cols], function(x) {
    x_pos <- x[x > 0 & !is.na(x)]
    if (length(x_pos) == 0) {
      return(1e-10)
    }
    return(min(x_pos) / 2)
  })

  # Build a log10-transformed matrix avoiding log(0)
  log10_mat <- as.data.frame(group_data[species_cols])
  for (col_name in species_cols) {
    x <- log10(ifelse(group_data[[col_name]] > 0, group_data[[col_name]], min_nonzero[[col_name]]))
    log10_mat[[col_name]] <- x
  }

  # Prevalence and mean log10 abundance by group/species
  prevalence_df <- group_data %>%
    group_by(group) %>%
    summarise(across(all_of(species_cols), ~ mean(.x > 0, na.rm = TRUE) * 100), .groups = "drop") %>%
    pivot_longer(-group, names_to = "species", values_to = "prevalence")

  log10_abun_df <- bind_cols(group = group_data$group, log10_mat) %>%
    group_by(group) %>%
    summarise(across(all_of(species_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    pivot_longer(-group, names_to = "species", values_to = "log10_abundance")

  stats <- prevalence_df %>%
    left_join(log10_abun_df, by = c("group", "species")) %>%
    mutate(species = factor(species, levels = taxa_list))

  return(stats)
}

# Calculate stats for ST131 detection
st131_stats <- calculate_group_stats(taxa_data, st131_detect) %>%
  mutate(
    group_type = "ST131 Detection",
    group_label = case_when(
      group == "yes" ~ "ST131+",
      group == "no" ~ "ST131-"
    )
  )

# Calculate stats for carrier status
carrier_stats <- calculate_group_stats(taxa_data, carrier) %>%
  mutate(
    group_type = "Carrier Status",
    group_label = case_when(
      group == "non_carrier" ~ "Non-carrier",
      group == "int" ~ "Intermittent",
      group == "per" ~ "Persistent"
    )
  )

# Combine data and set group levels with a blank separator between ST131 and carriers
group_levels <- c("ST131+", "ST131-", " ", "Non-carrier", "Intermittent", "Persistent")
combined_stats <- bind_rows(st131_stats, carrier_stats) %>%
  mutate(
    group_type = factor(group_type, levels = c("ST131 Detection", "Carrier Status")),
    group_label = factor(group_label, levels = group_levels)
  )

# Determine color scale limits and midpoint from actual data
ab_min <- suppressWarnings(min(combined_stats$log10_abundance, na.rm = TRUE))
ab_max <- suppressWarnings(max(combined_stats$log10_abundance, na.rm = TRUE))
mid_val <- if (is.finite(ab_min) && is.finite(ab_max) && ab_min < 0 && ab_max > 0) 0 else (ab_min + ab_max) / 2
vline_positions <- which(group_levels != " ")

# Create the bubble plot (design matched to pwy_bubble_plot.R)
bubble_plot <- ggplot(combined_stats, aes(x = group_label, y = species)) +
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
    low = "#8C7BD9", mid = "#FFFFFF", high = "#F58E7E",
    midpoint = mid_val, limits = c(ab_min, ab_max), oob = scales::squish
  ) +
  scale_x_discrete(position = "top", limits = group_levels, drop = FALSE) +
  coord_cartesian(clip = "off") +
  labs(
    x = element_blank(),
    y = "Species"
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

# Save the plot
ggsave("taxa_bubble_plot.svg", bubble_plot, width = 12, height = 8, dpi = 300)

# Print the plot
print(bubble_plot)

# Print summary statistics
cat("Summary of prevalence and abundance by group:\n")
print(combined_stats %>%
  group_by(group_type, group_label) %>%
  summarise(
    mean_prevalence = mean(prevalence, na.rm = TRUE),
    mean_log10_abundance = mean(log10_abundance, na.rm = TRUE),
    .groups = "drop"
  ))

cat("Plot saved as taxa_bubble_plot.svg\n")
