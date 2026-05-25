library(tidyverse)
library(readxl)
library(ggplot2)
library(dplyr)
library(patchwork)

# Load CARD results
card_st131_results <- read.csv("LR_results/sigcard_st131.csv")
card_carrier_results <- read.csv("LR_results/sigcard_carrier.csv")
card_persistent_vs_intermittent_results <- read.csv("LR_results/sigcard_per_vs_int.csv")

# Load CARD reference names
card_names <- read_excel("CARD_reference_name_ID.xlsx") %>%
  dplyr::select(CARD, ID, CARD_ID) %>%
  dplyr::rename(Gene_Name = ID, CARD_ID_Value = CARD_ID) %>%
  # Remove underscores from CARD column to match LR results format
  dplyr::mutate(CARD_clean = str_replace_all(CARD, "_", ""))

# Add display names to results (one-line: "Gene_Name (CARD_ID)")
card_st131_results <- card_st131_results %>%
  # Remove underscores from CARD column for matching
  mutate(CARD_clean = str_replace_all(CARD, "_", "")) %>%
  left_join(card_names, by = c("CARD_clean" = "CARD_clean")) %>%
  mutate(display_name = ifelse(is.na(Gene_Name), CARD, paste0(Gene_Name, " (", CARD_ID_Value, ")"))) %>%
  select(-CARD_clean) # Remove temporary column

card_carrier_results <- card_carrier_results %>%
  # Remove underscores from CARD column for matching
  mutate(CARD_clean = str_replace_all(CARD, "_", "")) %>%
  left_join(card_names, by = c("CARD_clean" = "CARD_clean")) %>%
  mutate(display_name = ifelse(is.na(Gene_Name), CARD, paste0(Gene_Name, " (", CARD_ID_Value, ")"))) %>%
  select(-CARD_clean) # Remove temporary column

card_persistent_vs_intermittent_results <- card_persistent_vs_intermittent_results %>%
  # Remove underscores from CARD column for matching
  mutate(CARD_clean = str_replace_all(CARD, "_", "")) %>%
  left_join(card_names, by = c("CARD_clean" = "CARD_clean")) %>%
  mutate(display_name = ifelse(is.na(Gene_Name), CARD, paste0(Gene_Name, " (", CARD_ID_Value, ")"))) %>%
  select(-CARD_clean) # Remove temporary column

# ============================
# 1. Identify significant ARGs
# ============================

st131_sig_args <- card_st131_results %>%
  filter(term == "st131_detectyes", abs(estimate) > 0.3, padj < 0.05) %>%
  pull(display_name) %>%
  unique()

persistent_sig_args <- card_carrier_results %>%
  filter(term == "carrierper", abs(estimate) > 0.3, padj < 0.05) %>%
  pull(display_name) %>%
  unique()

intermittent_sig_args <- card_carrier_results %>%
  filter(term == "carrierint", abs(estimate) > 0.3, padj < 0.05) %>%
  pull(display_name) %>%
  unique()

persistent_vs_int_sig_args <- card_persistent_vs_intermittent_results %>%
  filter(term == "carrierper", abs(estimate) > 0.3, padj < 0.05) %>%
  pull(display_name) %>%
  unique()

all_sig_args <- unique(c(
  st131_sig_args,
  persistent_sig_args,
  intermittent_sig_args,
  persistent_vs_int_sig_args
))

if (length(all_sig_args) == 0) {
  cat("No significant ARGs found for plotting\n")
} else {
  # ============================
  # 2. Build per-comparison data
  # ============================

  # ST131 detection
  st131_data <- card_st131_results %>%
    filter(term == "st131_detectyes", display_name %in% all_sig_args) %>%
    group_by(display_name) %>%
    slice_min(padj, n = 1) %>% # most significant if duplicates
    ungroup() %>%
    mutate(
      comparison = "ST131 Detection",
      group = case_when(
        display_name %in% st131_sig_args & estimate > 0 ~ "Enriched",
        display_name %in% st131_sig_args & estimate < 0 ~ "Depleted",
        TRUE ~ "No significant"
      ),
      alpha_value = ifelse(display_name %in% st131_sig_args, 0.8, 0.3)
    ) %>%
    select(display_name, estimate, group, comparison, alpha_value)

  # Persistent vs Non-carrier
  persistent_data <- card_carrier_results %>%
    filter(term == "carrierper", display_name %in% all_sig_args) %>%
    group_by(display_name) %>%
    slice_min(padj, n = 1) %>%
    ungroup() %>%
    mutate(
      comparison = "Persistent vs Non-carrier",
      group = case_when(
        display_name %in% persistent_sig_args & estimate > 0 ~ "Enriched",
        display_name %in% persistent_sig_args & estimate < 0 ~ "Depleted",
        TRUE ~ "No significant"
      ),
      alpha_value = ifelse(display_name %in% persistent_sig_args, 0.8, 0.3)
    ) %>%
    select(display_name, estimate, group, comparison, alpha_value)

  # Intermittent vs Non-carrier
  intermittent_data <- card_carrier_results %>%
    filter(term == "carrierint", display_name %in% all_sig_args) %>%
    group_by(display_name) %>%
    slice_min(padj, n = 1) %>%
    ungroup() %>%
    mutate(
      comparison = "Intermittent vs Non-carrier",
      group = case_when(
        display_name %in% intermittent_sig_args & estimate > 0 ~ "Enriched",
        display_name %in% intermittent_sig_args & estimate < 0 ~ "Depleted",
        TRUE ~ "No significant"
      ),
      alpha_value = ifelse(display_name %in% intermittent_sig_args, 0.8, 0.3)
    ) %>%
    select(display_name, estimate, group, comparison, alpha_value)

  # Persistent vs Intermittent
  persistent_vs_int_data <- card_persistent_vs_intermittent_results %>%
    filter(term == "carrierper", display_name %in% all_sig_args) %>%
    group_by(display_name) %>%
    slice_min(padj, n = 1) %>%
    ungroup() %>%
    mutate(
      comparison = "Persistent vs Intermittent",
      group = case_when(
        display_name %in% persistent_vs_int_sig_args & estimate > 0 ~ "Enriched",
        display_name %in% persistent_vs_int_sig_args & estimate < 0 ~ "Depleted",
        TRUE ~ "No significant"
      ),
      alpha_value = ifelse(display_name %in% persistent_vs_int_sig_args, 0.8, 0.3)
    ) %>%
    select(display_name, estimate, group, comparison, alpha_value)

  # ============================
  # 3. Combine and filter ARGs
  # ============================

  combined_data <- bind_rows(
    st131_data,
    persistent_data,
    intermittent_data,
    persistent_vs_int_data
  ) %>%
    mutate(
      comparison = factor(
        comparison,
        levels = c(
          "ST131 Detection",
          "Persistent vs Non-carrier",
          "Persistent vs Intermittent",
          "Intermittent vs Non-carrier"
        )
      )
    ) %>%
    group_by(display_name) %>%
    filter(any(group != "No significant")) %>% # keep only ARGs coloured in ≥1 comparison
    ungroup()

  # Global ARG order (same across all four panels)
  arg_order <- combined_data %>%
    group_by(display_name) %>%
    summarise(max_abs = max(abs(estimate), na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(max_abs)) %>%
    pull(display_name)

  # Apply global order
  combined_data <- combined_data %>%
    mutate(display_name = factor(display_name, levels = rev(arg_order)))

  level_order <- levels(combined_data$comparison)

  # ============================
  # 4. Build each panel explicitly
  # ============================

  # Panel 1: ST131 Detection (show y-axis, no legend)
  df1 <- combined_data %>%
    filter(comparison == level_order[1])

  p1 <- ggplot(
    df1,
    aes(x = display_name, y = estimate, fill = group, alpha = alpha_value)
  ) +
    geom_col() +
    scale_fill_manual(
      values = c("Enriched" = "#b03535", "Depleted" = "#2a9df4", "No significant" = "grey70"),
      name = "Effect"
    ) +
    scale_alpha_identity() +
    coord_flip() +
    labs(
      title = level_order[1],
      x = "Antibiotic Resistance Gene",
      y = "Log Fold Change (log10 scale)",
      fill = "Effect"
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 14),
      strip.text = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )

  # Panel 2: Persistent vs Non_carrier (hide y text, keep legend)
  df2 <- combined_data %>%
    filter(comparison == level_order[2])

  p2 <- ggplot(
    df2,
    aes(x = display_name, y = estimate, fill = group, alpha = alpha_value)
  ) +
    geom_col() +
    scale_fill_manual(
      values = c("Enriched" = "#b03535", "Depleted" = "#2a9df4", "No significant" = "grey70"),
      name = "Effect"
    ) +
    scale_alpha_identity() +
    coord_flip() +
    labs(
      title = level_order[2],
      x = "Antibiotic Resistance Gene",
      y = "Log Fold Change (log10 scale)",
      fill = "Effect"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      strip.text = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )

  # Panel 3: Persistent vs Intermittent (hide y, no legend)
  df3 <- combined_data %>%
    filter(comparison == level_order[3])

  p3 <- ggplot(
    df3,
    aes(x = display_name, y = estimate, fill = group, alpha = alpha_value)
  ) +
    geom_col() +
    scale_fill_manual(
      values = c("Enriched" = "#b03535", "Depleted" = "#2a9df4", "No significant" = "grey70"),
      name = "Effect"
    ) +
    scale_alpha_identity() +
    coord_flip() +
    labs(
      title = level_order[3],
      x = "Antibiotic Resistance Gene",
      y = "Log Fold Change (log10 scale)",
      fill = "Effect"
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      strip.text = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )

  # Panel 4: Intermittent vs Non_carrier (hide y, no legend)
  df4 <- combined_data %>%
    filter(comparison == level_order[4])

  p4 <- ggplot(
    df4,
    aes(x = display_name, y = estimate, fill = group, alpha = alpha_value)
  ) +
    geom_col() +
    scale_fill_manual(
      values = c("Enriched" = "#b03535", "Depleted" = "#2a9df4", "No significant" = "grey70"),
      name = "Effect"
    ) +
    scale_alpha_identity() +
    coord_flip() +
    labs(
      title = level_order[4],
      x = "Antibiotic Resistance Gene",
      y = "Log Fold Change (log10 scale)",
      fill = "Effect"
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      strip.text = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )

  # ============================
  # 5. Combine and save
  # ============================

  combined_plot <- (p1 | p2 | p3 | p4)

  ggsave("arg_differential_abundance_combined.svg", combined_plot,
    width = 20, height = 10, dpi = 300
  )
  print(combined_plot)
}
