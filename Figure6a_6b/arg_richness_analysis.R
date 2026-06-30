library(tidyverse)
library(readxl)
library(ggplot2)
library(dplyr)
library(rstatix)

meta_per <- read_excel("metadata.xlsx", sheet = "R") %>%
  dplyr::select(
    SN, rawseqID, sample_code, subject_code, household, reads, timepoint, carrier, st131_detect, age,
    sex, abx_6months, index_pt, st131_qpcr_trpa_pabb, st131_mnth, st131_mnth_isolate_count, st131_wgs,
    st131pos_density_wgs, raw_counts, old_SN
  ) %>%
  filter(!is.na(st131_detect)) %>%
  filter(carrier %in% c("non_carrier", "per", "int"), !is.na(carrier)) %>%
  filter(st131_detect %in% c("yes", "no")) %>%
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

# Load CARD (ARG) data
card_data <- read_excel("metadata.xlsx", sheet = "card") %>%
  pivot_longer(cols = -CARD, names_to = "SN", values_to = "abundance") %>%
  filter(SN %in% c(meta_per$rawseqID, meta_per$SN)) %>%
  mutate(presence = ifelse(abundance > 0, 1, 0))

arg_richness <- card_data %>%
  group_by(SN) %>%
  summarise(ARG_richness = sum(presence, na.rm = TRUE), .groups = "drop")

arg_meta <- meta_per %>%
  left_join(arg_richness, by = "SN") %>%
  left_join(arg_richness, by = c("rawseqID" = "SN"), suffix = c("", "_alt")) %>%
  mutate(ARG_richness = coalesce(ARG_richness, ARG_richness_alt)) %>%
  dplyr::select(-ARG_richness_alt) %>%
  filter(!is.na(ARG_richness))


arg_meta_st131 <- arg_meta
arg_meta_carrier <- arg_meta %>%
  filter(timepoint == "Baseline")

# KW tests
# ST131 detection - uses both baseline and repeat samples
kruskal_st131 <- arg_meta_st131 %>%
  kruskal_test(ARG_richness ~ st131_detect) %>%
  add_significance("p")

# Carriage status - uses only baseline samples
kruskal_carriage <- arg_meta_carrier %>%
  kruskal_test(ARG_richness ~ carrier) %>%
  add_significance("p")

# Dunn test for carriage status (post-hoc for Kruskal-Wallis) - uses only baseline samples
dunn_carriage <- arg_meta_carrier %>%
  dunn_test(ARG_richness ~ carrier, p.adjust.method = "BH") %>%
  add_significance("p.adj")

print(kruskal_st131)
print(kruskal_carriage)
print(dunn_carriage)

# plots
# A) ARG richness by carriage status (baseline samples only)
# Calculate sample sizes for labels
carrier_counts <- arg_meta_carrier %>%
  group_by(carrier) %>%
  summarise(n = n(), .groups = "drop")

carrier_labels <- setNames(
  paste0(c("non_carrier", "intermittent", "persistent"), "(n=", carrier_counts$n, ")"),
  c("non_carrier", "int", "per")
)

carriage_plot <- ggplot(arg_meta_carrier, aes(x = reorder(carrier, -ARG_richness), y = ARG_richness, fill = carrier)) +
  geom_boxplot(show.legend = FALSE, outlier.shape = NA, alpha = 0.8, width = 0.8) +
  geom_jitter(show.legend = FALSE, width = 0.25, size = 3, shape = 21, colour = "black", alpha = 0.8) +
  labs(x = "Carriage status", y = "ARG richness") +
  scale_fill_manual(values = c("#9191BA", "#accc87", "#b03535")) + # Blue, Pink, Orange
  scale_x_discrete(labels = carrier_labels) +
  theme_minimal() +
  theme(
    axis.text = element_text(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  ) +
  annotate("text",
    x = 1, y = max(arg_meta_carrier$ARG_richness),
    label = paste("Kruskal-Wallis P =", round(kruskal_carriage$p, 4)), size = 4
  )

# Add Dunn test significance stars for carriage status
if (kruskal_carriage$p < 0.05) {
  # Get significant Dunn comparisons
  sig_dunn <- dunn_carriage %>% filter(p.adj < 0.05)

  if (nrow(sig_dunn) > 0) {
    # Calculate y positions for significance lines
    y_max <- max(arg_meta_carrier$ARG_richness, na.rm = TRUE)
    y_range <- y_max - min(arg_meta_carrier$ARG_richness, na.rm = TRUE)
    y_base <- y_max + (y_range * 0.1)
    y_spacing <- y_range * 0.15

    # Add significance lines and stars
    for (i in 1:nrow(sig_dunn)) {
      comparison <- sig_dunn[i, ]
      p_adj <- comparison$p.adj
      y_pos <- y_base + (i - 1) * y_spacing

      # Determine significance stars
      if (p_adj < 0.001) {
        stars <- "***"
      } else if (p_adj < 0.01) {
        stars <- "**"
      } else if (p_adj < 0.05) {
        stars <- "*"
      } else {
        next
      }

      # Get x positions for the comparison
      group1 <- comparison$group1
      group2 <- comparison$group2

      # Map group names to x positions (based on reorder: highest richness = 1, lowest = 3)
      # The actual order from the plot is: per (highest richness), non_carrier, int (lowest richness)
      x_positions <- c("per" = 1, "non_carrier" = 2, "int" = 3)
      x1 <- x_positions[group1]
      x2 <- x_positions[group2]

      if (!is.na(x1) && !is.na(x2)) {
        # Add significance line and stars
        carriage_plot <- carriage_plot +
          annotate("segment",
            x = x1, xend = x2, y = y_pos, yend = y_pos,
            color = "black", linewidth = 0.5
          ) +
          annotate("text",
            x = (x1 + x2) / 2, y = y_pos + (y_range * 0.05),
            label = stars, size = 4, fontface = "bold"
          )
      }
    }
  }
}

print(carriage_plot)

# B) ARG richness by ST131 detection (both baseline and repeat samples)
# Calculate sample sizes for labels
st131_counts <- arg_meta_st131 %>%
  group_by(st131_detect) %>%
  summarise(n = n(), .groups = "drop")

st131_labels <- setNames(
  paste0(c("no", "yes"), "(n=", st131_counts$n, ")"),
  c("no", "yes")
)

st131_plot <- ggplot(arg_meta_st131, aes(x = reorder(st131_detect, -ARG_richness), y = ARG_richness, fill = st131_detect)) +
  geom_boxplot(show.legend = FALSE, outlier.shape = NA, alpha = 0.8, width = 0.8) +
  geom_jitter(show.legend = FALSE, width = 0.25, size = 3, shape = 21, colour = "black", alpha = 0.8) +
  labs(x = "ST131 detection status", y = "ARG richness") +
  scale_fill_manual(values = c("#accc87", "#b03535")) + # Green, Red
  scale_x_discrete(labels = st131_labels) +
  theme_minimal() +
  theme(
    axis.text = element_text(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  ) +
  annotate("text",
    x = 0.8, y = max(arg_meta_st131$ARG_richness),
    label = paste("Kruskal-Wallis P =", round(kruskal_st131$p, 4)), size = 4
  )

# Add Kruskal-Wallis significance stars for ST131 detection
if (kruskal_st131$p < 0.05) {
  # Calculate y position for significance line
  y_max <- max(arg_meta_st131$ARG_richness, na.rm = TRUE)
  y_range <- y_max - min(arg_meta_st131$ARG_richness, na.rm = TRUE)
  y_line <- y_max + (y_range * 0.1)
  y_stars <- y_line + (y_range * 0.05)

  # Determine significance stars
  if (kruskal_st131$p < 0.001) {
    stars <- "***"
  } else if (kruskal_st131$p < 0.01) {
    stars <- "**"
  } else if (kruskal_st131$p < 0.05) {
    stars <- "*"
  }

  # Add significance line and stars
  st131_plot <- st131_plot +
    annotate("segment",
      x = 1, xend = 2, y = y_line, yend = y_line,
      color = "black", linewidth = 0.5
    ) +
    annotate("text",
      x = 1.5, y = y_stars,
      label = stars, size = 4, fontface = "bold"
    )
}

print(st131_plot)

ggsave("arg_richness_carriage_status.svg", carriage_plot, width = 10, height = 8, dpi = 300)
ggsave("arg_richness_st131_detection.svg", st131_plot, width = 10, height = 8, dpi = 300)

# Add significance annotation
# if(kruskal_st131$p < 0.05) {
#   st131_plot <- st131_plot +
#     geom_segment(aes(x = 1, xend = 2, y = max(arg_meta$ARG_richness) * 0.85, yend = max(arg_meta$ARG_richness) * 0.85),
#                  inherit.aes = FALSE, linewidth = 0.5, color = "black") +
#     annotate("text", x = 1.5, y = max(arg_meta$ARG_richness) * 0.9,
#              label = ifelse(kruskal_st131$p < 0.001, "***",
#                            ifelse(kruskal_st131$p < 0.01, "**", "*")),
#              size = 6, fontface = "bold", color = "black")
# } else {
#   st131_plot <- st131_plot +
#     annotate("text", x = 1.5, y = max(arg_meta$ARG_richness) * 0.85,
#              label = paste("Kruskal-Wallis P =", round(kruskal_st131$p, 4)), size = 4)
# }

# print(st131_plot)

ggsave("arg_richness_carriage_status.png", carriage_plot, width = 10, height = 8, dpi = 300)
ggsave("arg_richness_st131_detection.png", st131_plot, width = 10, height = 8, dpi = 300)

# Save both datasets
#write.csv(arg_meta_st131, "arg_richness_data_st131.csv", row.names = FALSE) # ST131 analysis (baseline + repeat)
#write.csv(arg_meta_carrier, "arg_richness_data_carrier.csv", row.names = FALSE) # Carriage analysis (baseline only)
#write.csv(kruskal_st131, "arg_richness_st131_statistics.csv", row.names = FALSE)
#write.csv(kruskal_carriage, "arg_richness_carriage_statistics.csv", row.names = FALSE)
#write.csv(dunn_carriage, "arg_richness_dunn_carriage_statistics.csv", row.names = FALSE)
