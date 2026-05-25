library(tidyverse)
library(readxl)
library(patchwork)

pwy_names <- read_excel("meta_2025_v7.xlsx", sheet = "pwynames") %>%
  dplyr::select(PWY, Names)

create_barplot <- function(df, x_var, y_var, fill_var, title = NULL, x_lab = "", y_lab = "Log Fold Change (log10 scale)", width = NULL) {
  p <- ggplot(df, aes(x = !!sym(x_var), y = !!sym(y_var), fill = !!sym(fill_var)))

  if (is.null(width)) {
    p <- p + geom_col()
  } else {
    p <- p + geom_col(width = width)
  }

  p + scale_fill_manual(
    values = c("Enriched" = "#b03535", "Depleted" = "#2a9df4", "Not Significant" = "grey70"),
    name = "Effect", drop = FALSE
  ) +
    coord_flip() +
    labs(x = x_lab, y = y_lab, title = title) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 12),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5)
}

# Combined Pathway Carrier Plot
prepare_comparison_data <- function(file, term, label, pwy_list) {
  read.csv(file) %>%
    filter(term == !!term, PWY %in% pwy_list) %>%
    left_join(pwy_names, by = "PWY") %>%
    mutate(
      display_name = ifelse(is.na(Names), PWY, Names),
      comparison = label,
      log_fold_change = estimate,
      effect_direction = factor(ifelse(padj < 0.05, ifelse(estimate > 0, "Enriched", "Depleted"), "Not Significant"),
        levels = c("Enriched", "Depleted", "Not Significant")
      ),
      alpha_value = ifelse(padj < 0.05, 0.8, 0.3)
    ) %>%
    dplyr::select(display_name, comparison, log_fold_change, effect_direction, alpha_value)
}

carrier_res <- read.csv("LR_results/sigpwy_carrier.csv")
per_int_res <- read.csv("LR_results/sigpwy_persistent_vs_intermittent.csv")

sig_pwy <- unique(c(
  carrier_res %>% filter(term == "carrierper", abs(estimate) > 0.3, padj < 0.05) %>% pull(PWY),
  carrier_res %>% filter(term == "carrierint", abs(estimate) > 0.3, padj < 0.05) %>% pull(PWY),
  per_int_res %>% filter(term == "carrierper", abs(estimate) > 0.3, padj < 0.05) %>% pull(PWY)
))

combined_df <- bind_rows(
  prepare_comparison_data("LR_results/sigpwy_carrier.csv", "carrierper", "Persistent vs Non-carrier", sig_pwy),
  prepare_comparison_data("LR_results/sigpwy_persistent_vs_intermittent.csv", "carrierper", "Persistent vs Intermittent", sig_pwy),
  prepare_comparison_data("LR_results/sigpwy_carrier.csv", "carrierint", "Intermittent vs Non-carrier", sig_pwy)
)

pwy_order <- combined_df %>%
  group_by(display_name) %>%
  summarize(m = max(abs(log_fold_change))) %>%
  arrange(m) %>%
  pull(display_name)
combined_df$display_name <- factor(combined_df$display_name, levels = pwy_order)

build_panel <- function(df, title, show_y = TRUE) {
  p <- create_barplot(df, "display_name", "log_fold_change", "effect_direction", title = title, x_lab = "Pathway") +
    aes(alpha = alpha_value) + scale_alpha_identity()
  if (!show_y) p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank())
  p
}

p1 <- build_panel(combined_df %>% filter(comparison == "Persistent vs Non-carrier"), "Persistent vs Non-carrier") + theme(legend.position = "none")
p2 <- build_panel(combined_df %>% filter(comparison == "Persistent vs Intermittent"), "Persistent vs Intermittent", FALSE)
p3 <- build_panel(combined_df %>% filter(comparison == "Intermittent vs Non-carrier"), "Intermittent vs Non-carrier", FALSE) + theme(legend.position = "none")

ggsave("combined_pathway_carrier.svg", (p1 | p2 | p3), width = 16, height = 10)

# 2. Pathway ST131 Detection Plot
pwy_st131 <- read.csv("LR_results/sigpwy_st131.csv") %>%
  filter(term == "st131_detectyes", abs(estimate) > 0.3, padj < 0.05) %>%
  left_join(pwy_names, by = "PWY") %>%
  mutate(
    display_name = factor(ifelse(is.na(Names) | Names == "", PWY, Names)),
    effect_direction = ifelse(estimate > 0, "Enriched", "Depleted"),
    log_fold_change = estimate
  ) %>%
  arrange(abs(log_fold_change)) %>%
  mutate(display_name = fct_inorder(display_name))

ggsave("pathway_st131detect_barplot.svg", create_barplot(pwy_st131, "display_name", "log_fold_change", "effect_direction", x_lab = "Pathway", width = 0.7) + aes(alpha = 0.8) + scale_alpha_identity(), width = 10, height = 8)

# 3. Taxa ST131 Detection Plot
taxa_st131 <- read.csv("LR_results/sigtaxa_st131.csv") %>%
  filter(term == "st131_detectyes", abs(estimate) > 0.3, padj < 0.05) %>%
  mutate(
    display_name = factor(Species),
    effect_direction = ifelse(estimate > 0, "Enriched", "Depleted"),
    log_fold_change = estimate
  ) %>%
  arrange(abs(log_fold_change)) %>%
  mutate(display_name = fct_inorder(display_name))

ggsave("taxa_st131detect_barplot.svg", create_barplot(taxa_st131, "display_name", "log_fold_change", "effect_direction", x_lab = "Taxa", width = 0.7) + aes(alpha = 0.8) + scale_alpha_identity(), width = 10, height = 8)
