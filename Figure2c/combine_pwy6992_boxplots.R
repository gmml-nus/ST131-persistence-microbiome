library(tidyverse)
library(readxl)
library(ggplot2)
library(patchwork)
library(rstatix)
library(ggpubr)

meta <- read_excel("metadata.xlsx", sheet = "R") %>%
  dplyr::select(SN, rawseqID, sample_code, subject_code, household, reads, timepoint, carrier, st131_detect, age, 
                sex, abx_6months, index_pt, st131_qpcr_trpa_pabb, st131_mnth, st131_mnth_isolate_count, st131_wgs,	
                st131pos_density_wgs, raw_counts, old_SN) %>%
  filter(st131_detect %in% c("yes", "no")) %>%
  mutate(st131_detect = factor(st131_detect, levels = c("no", "yes"), labels = c("ST131-", "ST131+")))

meta_all <- meta
meta_abx_no <- meta %>% filter(abx_6months == "no")

# load PWY data
pwy_raw <- read_excel("metadata.xlsx", sheet = "pwy")
pwy_abundance <- pwy_raw %>%
  column_to_rownames("PWY") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("SN")

# Calculate TPM
pathway_cols <- setdiff(colnames(pwy_abundance), c("SN", "UNMAPPED", "UNINTEGRATED", "UNGROUPED"))

sample_depths <- pwy_abundance %>%
  select(SN, all_of(pathway_cols)) %>%
  rowwise() %>%
  mutate(total_rpk = sum(c_across(all_of(pathway_cols)), na.rm = TRUE)) %>%
  ungroup() %>%
  select(SN, total_rpk)

# We only care about PWY_344 (PWY-6992)
pwy_data <- pwy_abundance %>%
  select(SN, PWY_344) %>%
  gather(PWY, Raw_RPK, -SN) %>%
  left_join(sample_depths, by = "SN") %>%
  mutate(
    TPM = (as.numeric(Raw_RPK) / total_rpk) * 1e6,
    plot_value = log2(TPM + 1)
  )

create_st131_boxplot <- function(meta_df, title_suffix) {
  d <- pwy_data %>%
    inner_join(meta_df %>% select(SN, st131_detect), by = "SN")
  
  st131_counts <- table(d$st131_detect)
  label_no <- paste0("ST131-\n(n=", as.numeric(st131_counts["ST131-"]), ")")
  label_yes <- paste0("ST131+\n(n=", as.numeric(st131_counts["ST131+"]), ")")
  
  d <- d %>%
    mutate(st131_label = case_when(
      st131_detect == "ST131-" ~ label_no,
      st131_detect == "ST131+" ~ label_yes
    )) %>%
    mutate(st131_label = factor(st131_label, levels = c(label_no, label_yes)))
  
  # Wilcoxon test
  stat_res <- d %>% 
    wilcox_test(plot_value ~ st131_label) %>%
    add_significance("p")
  
  pval_label <- paste0("p = ", signif(stat_res$p, 3))
  
  pval_text_df <- data.frame(
    Feature = paste0("PWY-6992\n(", title_suffix, ")"),
    pval_label = pval_label,
    x_pos = Inf,
    y_pos = Inf
  )
  
  d$Feature <- paste0("PWY-6992\n(", title_suffix, ")")
  
  p <- ggplot(d, aes(x = st131_label, y = plot_value, fill = st131_label)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "grey20") +
    geom_jitter(width = 0.2, alpha = 0.6, size = 1.5, color = "grey30") +
    scale_fill_manual(values = setNames(c("#accc87", "#b03535"), c(label_no, label_yes))) +
    facet_wrap(~ Feature) +
    geom_text(data = pval_text_df, aes(x = x_pos, y = y_pos, label = pval_label), 
              inherit.aes = FALSE, hjust = 1.1, vjust = 1.5, size = 3.5, fontface = "italic") +
    labs(
      x = "", 
      y = "Log2(TPM)",
      fill = "ST131 Detection"
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.text = element_text(size = 9, face = "bold")
    )
  
  if (stat_res$p < 0.05) {
    wilcox_res_pos <- stat_res %>% add_xy_position(x = "st131_label")
    p <- p + stat_pvalue_manual(
      wilcox_res_pos,
      label = "p.signif",
      hide.ns = TRUE,
      tip.length = 0.01
    )
  }
  
  return(p)
}

p_all <- create_st131_boxplot(meta_all, "All Samples")
p_abx_no <- create_st131_boxplot(meta_abx_no, "abx_6months == 'no'")

combined_plot <- (p_all | p_abx_no) +
  plot_annotation(title = "PWY-6992 Abundance by ST131 Detection",
                  theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))

ggsave("combined_pwy6992_st131_boxplots.png", combined_plot, width = 8, height = 5, dpi = 300)
ggsave("combined_pwy6992_st131_boxplots.svg", combined_plot, width = 8, height = 5, dpi = 300)
