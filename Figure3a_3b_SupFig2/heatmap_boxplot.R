library(tidyverse)
library(readxl)
library(patchwork)
library(rstatix)
library(ggpubr)


create_boxplot_stats <- function(data, group_col, title, y_label, colors) {
  group_sym <- sym(group_col)
  num_groups <- length(unique(data[[group_col]]))
  data <- data %>% mutate(Feature = title)

  if (num_groups > 2) {
    stat_res <- data %>%
      kruskal_test(as.formula(paste("plot_value ~", group_col))) %>%
      add_significance("p")
    pval_label <- paste0("p = ", signif(stat_res$p, 3))
  } else {
    stat_res <- data %>%
      wilcox_test(as.formula(paste("plot_value ~", group_col))) %>%
      add_significance("p")
    pval_label <- paste0("p = ", signif(stat_res$p, 3))
  }

  pval_text_df <- data.frame(
    Feature = title,
    pval_label = pval_label,
    x_pos = Inf,
    y_pos = Inf
  )

  p <- ggplot(data, aes(x = !!group_sym, y = plot_value, fill = !!group_sym)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "grey20") +
    geom_jitter(width = 0.2, alpha = 0.6, size = 1.5, color = "grey30") +
    scale_fill_manual(values = colors) +
    facet_wrap(~ Feature) +
    geom_text(
      data = pval_text_df, aes(x = x_pos, y = y_pos, label = pval_label),
      inherit.aes = FALSE, hjust = 1.1, vjust = 1.5, size = 3.5, fontface = "italic"
    ) +
    labs(x = "", y = y_label, fill = tools::toTitleCase(gsub("_", " ", group_col))) +
    theme_bw() +
    theme(legend.position = "none", strip.text = element_text(size = 9, face = "bold"))

  if (stat_res$p < 0.05) {
    if (num_groups > 2) {
      dunn_res <- data %>%
        dunn_test(as.formula(paste("plot_value ~", group_col)), p.adjust.method = "BH") %>%
        add_significance("p.adj") %>%
        add_xy_position(x = group_col, dodge = 0.8)
      sig_dunn <- dunn_res %>% filter(p.adj < 0.05)
      if (nrow(sig_dunn) > 0) {
        p <- p + stat_pvalue_manual(
          sig_dunn,
          label = "p.adj.signif",
          hide.ns = TRUE,
          tip.length = 0.01,
          step.increase = 0.05
        )
      }
    } else {
      wilcox_res_pos <- stat_res %>% add_xy_position(x = group_col)
      p <- p + stat_pvalue_manual(
        wilcox_res_pos,
        label = "p.signif",
        hide.ns = TRUE,
        tip.length = 0.01
      )
    }
  }

  p
}

add_hd_exposure_counts <- function(meta_df) {
  meta_df <- meta_df %>%
    mutate(
      hd_exposure = factor(
        hd_exposure,
        levels = c("no_exposure", "low_exposure", "high_exposure"),
        labels = c("No Exposure", "Low Exposure", "High Exposure")
      )
    )
  hd_counts <- table(meta_df$hd_exposure)
  meta_df %>%
    mutate(
      hd_exposure = factor(
        hd_exposure,
        levels = names(hd_counts),
        labels = paste0(names(hd_counts), "\n(n=", as.numeric(hd_counts), ")")
      )
    )
}

build_household_exposure <- function(meta_df) {
  household_exposure <- meta_df %>%
    group_by(household) %>%
    summarise(
      hd_exposure = case_when(
        any(carrier == "per", na.rm = TRUE) ~ "high_exposure",
        any(carrier == "int", na.rm = TRUE) ~ "low_exposure",
        all(carrier == "non_carrier", na.rm = TRUE) ~ "no_exposure",
        TRUE ~ NA_character_
      ),
      .groups = "drop"
    )
  meta_df %>%
    left_join(household_exposure, by = "household") %>%
    filter(carrier %in% c("non_carrier", "per", "int"), !is.na(carrier)) %>%
    filter(st131_detect %in% c("yes", "no"), !is.na(st131_detect)) %>%
    filter(!is.na(hd_exposure))
}

pwy_names <- read_excel("metadata.xlsx", sheet = "pwynames") %>%
  select(PWY, Names)

pwy_hd_res <- read.csv("LR_results/sigpwy_hd_exposure.csv")
pwy_hd_hl_res <- read.csv("LR_results/sigpwy_hd_exposure_high_vs_low.csv")
sig_pwy_hd <- unique(c(
  pwy_hd_res %>% filter(term == "hd_exposurehigh_exposure", abs(estimate) > 0.3, padj < 0.05) %>% pull(PWY),
  pwy_hd_res %>% filter(term == "hd_exposurelow_exposure", abs(estimate) > 0.3, padj < 0.05) %>% pull(PWY),
  pwy_hd_hl_res %>% filter(term == "hd_exposurehigh_exposure", abs(estimate) > 0.3, padj < 0.05) %>% pull(PWY)
))

if (length(sig_pwy_hd) > 0) {
  meta_hd_box <- read_excel("metadata.xlsx", sheet = "R") %>%
    build_household_exposure() %>%
    add_hd_exposure_counts()

  pwy_raw_box <- read_excel("metadata.xlsx", sheet = "pwy")
  pwy_abundance_box <- pwy_raw_box %>%
    column_to_rownames("PWY") %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("SN")

  pathway_cols_box <- setdiff(colnames(pwy_abundance_box), c("SN", "UNMAPPED", "UNINTEGRATED", "UNGROUPED"))
  sample_depths_box <- pwy_abundance_box %>%
    select(SN, all_of(pathway_cols_box)) %>%
    rowwise() %>%
    mutate(total_rpk = sum(c_across(all_of(pathway_cols_box)), na.rm = TRUE)) %>%
    ungroup() %>%
    select(SN, total_rpk)

  sig_pwy_hd_avail <- intersect(sig_pwy_hd, colnames(pwy_abundance_box))

  if (length(sig_pwy_hd_avail) > 0) {
    pwy_long_box <- pwy_abundance_box %>%
      select(SN, all_of(sig_pwy_hd_avail)) %>%
      gather(PWY, Raw_RPK, -SN) %>%
      filter(SN %in% meta_hd_box$SN) %>%
      left_join(sample_depths_box, by = "SN") %>%
      mutate(
        TPM = (as.numeric(Raw_RPK) / total_rpk) * 1e6,
        plot_value = log2(TPM + 1)
      ) %>%
      left_join(meta_hd_box %>% select(SN, hd_exposure), by = "SN") %>%
      left_join(pwy_names, by = "PWY") %>%
      mutate(display_name = if_else(is.na(Names) | Names == "", PWY, Names))

    hd_colors_box <- c("#9191BA", "#accc87", "#b03535")
    names(hd_colors_box) <- levels(meta_hd_box$hd_exposure)

    plot_list_hd_pwy <- map(sig_pwy_hd_avail, function(pwy_id) {
      d <- pwy_long_box %>% filter(PWY == pwy_id)
      feat_name <- unique(d$display_name)[1]
      create_boxplot_stats(d, "hd_exposure", str_wrap(feat_name, width = 40), "Log2(TPM)", hd_colors_box)
    })
    names(plot_list_hd_pwy) <- sig_pwy_hd_avail

    combined_hd_pwy_box <- wrap_plots(plot_list_hd_pwy, ncol = min(2, length(plot_list_hd_pwy))) +
      plot_annotation(
        title = "Pathway abundance by household exposure",
        theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
      ) +
      plot_layout(guides = "collect") &
      theme(legend.position = "bottom")

    ggsave("pathway_hd_exposure_boxplots.svg", combined_hd_pwy_box, width = 12, height = 6, dpi = 300)
  }
}

taxa_hd_res <- read.csv("LR_results/sigtaxa_hd_exposure.csv")
taxa_hd_hl_res <- read.csv("LR_results/sigtaxa_hd_exposure_high_vs_low.csv")
sig_taxa_hd <- unique(c(
  taxa_hd_res %>% filter(term == "hd_exposurehigh_exposure", abs(estimate) > 0.3, padj < 0.05) %>% pull(Species),
  taxa_hd_res %>% filter(term == "hd_exposurelow_exposure", abs(estimate) > 0.3, padj < 0.05) %>% pull(Species),
  taxa_hd_hl_res %>% filter(term == "hd_exposurehigh_exposure", abs(estimate) > 0.3, padj < 0.05) %>% pull(Species)
))

if (length(sig_taxa_hd) > 0) {
  combined_taxa_hd <- bind_rows(
    taxa_hd_res %>% filter(term == "hd_exposurehigh_exposure", Species %in% sig_taxa_hd) %>%
      mutate(display_name = Species, log_fold_change = estimate),
    taxa_hd_hl_res %>% filter(term == "hd_exposurehigh_exposure", Species %in% sig_taxa_hd) %>%
      mutate(display_name = Species, log_fold_change = estimate),
    taxa_hd_res %>% filter(term == "hd_exposurelow_exposure", Species %in% sig_taxa_hd) %>%
      mutate(display_name = Species, log_fold_change = estimate)
  )
  taxa_hd_order <- combined_taxa_hd %>%
    group_by(display_name) %>%
    summarize(m = max(abs(log_fold_change))) %>%
    arrange(m) %>%
    pull(display_name)

  meta_taxa_box <- read_excel("metadata.xlsx", sheet = "R") %>%
    build_household_exposure() %>%
    select(SN, rawseqID, carrier, st131_detect, hd_exposure) %>%
    add_hd_exposure_counts()

  taxa_data_hd <- read_excel("metadata.xlsx", sheet = "taxa") %>%
    select(-Kingdom, -Phyla, -Class, -Order, -Family, -Genus) %>%
    filter(Species %in% sig_taxa_hd) %>%
    pivot_longer(cols = -Species, names_to = "rawseqID", values_to = "relab") %>%
    mutate(
      relab = relab / 100,
      plot_value = log10(relab + 1e-6)
    ) %>%
    filter(rawseqID %in% meta_taxa_box$rawseqID) %>%
    left_join(meta_taxa_box %>% select(rawseqID, hd_exposure), by = "rawseqID")

  hd_colors_taxa <- c("#9191BA", "#accc87", "#b03535")
  names(hd_colors_taxa) <- levels(meta_taxa_box$hd_exposure)

  plot_list_taxa_hd <- list()
  for (taxon in taxa_hd_order) {
    d <- taxa_data_hd %>% filter(Species == taxon)
    if (nrow(d) < 3L) next
    plot_list_taxa_hd[[taxon]] <- create_boxplot_stats(d, "hd_exposure", taxon, "Log10(Relative Abundance)", hd_colors_taxa)
  }

  if (length(plot_list_taxa_hd) > 0) {
    taxa_box_order <- rev(taxa_hd_order[taxa_hd_order %in% names(plot_list_taxa_hd)])
    combined_taxa_hd_box <- wrap_plots(plot_list_taxa_hd[taxa_box_order], ncol = 4) +
      plot_annotation(
        title = "Taxa abundance by household exposure",
        theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
      ) +
      plot_layout(guides = "collect") &
      theme(legend.position = "bottom")

    n_panels <- length(plot_list_taxa_hd)
    ggsave(
      "taxa_hd_exposure_boxplots.svg",
      combined_taxa_hd_box,
      width = 14,
      height = max(10, 2 + 2.2 * ceiling(n_panels / 4)),
      dpi = 300
    )
  }

  meta_taxa_hm <- read_excel("metadata.xlsx", sheet = "R") %>%
    build_household_exposure() %>%
    mutate(
      hd_exposure = factor(
        hd_exposure,
        levels = c("no_exposure", "low_exposure", "high_exposure"),
        labels = c("No exposure", "Low exposure", "High exposure")
      )
    )

  taxa_long_hm <- read_excel("metadata.xlsx", sheet = "taxa") %>%
    select(-Kingdom, -Phyla, -Class, -Order, -Family, -Genus) %>%
    filter(Species %in% sig_taxa_hd) %>%
    pivot_longer(cols = -Species, names_to = "rawseqID", values_to = "relab_raw") %>%
    mutate(
      relab_raw = as.numeric(relab_raw),
      relab = relab_raw / 100,
      plot_value = log10(relab + 1e-6)
    ) %>%
    inner_join(meta_taxa_hm %>% select(rawseqID, hd_exposure), by = "rawseqID")

  hm_summary <- taxa_long_hm %>%
    group_by(Species, hd_exposure) %>%
    summarize(mean_log10_relab = mean(plot_value, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = hd_exposure, values_from = mean_log10_relab)

  col_order <- c("No exposure", "Low exposure", "High exposure")
  for (mc in setdiff(col_order, colnames(hm_summary))) {
    hm_summary[[mc]] <- NA_real_
  }
  hm_summary <- hm_summary %>%
    select(Species, all_of(col_order)) %>%
    mutate(Species_ord = factor(Species, levels = taxa_hd_order)) %>%
    arrange(desc(as.integer(Species_ord))) %>%
    filter(!is.na(Species_ord)) %>%
    select(-Species_ord)

  mat <- hm_summary %>%
    column_to_rownames("Species") %>%
    as.matrix()

  min_val <- suppressWarnings(min(mat, na.rm = TRUE))
  max_val <- suppressWarnings(max(mat, na.rm = TRUE))
  if (!is.finite(min_val)) min_val <- -6
  if (!is.finite(max_val) || max_val <= min_val) max_val <- min_val + 0.01
  brks <- min_val + (0:3 / 3) * (max_val - min_val)
  col_fun <- circlize::colorRamp2(brks, c("#053061", "#2166ac", "#f46d43", "#a50026"))

  mat_t <- t(mat)
  ht_t <- ComplexHeatmap::Heatmap(
    mat_t,
    name = "Mean\nlog10(rel+ab)",
    col = col_fun,
    cluster_rows = nrow(mat_t) > 1,
    cluster_columns = ncol(mat_t) > 1,
    clustering_method_rows = "complete",
    clustering_method_columns = "complete",
    row_names_gp = grid::gpar(fontsize = 11),
    column_names_gp = grid::gpar(fontsize = 7),
    column_names_rot = 90,
    heatmap_legend_param = list(title = "Mean log10\n(relab + 1e-6)"),
    border = TRUE
  )

  w_t <- max(13, 0.14 * ncol(mat_t) + 4)
  h_t <- max(5, 0.35 * nrow(mat_t) + 3)
  pdf("taxa_hd_exposure_log10_relab_heatmap_transposed.pdf", width = w_t, height = h_t)
  ComplexHeatmap::draw(ht_t)
  dev.off()
}
