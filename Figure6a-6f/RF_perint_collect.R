library(pROC)
library(ggplot2)
library(dplyr)

set.seed(123)

color_clinical <- "#E64B35FF"
color_pathway <- "#4DBBD5FF"
color_taxa <- "#00A087FF"

source("random_forest/RF_perint_pwy.R")
roc_clinical_pwy <- roc_clinical
roc_pathway_pwy <- roc_pathway

source("random_forest/RF_perint_taxa.R")
roc_taxa_taxa <- if (exists("roc_taxa")) roc_taxa else roc_pathway

clin_ci <- ci.auc(roc_clinical_pwy)
pwy_ci <- ci.auc(roc_pathway_pwy)
taxa_ci <- ci.auc(roc_taxa_taxa)

p_pwy_vs_clin <- tryCatch(pROC::roc.test(roc_pathway_pwy, roc_clinical_pwy, method = "delong")$p.value, error = function(e) NA)
p_taxa_vs_clin <- tryCatch(pROC::roc.test(roc_taxa_taxa, roc_clinical_pwy, method = "delong")$p.value, error = function(e) NA)

subtitle_text <- sprintf(
  "DeLong p: Pathway vs Clinical = %.3g; Taxa vs Clinical = %.3g",
  p_pwy_vs_clin, p_taxa_vs_clin
)

df_clinical <- coords(roc_clinical_pwy, "all", ret = c("threshold", "specificity", "sensitivity")) %>%
  rename(Specificity = specificity, Sensitivity = sensitivity, Threshold = threshold) %>%
  arrange(desc(Threshold)) %>%
  select(-Threshold) %>%
  mutate(Model = "Clinical")

df_pathway <- coords(roc_pathway_pwy, "all", ret = c("threshold", "specificity", "sensitivity")) %>%
  rename(Specificity = specificity, Sensitivity = sensitivity, Threshold = threshold) %>%
  arrange(desc(Threshold)) %>%
  select(-Threshold) %>%
  mutate(Model = "Pathway")

df_taxa <- coords(roc_taxa_taxa, "all", ret = c("threshold", "specificity", "sensitivity")) %>%
  rename(Specificity = specificity, Sensitivity = sensitivity, Threshold = threshold) %>%
  arrange(desc(Threshold)) %>%
  select(-Threshold) %>%
  mutate(Model = "Taxa")

roc_data <- bind_rows(df_clinical, df_pathway, df_taxa)

auc_texts <- list(
  Clinical = sprintf("AUC = %.3f (95%% CI: %.3f-%.3f)", as.numeric(auc(roc_clinical_pwy)), clin_ci[1], clin_ci[3]),
  Pathway  = sprintf("AUC = %.3f (95%% CI: %.3f-%.3f)", as.numeric(auc(roc_pathway_pwy)), pwy_ci[1], pwy_ci[3]),
  Taxa     = sprintf("AUC = %.3f (95%% CI: %.3f-%.3f)", as.numeric(auc(roc_taxa_taxa)), taxa_ci[1], taxa_ci[3])
)

color_mapping <- c("Clinical" = color_clinical, "Pathway" = color_pathway, "Taxa" = color_taxa)

roc_plot <- ggplot(roc_data, aes(x = 1 - Specificity, y = Sensitivity, color = Model)) +
  geom_path(linewidth = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", alpha = 0.7) +
  scale_color_manual(values = color_mapping) +
  annotate("text", x = 0.52, y = 0.12, label = auc_texts$Clinical, color = color_clinical, size = 3, hjust = 0) +
  annotate("text", x = 0.52, y = 0.08, label = auc_texts$Pathway, color = color_pathway, size = 3, hjust = 0) +
  annotate("text", x = 0.52, y = 0.04, label = auc_texts$Taxa, color = color_taxa, size = 3, hjust = 0) +
  labs(
    title = "Random Forest ROC Curves",
    subtitle = subtitle_text,
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  ) +
  guides(color = guide_legend(title = "Model")) +
  coord_equal()

ggsave("combined_ROC_perint.svg", roc_plot, width = 8, height = 8, dpi = 300)

cat("Combined ROC saved to combined_ROC_from_originals.svg\n")
