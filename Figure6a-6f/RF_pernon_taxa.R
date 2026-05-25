library(randomForest)
library(caret)
library(pROC)
library(readxl)
library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)
library(patchwork)
library(rstatix)
library(ggpubr)
setwd("/Users/zyggyzdys/Desktop/st131")
set.seed(123)


meta <- read_excel("meta_2025_v7.xlsx", sheet = "R") %>%
  dplyr::select(SN, rawseqID, sample_code, subject_code, household, reads, timepoint, carrier,	st131_detect, age,
                sex, abx_6months,	index_pt,	st131_qpcr_trpa_pabb, st131_mnth,	st131_mnth_isolate_count, st131_wgs,
                st131pos_density_wgs,	raw_counts,	old_SN) %>%
  mutate(st131_detect=factor(st131_detect, levels=c("yes", "no")),
         index_pt=factor(index_pt, levels=c("yes", "no")),
         carrier = factor(carrier, levels =c("non_carrier", "per", "int")),
         abx_6months=factor(abx_6months, levels=c("yes","no")),
         timepoint=factor(timepoint, levels=c("Baseline", "Repeat")),
         sex=factor(sex, levels=c("M", "F")),
         subject_code=factor(subject_code),
         household=factor(household)) %>%
  filter(st131_detect %in% c("yes", "no"))

carriage_data <- meta %>%
  filter(carrier %in% c("per", "non_carrier")) %>%
  select(subject_code, SN, carrier, age, sex, abx_6months, timepoint) %>%
  mutate(
    carriage_status = factor(carrier, levels = c("non_carrier", "per"), labels = c("Non_Carrier", "Persistent")),
    age_scaled = scale(age)[,1]
  ) %>%
  select(subject_code, SN, carriage_status, age_scaled, sex, abx_6months, timepoint) %>%
  na.omit() %>%
  left_join(meta %>% dplyr::select(SN, rawseqID), by = "SN")

taxa_raw <- read_excel("meta_2025_v7.xlsx", sheet = "taxa")

taxa_abundance <- taxa_raw %>%
  dplyr::select(-Kingdom, -Phyla, -Class, -Order, -Family, -Genus) %>%
  tidyr::pivot_longer(cols = -Species, names_to = "rawseqID", values_to = "relab") %>%
  mutate(relab = as.numeric(relab) / 100) %>%
  tidyr::pivot_wider(names_from = Species, values_from = relab)

taxa_data <- carriage_data %>%
  inner_join(taxa_abundance, by = "rawseqID") %>%
  select(-SN, -timepoint, -rawseqID)

clinical_data_original <- taxa_data %>%
  select(subject_code, carriage_status, age_scaled, sex, abx_6months)

taxa_data <- taxa_data %>% na.omit()

taxa_cols_in_data <- setdiff(colnames(taxa_data), c("subject_code", "carriage_status", "age_scaled", "sex", "abx_6months"))

taxa_data[taxa_cols_in_data] <- lapply(taxa_data[taxa_cols_in_data], function(x) {
  if (!is.numeric(x)) {
    x <- as.numeric(as.character(x))
  }
  return(x)
})

taxa_data <- taxa_data %>%
  select_if(~!all(is.na(.)))

taxa_cols_in_data <- setdiff(colnames(taxa_data), c("subject_code", "carriage_status", "age_scaled", "sex", "abx_6months"))

taxa_transformed <- log10(taxa_data[taxa_cols_in_data] + 1e-6)

taxa_scaled <- scale(taxa_transformed)

clinical_vars <- clinical_data_original

taxa_data_only <- data.frame(
  subject_code = taxa_data$subject_code,
  carriage_status = taxa_data$carriage_status,
  taxa_scaled
)

cat("Checking for NaN values in taxa data...\n")
nan_count <- sum(is.nan(as.matrix(taxa_data_only[,-c(1,2)])))
cat("NaN values found:", nan_count, "\n")

if(nan_count > 0) {
  cat("Replacing NaN values with 0...\n")
  taxa_data_only[,-c(1,2)][is.nan(as.matrix(taxa_data_only[,-c(1,2)]))] <- 0
}

complete_clinical_rows <- complete.cases(clinical_vars)
complete_pathway_rows <- complete.cases(taxa_data_only)
complete_both <- complete_clinical_rows & complete_pathway_rows

clinical_vars <- clinical_vars[complete_both, ]
taxa_data_only <- taxa_data_only[complete_both, ]


clinical_vars_clean <- clinical_vars %>% na.omit()

set.seed(123)
subjects <- factor(clinical_vars_clean$subject_code)
k <- 5; repeats <- 5
cv_index <- list(); cv_indexOut <- list(); counter <- 1
for (r in 1:repeats) {
  set.seed(1000 + r)
  folds <- caret::groupKFold(subjects, k = k)
  for (i in seq_along(folds)) {
    train_idx <- folds[[i]]
    cv_index[[counter]] <- train_idx
    cv_indexOut[[counter]] <- setdiff(seq_len(nrow(clinical_vars_clean)), train_idx)
    counter <- counter + 1
  }
}

cv_control <- trainControl(
  method = "cv", number = k,
  classProbs = TRUE, summaryFunction = twoClassSummary,
  savePredictions = "final",
  index = cv_index, indexOut = cv_indexOut
)

clinical_model <- train(
  carriage_status ~ age_scaled + sex + abx_6months,
  data = clinical_vars_clean,
  method = "rf", trControl = cv_control, metric = "ROC",
  tuneGrid = expand.grid(mtry = 1:3), ntree = 1000
)


complete_subjects <- rownames(clinical_vars_clean)
taxa_data_matched <- taxa_data_only[complete_subjects, ]

target_features <- 15
n_total_features <- ncol(taxa_data_matched) - 2

current_taxa <- setdiff(colnames(taxa_data_matched), c("subject_code", "carriage_status"))

elimination_steps <- c(100, 50)

for(step_size in elimination_steps) {
  if(length(current_taxa) > step_size) {
    cat("\n=== Reducing from", length(current_taxa), "to", step_size, "taxa ===\n")

    current_data <- taxa_data_matched[, c("carriage_status", current_taxa)]

    temp_rf <- randomForest(
      carriage_status ~ .,
      data = current_data,
      ntree = 200,
      importance = TRUE
    )

    importance_scores <- importance(temp_rf)[, "MeanDecreaseGini"]
    sorted_taxa <- names(sort(importance_scores, decreasing = TRUE))

    current_taxa <- sorted_taxa[1:step_size]

    cat("Kept top", length(current_taxa), "taxa based on importance\n")
    cat("Top 5 taxa:", paste(current_taxa[1:5], collapse = ", "), "\n")
  }
}

rfe_taxa <- current_taxa

while(length(rfe_taxa) > target_features) {
  cat("Current features:", length(rfe_taxa), "→ removing worst...\n")

  rfe_data <- taxa_data_matched[, c("carriage_status", rfe_taxa)]

  rfe_model <- randomForest(
    carriage_status ~ .,
    data = rfe_data,
    ntree = 200,
    importance = TRUE
  )

  importance_scores <- importance(rfe_model)[, "MeanDecreaseGini"]
  worst_taxon <- names(which.min(importance_scores))
  rfe_taxa <- rfe_taxa[rfe_taxa != worst_taxon]

  cat("Removed:", worst_taxon, "(importance:", round(min(importance_scores), 4), ")\n")
}

cat("Final", target_features, "features selected\n")

cat("Evaluating", target_features, "-feature model with CV...\n")
final_data <- taxa_data_matched[, c("carriage_status", rfe_taxa)]

cv_aucs <- numeric()
for (i in seq_along(cv_index)) {
  train_idx <- cv_index[[i]]
  test_idx <- cv_indexOut[[i]]
  rf_model <- randomForest(carriage_status ~ ., data = final_data[train_idx, ], ntree = 200)
  prob <- predict(rf_model, final_data[test_idx, ], type = "prob")[, "Persistent"]
  actual <- final_data$carriage_status[test_idx]
  if (length(unique(actual)) == 2) {
    cv_roc <- roc(actual, prob, quiet = TRUE)
    cv_aucs <- c(cv_aucs, as.numeric(cv_roc$auc))
  }
}

if(length(cv_aucs) > 0) {
  mean_auc <- mean(cv_aucs)
  sd_auc <- sd(cv_aucs)
  cat("Final", target_features, "-feature model: CV AUC =", round(mean_auc, 3), "±", round(sd_auc, 3), "\n")
}

current_taxa <- rfe_taxa
optimal_n_features <- target_features

rfe_result <- list(
  optsize = target_features,
  optVariables = current_taxa,
  results = data.frame(
    Variables = target_features,
    ROC = ifelse(length(cv_aucs) > 0, mean_auc, NA),
    ROCSD = ifelse(length(cv_aucs) > 0, sd_auc, NA)
  )
)

selected_features <- rfe_result$optVariables
taxa_data_rfe <- taxa_data_matched %>%
  select(carriage_status, all_of(selected_features))


n_selected_features <- length(selected_features)
max_mtry_rfe <- min(floor(sqrt(n_selected_features)), 20)
mtry_values <- unique(c(1, floor(sqrt(n_selected_features)), max_mtry_rfe))

taxa_model <- train(
  carriage_status ~ .,
  data = taxa_data_rfe,
  method = "rf",
  trControl = cv_control,
  metric = "ROC",
  tuneGrid = expand.grid(mtry = mtry_values),
  ntree = 1000
)


clinical_preds <- clinical_model$pred[clinical_model$pred$mtry == clinical_model$bestTune$mtry, ]
taxa_preds <- taxa_model$pred[taxa_model$pred$mtry == taxa_model$bestTune$mtry, ]

roc_clinical <- roc(clinical_preds$obs, clinical_preds$Persistent,
                  levels = c("Non_Carrier", "Persistent"),
                  direction = "<", quiet = TRUE)
roc_taxa <- roc(taxa_preds$obs, taxa_preds$Persistent,
                 levels = c("Non_Carrier", "Persistent"),
                 direction = "<", quiet = TRUE)

delong_test <- roc.test(roc_taxa, roc_clinical, method = "delong")
p_value_taxa <- delong_test$p.value


clinical_auc <- round(auc(roc_clinical), 3)
taxa_auc <- round(auc(roc_taxa), 3)

clinical_ci <- ci.auc(roc_clinical)
taxa_ci <- ci.auc(roc_taxa)

clinical_coords <- coords(roc_clinical, "all", ret = c("threshold", "specificity", "sensitivity"))
taxa_coords <- coords(roc_taxa, "all", ret = c("threshold", "specificity", "sensitivity"))

roc_data <- rbind(
  data.frame(
    Specificity = clinical_coords$specificity,
    Sensitivity = clinical_coords$sensitivity,
    Threshold = clinical_coords$threshold,
    Model = "Clinical"
  ) %>%
    arrange(desc(Threshold)) %>%
    select(-Threshold),
  data.frame(
    Specificity = taxa_coords$specificity,
    Sensitivity = taxa_coords$sensitivity,
    Threshold = taxa_coords$threshold,
    Model = "Taxa"
  ) %>%
    arrange(desc(Threshold)) %>%
    select(-Threshold)
)

color_mapping <- c("Clinical" = "#E64B35FF", "Taxa" = "#4DBBD5FF")

clinical_auc_text <- sprintf("AUC = %.3f (95%% CI: %.3f-%.3f)",
                           clinical_auc, round(clinical_ci[1], 3), round(clinical_ci[3], 3))
taxa_auc_text <- sprintf("AUC = %.3f (95%% CI: %.3f-%.3f)",
                          taxa_auc, round(taxa_ci[1], 3), round(taxa_ci[3], 3))

roc_plot <- ggplot(roc_data, aes(x = 1 - Specificity, y = Sensitivity, color = Model)) +
  geom_path(linewidth = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", alpha = 0.7) +
  scale_color_manual(values = color_mapping) +
  annotate("text", x = 0.5, y = 0.1, label = clinical_auc_text,
           color = "#E64B35FF", size = 3, hjust = 0) +
  annotate("text", x = 0.5, y = 0.05, label = taxa_auc_text,
           color = "#4DBBD5FF", size = 3, hjust = 0) +
  labs(
    title = "Random Forest ROC Curves",
    subtitle = paste0("delong test p value: ", p_value_taxa),
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

print(roc_plot)

ggsave("random_forest.png", roc_plot, width = 8, height = 8, dpi = 300)
ggsave("random_forest.svg", roc_plot, width = 8, height = 8, dpi = 300)


taxa_importance <- varImp(taxa_model)
importance_scores <- taxa_importance$importance

cat("\nTraining full model with all taxa to get global importance ranks...\n")
full_taxa_model <- randomForest(
  carriage_status ~ .,
  data = taxa_data_matched[, c("carriage_status", setdiff(colnames(taxa_data_matched), c("subject_code", "carriage_status")))],
  ntree = 300,
  importance = TRUE
)

full_importance_raw <- importance(full_taxa_model)[, "MeanDecreaseGini"]

selected_taxa_importance <- data.frame(
  Taxa = selected_features,
  stringsAsFactors = FALSE
) %>%
  mutate(
    RFE_Model_Raw = importance_scores$Overall[match(Taxa, rownames(importance_scores))],
    Full_Model_Raw = full_importance_raw[Taxa],
    Full_Model_Rank = rank(-full_importance_raw)[Taxa],
    Full_Model_Percentile = (length(full_importance_raw) - Full_Model_Rank + 1) / length(full_importance_raw) * 100
  ) %>%
  arrange(desc(Full_Model_Raw))

raw_importance <- importance_scores$Overall

importance_df <- data.frame(
  Taxa = rownames(importance_scores),
  Raw_Gini_RFE = raw_importance,
  Raw_Gini_Full = full_importance_raw[rownames(importance_scores)],
  Importance_100 = (raw_importance / max(raw_importance)) * 100,
  Importance_100_Full = (full_importance_raw[rownames(importance_scores)] / max(full_importance_raw)) * 100,
  Importance_0_5_Full = (full_importance_raw[rownames(importance_scores)] - min(full_importance_raw)) / (max(full_importance_raw) - min(full_importance_raw)) * 5,
  Full_Rank = rank(-full_importance_raw)[rownames(importance_scores)],
  Full_Percentile = (length(full_importance_raw) - rank(-full_importance_raw)[rownames(importance_scores)] + 1) / length(full_importance_raw) * 100,
  stringsAsFactors = FALSE
) %>%
  arrange(desc(Raw_Gini_Full)) %>%
  mutate(
    display_name = Taxa,
    display_name = ifelse(nchar(display_name) > 40,
                         paste0(substr(display_name, 1, 37), "..."),
                         display_name),
    selected_by_rfe = Taxa %in% selected_features,
    Importance = Importance_100_Full
  )

cat("RFE-selected taxa features and their importance:\n")
cat("\n=== RAW IMPORTANCE SCORES (NO NORMALIZATION) ===\n")

importance_display <- importance_df %>%
  select(Taxa, display_name, Raw_Gini_Full) %>%
  rename(full_name = display_name) %>%
  arrange(desc(Raw_Gini_Full)) %>%
  rename(Raw_Importance = Raw_Gini_Full)

print(importance_display)

cat("\n=== GLOBAL IMPORTANCE RANKING (Among All Taxa) ===\n")

global_ranking <- selected_taxa_importance %>%
  mutate(full_taxa_name = Taxa) %>%
  select(Taxa, full_taxa_name, Full_Model_Rank, Full_Model_Percentile, Full_Model_Raw) %>%
  arrange(Full_Model_Rank)

print(global_ranking)

n_selected <- length(selected_features)
plot_data <- importance_display %>%
  slice_head(n = n_selected)

importance_plot <- ggplot(plot_data, aes(y = reorder(full_name, Raw_Importance), x = Raw_Importance)) +
  geom_col(fill = "#AC9ECEFF", alpha = 0.8, width = 0.7) +
  labs(
    title = "Feature Importance",
    subtitle = paste0("RFE-Selected Taxa Features (n=", n_selected, ")"),
    y = "Taxa Features",
    x = "Importance Score (Mean Decrease Gini)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 10),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)))

abundance_data <- taxa_data %>%
  select(carriage_status, all_of(selected_features)) %>%
  pivot_longer(cols = -carriage_status, names_to = "Taxa", values_to = "Relative_Abundance") %>%
  left_join(plot_data %>% select(Taxa, full_name, Raw_Importance), by = "Taxa") %>%
  filter(!is.na(full_name))

abundance_data$full_name <- factor(abundance_data$full_name,
                                   levels = levels(reorder(plot_data$full_name, plot_data$Raw_Importance)))

abundance_data$Log10_Abundance <- log10(abundance_data$Relative_Abundance + 1e-6)

dunn_res <- abundance_data %>%
  group_by(full_name) %>%
  dunn_test(Log10_Abundance ~ carriage_status, p.adjust.method = "fdr") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "full_name", dodge = 0.75)

boxplot_plot <- ggplot(abundance_data, aes(x = full_name, y = Log10_Abundance, fill = carriage_status)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7, position = position_dodge(width = 0.75)) +
  stat_pvalue_manual(dunn_res, label = "p.adj.signif", coord.flip = TRUE, hide.ns = TRUE, tip.length = 0.01) +
  coord_flip() +
  scale_fill_manual(values = c("Non_Carrier" = "#4DBBD5FF", "Persistent" = "#E64B35FF")) +
  labs(
    title = "Relative Abundance",
    subtitle = "Log10(Abundance)",
    y = "Log10 Relative Abundance",
    x = "",
    fill = "Status"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 10),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom"
  )

combined_plot <- importance_plot + boxplot_plot + plot_layout(widths = c(1, 1))

print(combined_plot)

ggsave("taxa_pernon_feature_importance_combined.png", combined_plot, width = 12, height = 8, dpi = 300)
ggsave("taxa_pernon_feature_importance_combined.svg", combined_plot, width = 12, height = 8, dpi = 300)
