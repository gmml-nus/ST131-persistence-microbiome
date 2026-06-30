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
set.seed(123)

meta <- read_excel("meta_2025_v7.xlsx", sheet = "R") %>%
  dplyr::select(
    SN, rawseqID, sample_code, subject_code, household, reads, timepoint, carrier, st131_detect, age,
    sex, abx_6months, index_pt, st131_qpcr_trpa_pabb, st131_mnth, st131_mnth_isolate_count, st131_wgs,
    st131pos_density_wgs, raw_counts, old_SN
  ) %>%
  mutate(
    st131_detect = factor(st131_detect, levels = c("yes", "no")),
    index_pt = factor(index_pt, levels = c("yes", "no")),
    carrier = factor(carrier, levels = c("non_carrier", "per", "int")),
    abx_6months = factor(abx_6months, levels = c("yes", "no")),
    timepoint = factor(timepoint, levels = c("Baseline", "Repeat")),
    sex = factor(sex, levels = c("M", "F")),
    subject_code = factor(subject_code),
    household = factor(household)
  ) %>%
  filter(st131_detect %in% c("yes", "no"))

carriage_data <- meta %>%
  filter(carrier %in% c("int", "per")) %>%
  select(subject_code, SN, carrier, age, sex, abx_6months, timepoint) %>%
  mutate(
    carriage_status = factor(carrier, levels = c("int", "per"), labels = c("Intermittent", "Persistent")),
    age_scaled = scale(age)[, 1]
  ) %>%
  select(subject_code, SN, carriage_status, age_scaled, sex, abx_6months, timepoint) %>%
  na.omit()

pwy_raw <- read_excel("meta_2025_v7.xlsx", sheet = "pwy")

pathway_abundance <- pwy_raw %>%
  column_to_rownames("PWY") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("SN")

pathway_cols <- setdiff(colnames(pathway_abundance), c("SN", "UNMAPPED", "UNINTEGRATED"))

sample_depths <- pathway_abundance %>%
  select(SN, all_of(pathway_cols)) %>%
  rowwise() %>%
  mutate(total_rpk = sum(c_across(all_of(pathway_cols)), na.rm = TRUE)) %>%
  ungroup() %>%
  select(SN, total_rpk)

pathway_abundance <- pathway_abundance %>%
  select(SN, all_of(pathway_cols))

pathway_data <- carriage_data %>%
  inner_join(pathway_abundance, by = "SN") %>%
  select(-SN, -timepoint)

clinical_data_original <- pathway_data %>%
  select(subject_code, carriage_status, age_scaled, sex, abx_6months)

pathway_data <- pathway_data %>% na.omit()

pathway_cols_in_data <- setdiff(colnames(pathway_data), c("subject_code", "carriage_status", "age_scaled", "sex", "abx_6months"))

pathway_data[pathway_cols_in_data] <- lapply(pathway_data[pathway_cols_in_data], function(x) {
  if (!is.numeric(x)) {
    x <- as.numeric(as.character(x))
  }
  return(x)
})

pathway_data <- pathway_data %>%
  select_if(~ !all(is.na(.)))

pathway_cols_in_data <- setdiff(colnames(pathway_data), c("subject_code", "carriage_status", "age_scaled", "sex", "abx_6months"))

pathway_transformed <- log10(pathway_data[pathway_cols_in_data] + 1e-6)

pathway_scaled <- scale(pathway_transformed)

clinical_vars <- clinical_data_original

pathway_data_only <- data.frame(
  subject_code = pathway_data$subject_code,
  carriage_status = pathway_data$carriage_status,
  pathway_scaled
)

cat("Checking for NaN values in pathway data...\n")
nan_count <- sum(is.nan(as.matrix(pathway_data_only[, -c(1, 2)])))
cat("NaN values found:", nan_count, "\n")

if (nan_count > 0) {
  cat("Replacing NaN values with 0...\n")
  pathway_data_only[, -c(1, 2)][is.nan(as.matrix(pathway_data_only[, -c(1, 2)]))] <- 0
}

complete_clinical_rows <- complete.cases(clinical_vars)
complete_pathway_rows <- complete.cases(pathway_data_only)
complete_both <- complete_clinical_rows & complete_pathway_rows

clinical_vars <- clinical_vars[complete_both, ]
pathway_data_only <- pathway_data_only[complete_both, ]


clinical_vars_clean <- clinical_vars %>% na.omit()

set.seed(123)
subjects <- factor(clinical_vars_clean$subject_code)
k <- 5
repeats <- 5
cv_index <- list()
cv_indexOut <- list()
counter <- 1
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
pathway_data_matched <- pathway_data_only[complete_subjects, ]

target_features <- 15
n_total_features <- ncol(pathway_data_matched) - 2

current_pathways <- setdiff(colnames(pathway_data_matched), c("subject_code", "carriage_status"))

elimination_steps <- c(100, 50)

for (step_size in elimination_steps) {
  if (length(current_pathways) > step_size) {
    cat("\n=== Reducing from", length(current_pathways), "to", step_size, "pathways ===\n")

    current_data <- pathway_data_matched[, c("carriage_status", current_pathways)]

    temp_rf <- randomForest(
      carriage_status ~ .,
      data = current_data,
      ntree = 200,
      importance = TRUE
    )

    importance_scores <- importance(temp_rf)[, "MeanDecreaseGini"]
    sorted_pathways <- names(sort(importance_scores, decreasing = TRUE))

    current_pathways <- sorted_pathways[1:step_size]

    cat("Kept top", length(current_pathways), "pathways based on importance\n")
    cat("Top 5 pathways:", paste(current_pathways[1:5], collapse = ", "), "\n")
  }
}

rfe_pathways <- current_pathways

while (length(rfe_pathways) > target_features) {
  cat("Current features:", length(rfe_pathways), "→ removing worst...\n")

  rfe_data <- pathway_data_matched[, c("carriage_status", rfe_pathways)]

  rfe_model <- randomForest(
    carriage_status ~ .,
    data = rfe_data,
    ntree = 200,
    importance = TRUE
  )

  importance_scores <- importance(rfe_model)[, "MeanDecreaseGini"]

  worst_pathway <- names(which.min(importance_scores))
  rfe_pathways <- rfe_pathways[rfe_pathways != worst_pathway]

  cat("Removed:", worst_pathway, "(importance:", round(min(importance_scores), 4), ")\n")
}

cat("Final", target_features, "features selected\n")

cat("Evaluating", target_features, "-feature model with CV...\n")
final_data <- pathway_data_matched[, c("carriage_status", rfe_pathways)]

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

if (length(cv_aucs) > 0) {
  mean_auc <- mean(cv_aucs)
  sd_auc <- sd(cv_aucs)
  cat("Final", target_features, "-feature model: CV AUC =", round(mean_auc, 3), "±", round(sd_auc, 3), "\n")
}

current_pathways <- rfe_pathways
optimal_n_features <- target_features

rfe_result <- list(
  optsize = target_features,
  optVariables = current_pathways,
  results = data.frame(
    Variables = target_features,
    ROC = ifelse(length(cv_aucs) > 0, mean_auc, NA),
    ROCSD = ifelse(length(cv_aucs) > 0, sd_auc, NA)
  )
)

selected_features <- rfe_result$optVariables
pathway_data_rfe <- pathway_data_matched %>%
  select(carriage_status, all_of(selected_features))


n_selected_features <- length(selected_features)
max_mtry_rfe <- min(floor(sqrt(n_selected_features)), 20)
mtry_values <- unique(c(1, floor(sqrt(n_selected_features)), max_mtry_rfe))

pathway_model <- train(
  carriage_status ~ .,
  data = pathway_data_rfe,
  method = "rf",
  trControl = cv_control,
  metric = "ROC",
  tuneGrid = expand.grid(mtry = mtry_values),
  ntree = 1000
)


clinical_preds <- clinical_model$pred[clinical_model$pred$mtry == clinical_model$bestTune$mtry, ]
pathway_preds <- pathway_model$pred[pathway_model$pred$mtry == pathway_model$bestTune$mtry, ]

roc_clinical <- roc(clinical_preds$obs, clinical_preds$Persistent,
  levels = c("Intermittent", "Persistent"),
  direction = "<", quiet = TRUE
)
roc_pathway <- roc(pathway_preds$obs, pathway_preds$Persistent,
  levels = c("Intermittent", "Persistent"),
  direction = "<", quiet = TRUE
)

delong_test <- roc.test(roc_pathway, roc_clinical, method = "delong")
p_value_pwy <- delong_test$p.value


clinical_auc <- round(auc(roc_clinical), 3)
pathway_auc <- round(auc(roc_pathway), 3)

clinical_ci <- ci.auc(roc_clinical)
pathway_ci <- ci.auc(roc_pathway)

clinical_coords <- coords(roc_clinical, "all", ret = c("threshold", "specificity", "sensitivity"))
pathway_coords <- coords(roc_pathway, "all", ret = c("threshold", "specificity", "sensitivity"))

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
    Specificity = pathway_coords$specificity,
    Sensitivity = pathway_coords$sensitivity,
    Threshold = pathway_coords$threshold,
    Model = "Pathway"
  ) %>%
    arrange(desc(Threshold)) %>%
    select(-Threshold)
)

color_mapping <- c("Clinical" = "#E64B35FF", "Pathway" = "#4DBBD5FF")

clinical_auc_text <- sprintf(
  "AUC = %.3f (95%% CI: %.3f-%.3f)",
  clinical_auc, round(clinical_ci[1], 3), round(clinical_ci[3], 3)
)
pathway_auc_text <- sprintf(
  "AUC = %.3f (95%% CI: %.3f-%.3f)",
  pathway_auc, round(pathway_ci[1], 3), round(pathway_ci[3], 3)
)


roc_plot <- ggplot(roc_data, aes(x = 1 - Specificity, y = Sensitivity, color = Model)) +
  geom_path(linewidth = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", alpha = 0.7) +
  scale_color_manual(values = color_mapping) +
  annotate("text",
    x = 0.5, y = 0.1, label = clinical_auc_text,
    color = "#E64B35FF", size = 3, hjust = 0
  ) +
  annotate("text",
    x = 0.5, y = 0.05, label = pathway_auc_text,
    color = "#4DBBD5FF", size = 3, hjust = 0
  ) +
  labs(
    title = "Random Forest ROC Curves: Persistent vs Intermittent",
    subtitle = paste0("All timepoints | Subject-level CV | DeLong test p = ", format.pval(p_value_pwy, digits = 3)),
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

ggsave("pathway_perint_ROC.png", roc_plot, width = 8, height = 8, dpi = 300)
ggsave("pathway_perint_ROC.svg", roc_plot, width = 8, height = 8, dpi = 300)


pwy_names <- read_excel("meta_2025_v7.xlsx", sheet = "pwynames") %>%
  dplyr::select(PWY, Names)

print(str(rfe_result$results))

if (!is.null(rfe_result$results) && nrow(rfe_result$results) > 0) {
  rfe_results_df <- rfe_result$results
  if (!"Variables" %in% colnames(rfe_results_df)) {
    names(rfe_results_df)[1] <- "Variables"
  }
} else {
  cat("Creating RFE results from our step-wise process...\n")
  rfe_results_df <- data.frame(
    Variables = target_features,
    ROC = ifelse(length(cv_aucs) > 0, mean_auc, 0.5),
    ROCSD = ifelse(length(cv_aucs) > 0, sd_auc, 0)
  )
}

pathway_importance <- varImp(pathway_model)
importance_scores <- pathway_importance$importance

cat("\nTraining full model with all pathways to get global importance ranks...\n")
full_pathway_model <- randomForest(
  carriage_status ~ .,
  data = pathway_data_matched[, c("carriage_status", setdiff(colnames(pathway_data_matched), c("subject_code", "carriage_status")))],
  ntree = 300,
  importance = TRUE
)

full_importance_raw <- importance(full_pathway_model)[, "MeanDecreaseGini"]

selected_pathway_importance <- data.frame(
  PWY = selected_features,
  stringsAsFactors = FALSE
) %>%
  mutate(
    RFE_Model_Raw = importance_scores$Overall[match(PWY, rownames(importance_scores))],
    Full_Model_Raw = full_importance_raw[PWY],
    Full_Model_Rank = rank(-full_importance_raw)[PWY],
    Full_Model_Percentile = (length(full_importance_raw) - Full_Model_Rank + 1) / length(full_importance_raw) * 100
  ) %>%
  arrange(desc(Full_Model_Raw))

raw_importance <- importance_scores$Overall

importance_df <- data.frame(
  PWY = rownames(importance_scores),
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
  left_join(pwy_names, by = "PWY") %>%
  mutate(
    display_name = ifelse(is.na(Names) | Names == "", PWY, Names),
    display_name = ifelse(nchar(display_name) > 40,
      paste0(substr(display_name, 1, 37), "..."),
      display_name
    ),
    selected_by_rfe = PWY %in% selected_features,
    Importance = Importance_100_Full
  )

cat("RFE-selected pathway features and their importance:\n")
cat("\n=== RAW IMPORTANCE SCORES (NO NORMALIZATION) ===\n")

importance_display <- importance_df %>%
  select(PWY, Raw_Gini_Full) %>%
  left_join(pwy_names, by = "PWY") %>%
  mutate(
    full_name = ifelse(is.na(Names) | Names == "", PWY, Names)
  ) %>%
  select(PWY, full_name, Raw_Gini_Full) %>%
  arrange(desc(Raw_Gini_Full)) %>%
  rename(Raw_Importance = Raw_Gini_Full)

print(importance_display)

cat("\n=== GLOBAL IMPORTANCE RANKING (Among All Pathways) ===\n")
global_ranking <- selected_pathway_importance %>%
  left_join(pwy_names, by = "PWY") %>%
  mutate(
    full_pathway_name = ifelse(is.na(Names) | Names == "", PWY, Names)
  ) %>%
  select(PWY, full_pathway_name, Full_Model_Rank, Full_Model_Percentile, Full_Model_Raw) %>%
  arrange(Full_Model_Rank)

print(global_ranking)

n_selected <- length(selected_features)
plot_data <- importance_display %>%
  slice_head(n = n_selected)

importance_plot <- ggplot(plot_data, aes(y = reorder(full_name, Raw_Importance), x = Raw_Importance)) +
  geom_col(fill = "#AC9ECEFF", alpha = 0.8, width = 0.7) +
  labs(
    title = "Feature Importance",
    subtitle = paste0("RFE-Selected Pathway Features (n=", n_selected, ")"),
    y = "Pathway Features",
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

abundance_data <- carriage_data %>%
  inner_join(pathway_abundance %>% select(SN, all_of(selected_features)), by = "SN") %>%
  inner_join(sample_depths, by = "SN") %>%
  select(carriage_status, total_rpk, all_of(selected_features)) %>%
  pivot_longer(cols = c(-carriage_status, -total_rpk), names_to = "PWY", values_to = "Raw_RPK") %>%
  mutate(Relative_Abundance = (Raw_RPK / total_rpk) * 1e6) %>%
  left_join(plot_data %>% select(PWY, full_name, Raw_Importance), by = "PWY") %>%
  filter(!is.na(full_name))

abundance_data$full_name <- factor(abundance_data$full_name,
  levels = levels(reorder(plot_data$full_name, plot_data$Raw_Importance))
)

abundance_data$Log2_Abundance <- log2(abundance_data$Relative_Abundance + 1)

dunn_res <- abundance_data %>%
  group_by(full_name) %>%
  dunn_test(Log2_Abundance ~ carriage_status, p.adjust.method = "fdr") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "full_name", dodge = 0.75)

write.csv(
  dunn_res %>%
    left_join(abundance_data %>% distinct(full_name, PWY), by = "full_name") %>%
    transmute(
      PWY,
      full_name,
      Group1 = group1,
      Group2 = group2,
      p_value = p,
      q_value = p.adj,
      Significance = p.adj.signif
    ),
  "pathway_perint_boxplot_pvalues.csv",
  row.names = FALSE
)

boxplot_plot <- ggplot(abundance_data, aes(x = full_name, y = Log2_Abundance, fill = carriage_status)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7, position = position_dodge(width = 0.75)) +
  stat_pvalue_manual(dunn_res, label = "p.adj.signif", coord.flip = TRUE, hide.ns = TRUE, tip.length = 0.01) +
  coord_flip() +
  scale_fill_manual(values = c("Intermittent" = "#4DBBD5FF", "Persistent" = "#E64B35FF")) +
  labs(
    title = "Pathway Abundance",
    subtitle = "Log2(TPM)",
    y = "Log2(TPM)",
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

ggsave("pathway_perint_feature_importance_combined.png", combined_plot, width = 14, height = 8, dpi = 300)
ggsave("pathway_perint_feature_importance_combined.svg", combined_plot, width = 14, height = 8, dpi = 300)
