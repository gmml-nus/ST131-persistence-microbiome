library(tidyverse)
library(readxl)
library(vegan)
library(broom.mixed)
library(lmerTest)

# Load and filter metadata
meta <- read_excel("metadata.xlsx", sheet = "R")

# Calculate hd_exposure at the household level
household_exposure <- meta %>%
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

meta_per <- meta %>%
  left_join(household_exposure, by = "household") %>%
  dplyr::select(SN, rawseqID, subject_code, carrier, st131_detect, hd_exposure) %>%
  mutate(
    st131_detect = factor(st131_detect, levels = c("no", "yes")),
    carrier = factor(carrier, levels = c("non_carrier", "per", "int")),
    hd_exposure = factor(hd_exposure, levels = c("no_exposure", "low_exposure", "high_exposure")),
    subject_code = factor(subject_code)
  ) %>%
  filter(carrier %in% c("non_carrier", "per", "int"), !is.na(carrier), st131_detect %in% c("yes", "no"), !is.na(hd_exposure)) %>%
  droplevels()

# Helper function to fit lmer or lm
fit_taxa_model <- function(data, response_col, fixed_part, random_group) {
  has_repeat <- any(table(data[[random_group]]) > 1)
  formula_lm <- as.formula(paste(response_col, "~", fixed_part))
  formula_lmer <- as.formula(paste(response_col, "~", fixed_part, "+ (1|", random_group, ")"))
  if (has_repeat) {
    tryCatch({ lmerTest::lmer(formula_lmer, data = data) }, error = function(e) { lm(formula_lm, data = data) })
  } else {
    lm(formula_lm, data = data)
  }
}

run_analysis <- function(sheet_name, id_col, output_prefix, meta_data) {
  df <- read_excel("metadata.xlsx", sheet = sheet_name)
  
  df_long <- df %>%
    # Remove taxonomic rank columns if they exist (except the one we want to keep)
    dplyr::select(-any_of(setdiff(c("Kingdom", "Phyla", "Class", "Order", "Family", "Genus"), id_col))) %>%
    gather(sample_id, relab, -all_of(id_col)) %>%
    mutate(relab = as.numeric(relab) / 100) %>%
    spread(all_of(id_col), relab) %>%
    filter(sample_id %in% c(meta_data$rawseqID, meta_data$SN))
  
  features <- setdiff(colnames(df_long), "sample_id")
  df_long[features] <- lapply(df_long[features], function(x) { x[is.na(x)] <- 0; as.numeric(x) })
  
  mat <- df_long %>% column_to_rownames(var = "sample_id")
  mat <- mat[rowSums(mat, na.rm = TRUE) > 0, , drop = FALSE]
  relab <- decostand(mat, MARGIN = 1, method = "total")
  
  # Prevalence filter (>5%) and log10
  f_features <- relab[, which(apply(relab, 2, function(x) length(which(x > 0)) / length(x)) > 0.05), drop = FALSE]
  ff <- apply(f_features, 2, function(x) { x[x == 0] <- min(x[x > 0]) / 2; x })
  log10_df <- as_tibble(log10(ff), rownames = "sample_id")
  
  # Join with metadata
  meta_keys <- meta_data %>%
    dplyr::select(rawseqID, SN, subject_code, carrier, st131_detect, hd_exposure) %>%
    pivot_longer(cols = c(rawseqID, SN), names_to = "key", values_to = "sample_id") %>%
    filter(!is.na(sample_id)) %>%
    distinct(sample_id, .keep_all = TRUE) %>%
    dplyr::select(-key)
  
  final_df <- log10_df %>% left_join(meta_keys, by = "sample_id") %>% drop_na(subject_code)
  final_features <- setdiff(colnames(log10_df), "sample_id")
  
  # Model 1: ST131 detection
  final_df$st131_detect <- relevel(factor(final_df$st131_detect), ref = "no")
  lm_st131 <- purrr::map_dfr(rlang::set_names(final_features), ~ broom.mixed::tidy(fit_taxa_model(final_df, .x, "st131_detect", "subject_code")), .id = id_col) %>%
    filter(effect == "fixed") %>% mutate(padj = p.adjust(p.value, method = "fdr"))
  write.csv(lm_st131, paste0("sig", output_prefix, "_st131.csv"), row.names = FALSE)
  
  # Model 2: Carrier status
  final_df$carrier <- relevel(factor(final_df$carrier), ref = "non_carrier")
  lm_carrier <- purrr::map_dfr(rlang::set_names(final_features), ~ broom.mixed::tidy(fit_taxa_model(final_df, .x, "carrier", "subject_code")), .id = id_col) %>%
    filter(effect == "fixed") %>% mutate(padj = p.adjust(p.value, method = "fdr"))
  write.csv(lm_carrier, paste0("sig", output_prefix, "_carrier.csv"), row.names = FALSE)
  
  # Model 3: Per vs Int
  per_int_df <- final_df %>% filter(carrier %in% c("per", "int")) %>% droplevels()
  per_int_df$carrier <- relevel(factor(per_int_df$carrier), ref = "int")
  lm_per_int <- purrr::map_dfr(rlang::set_names(final_features), ~ broom.mixed::tidy(fit_taxa_model(per_int_df, .x, "carrier", "subject_code")), .id = id_col) %>%
    filter(effect == "fixed") %>% mutate(padj = p.adjust(p.value, method = "fdr"))
  write.csv(lm_per_int, paste0("sig", output_prefix, "_persistent_vs_intermittent.csv"), row.names = FALSE)
  
  # Model 4: Household Exposure (High vs No, Low vs No)
  final_df$hd_exposure <- relevel(factor(final_df$hd_exposure), ref = "no_exposure")
  lm_hd_exposure <- purrr::map_dfr(rlang::set_names(final_features), ~ broom.mixed::tidy(fit_taxa_model(final_df, .x, "hd_exposure", "subject_code")), .id = id_col) %>%
    filter(effect == "fixed") %>% mutate(padj = p.adjust(p.value, method = "fdr"))
  write.csv(lm_hd_exposure, paste0("sig", output_prefix, "_hd_exposure.csv"), row.names = FALSE)
  
  # Model 5: Household Exposure (High vs Low)
  high_low_df <- final_df %>% filter(hd_exposure %in% c("high_exposure", "low_exposure")) %>% droplevels()
  high_low_df$hd_exposure <- relevel(factor(high_low_df$hd_exposure), ref = "low_exposure")
  lm_high_low <- purrr::map_dfr(rlang::set_names(final_features), ~ broom.mixed::tidy(fit_taxa_model(high_low_df, .x, "hd_exposure", "subject_code")), .id = id_col) %>%
    filter(effect == "fixed") %>% mutate(padj = p.adjust(p.value, method = "fdr"))
  write.csv(lm_high_low, paste0("sig", output_prefix, "_hd_exposure_high_vs_low.csv"), row.names = FALSE)
}

# Run the 3 feature sets (generates 3*3 = 9 CSVs)
run_analysis("taxa", "Species", "taxa", meta_per)
run_analysis("pwy", "PWY", "pwy", meta_per)
run_analysis("card", "CARD", "card", meta_per)
