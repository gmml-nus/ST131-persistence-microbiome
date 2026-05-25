library(tidyverse)
library(readxl)
library(vegan)
library(ape)
library(rstatix)

cbp1 <- c("#b03535", "#accc87", "#3549b0", "#008000", "#2a9df4", "#9f42f0", "#ffb50e", "#F39B7FB2", "#9191BA")

# Load metadata
meta_full <- read_excel("meta_2025_v7.xlsx", sheet = "R") %>%
  dplyr::select(rawseqID, timepoint, carrier, st131_detect, age, sex, abx_6months, household, raw_counts) %>%
  mutate(
    st131_detect = factor(st131_detect, levels = c("yes", "no", "NA")),
    carrier = factor(carrier, levels = c("non_carrier", "per", "int")),
    timepoint = factor(timepoint, levels = c("Baseline", "Repeat")),
    sex = factor(sex, levels = c("M", "F")),
    abx_6months = factor(abx_6months, levels = c("yes", "no")),
    household = factor(household)
  ) %>%
  filter(carrier %in% c("non_carrier", "per", "int"), !is.na(carrier)) %>%
  filter(st131_detect %in% c("yes", "no")) %>%
  filter(!is.na(timepoint))

# Load taxa
taxa_raw <- read_excel("meta_2025_v7.xlsx", sheet = "taxa") %>%
  select(-Kingdom, -Phyla, -Class, -Order, -Family, -Genus) %>%
  gather(rawseqID, relab, -Species) %>%
  mutate(relab = relab / 100) %>%
  spread(Species, relab) %>%
  filter(rawseqID %in% meta_full$rawseqID)

num_cols <- setdiff(colnames(taxa_raw), "rawseqID")
taxa_raw[num_cols] <- lapply(taxa_raw[num_cols], function(x) {
  x[is.na(x)] <- 0
  x
})

get_processed_data <- function(meta_subset, taxa_data) {
  taxa_sub <- taxa_data %>% filter(rawseqID %in% meta_subset$rawseqID)
  mat <- taxa_sub %>% column_to_rownames(var = "rawseqID")
  mat[] <- lapply(mat, as.numeric)
  mat <- mat[rowSums(mat, na.rm = TRUE) > 0, , drop = FALSE]
  relab <- decostand(mat, MARGIN = 1, method = "total")
  dist_mat <- vegdist(relab, "bray")
  return(list(relab = relab, dist = dist_mat))
}

# ANALYSIS: ALL SAMPLES (Baseline + Repeat)
data_all <- get_processed_data(meta_full, taxa_raw)

# Alpha Diversity (All Samples)
invsimp_all <- vegan::diversity(data_all$relab, index = "invsimpson")
alpha_all <- tibble(rawseqID = names(invsimp_all), invsimp = invsimp_all) %>%
  left_join(meta_full, by = "rawseqID")

# Plot: alpha_diversity_st131_detection.svg
wilcox_st131 <- alpha_all %>% wilcox_test(invsimp ~ st131_detect)
st131_p_label <- ifelse(wilcox_st131$p < 0.001, "p < 0.001",
  ifelse(wilcox_st131$p < 0.01, paste0("p = ", round(wilcox_st131$p, 3)),
    paste0("p = ", round(wilcox_st131$p, 2))
  )
)

plot_alpha_st131 <- alpha_all %>%
  ggplot(aes(x = reorder(st131_detect, -invsimp), y = invsimp, fill = st131_detect)) +
  geom_boxplot(show.legend = FALSE, outlier.shape = NA, alpha = 0.8, width = 0.8) +
  geom_jitter(show.legend = FALSE, width = 0.25, size = 3, shape = 21, colour = "black", alpha = 0.8) +
  labs(x = "ST131 Detection", y = "Inverse Simpson Index") +
  ylim(0, 40) +
  scale_x_discrete(labels = c("yes" = "yes (n=36)", "no" = "no (n=219)")) +
  scale_fill_manual(values = cbp1, breaks = c("yes", "no")) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 15), axis.title = element_text(size = 18),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = "black")
  ) +
  annotate("text", x = 2.2, y = 38, label = paste0("Wilcoxon rank-sum test\n", st131_p_label), size = 6, hjust = 1, vjust = 1)

ggsave("alpha_diversity_st131_detection.svg", plot_alpha_st131, width = 8, height = 6)

# PCoA (All Samples)
pcoa_all_res <- pcoa(data_all$dist)
pcoa_all_meta <- pcoa_all_res$vectors[, 1:5] %>%
  as_tibble() %>%
  add_column(rawseqID = rownames(data_all$relab)) %>%
  left_join(meta_full, by = "rawseqID")

set.seed(123)
adonis_st131_all <- adonis2(data_all$dist ~ st131_detect, data = pcoa_all_meta, by = "one", permutations = 9999, parallel = 7)

# Plot: pcoa_st131_all_samples.svg
plot_pcoa_st131 <- pcoa_all_meta %>%
  ggplot(aes(Axis.1, Axis.2, colour = st131_detect)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_colour_manual(values = c("yes" = "#b03535", "no" = "#accc87"), name = "ST131 Detection", labels = c("yes" = "Yes", "no" = "No")) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black", size = 15), axis.title = element_text(size = 18),
    legend.position = "bottom", legend.text = element_text(size = 15), legend.title = element_text(size = 18), plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
  ) +
  xlab(paste0("PC1 (", round(pcoa_all_res$values$Relative_eig[1] * 100, 2), "%)")) +
  ylab(paste0("PC2 (", round(pcoa_all_res$values$Relative_eig[2] * 100, 2), "%)")) +
  ggtitle(sprintf("PERMANOVA: R² = %.3f, p = %.3g", adonis_st131_all$R2[1], adonis_st131_all$`Pr(>F)`[1]))

ggsave("pcoa_st131_all_samples.svg", plot_pcoa_st131, width = 8, height = 6)

# ANALYSIS: BASELINE ONLY
meta_baseline <- meta_full %>% filter(timepoint == "Baseline")
data_baseline <- get_processed_data(meta_baseline, taxa_raw)

# Alpha Diversity (Baseline)
invsimp_base <- vegan::diversity(data_baseline$relab, index = "invsimpson")
alpha_base <- tibble(rawseqID = names(invsimp_base), invsimp = invsimp_base) %>%
  left_join(meta_baseline, by = "rawseqID")

# Plot: alpha_diversity_carrier_status.svg
kurskal_carrier <- alpha_base %>% kruskal_test(invsimp ~ carrier)
carrier_p_label <- ifelse(kurskal_carrier$p < 0.001, "p < 0.001",
  ifelse(kurskal_carrier$p < 0.01, paste0("p = ", round(kurskal_carrier$p, 3)),
    paste0("p = ", round(kurskal_carrier$p, 2))
  )
)

plot_alpha_carrier <- alpha_base %>%
  ggplot(aes(x = reorder(carrier, -invsimp), y = invsimp, fill = carrier)) +
  geom_boxplot(show.legend = FALSE, outlier.shape = NA, alpha = 0.8, width = 0.8) +
  geom_jitter(show.legend = FALSE, width = 0.25, size = 3, shape = 21, colour = "black", alpha = 0.8) +
  labs(x = "Carriage Status", y = "Inverse Simpson Index") +
  ylim(0, 40) +
  scale_x_discrete(labels = c("non_carrier" = "non-carrier (n=80)", "per" = "persistant(n=9)", "int" = "intermittent (n=30)")) +
  scale_fill_manual(values = cbp1, breaks = c("per", "int", "non_carrier")) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 15), axis.title = element_text(size = 18),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = "black")
  ) +
  annotate("text", x = 3.2, y = 38, label = paste0("Kruskal-Wallis test\n", carrier_p_label), size = 6, hjust = 1, vjust = 1)

ggsave("alpha_diversity_carrier_status.svg", plot_alpha_carrier, width = 8, height = 6)

# PERMANOVA (Baseline - Multivariate)
set.seed(123)
adonis_base_multi <- adonis2(data_baseline$dist ~ st131_detect + age + sex + abx_6months + carrier + household, data = alpha_base, by = "margin", permutations = 9999)

# Plot: permanova_effect_sizes.svg
adonis_effect <- data.frame(
  factor = c("st131_detect", "age", "sex", "abx_6months", "carrier", "household"),
  R = adonis_base_multi$R2[1:6], p_adj = p.adjust(adonis_base_multi$`Pr(>F)`[1:6], method = "BH")
) %>%
  filter(!is.na(R)) %>%
  mutate(star = cut(p_adj, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), labels = c("***", "**", "*", "")))

plot_permanova_base <- adonis_effect %>%
  ggplot(aes(x = reorder(factor, R), y = R)) +
  geom_bar(stat = "identity", fill = "#9191BA", alpha = .6, width = .8) +
  geom_text(aes(label = star), hjust = -0.1, size = 4, color = "black") +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "", y = "R² (Effect Size)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave("permanova_effect_sizes.svg", plot_permanova_base, width = 8, height = 6)

# PCoA (Baseline - Carrier)
pcoa_base_res <- pcoa(data_baseline$dist)
pcoa_base_meta <- pcoa_base_res$vectors[, 1:5] %>%
  as_tibble() %>%
  add_column(rawseqID = rownames(data_baseline$relab)) %>%
  left_join(meta_baseline, by = "rawseqID")

set.seed(123)
adonis_carrier_base <- adonis2(data_baseline$dist ~ carrier, data = pcoa_base_meta, by = "one", permutations = 9999, parallel = 7)

# Plot: pcoa_carrier_status.svg
plot_pcoa_carrier <- pcoa_base_meta %>%
  ggplot(aes(Axis.1, Axis.2, colour = carrier)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_colour_manual(values = c("per" = "#b03535", "int" = "#accc87", "non_carrier" = "#3549b0"), breaks = c("per", "int", "non_carrier"), name = "Carrier Status", labels = c("per" = "Persistent", "int" = "Intermittent", "non_carrier" = "Non-carrier")) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black", size = 15), axis.title = element_text(size = 18),
    legend.position = "bottom", legend.text = element_text(size = 15), legend.title = element_text(size = 18), plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
  ) +
  xlab(paste0("PC1 (", round(pcoa_base_res$values$Relative_eig[1] * 100, 2), "%)")) +
  ylab(paste0("PC2 (", round(pcoa_base_res$values$Relative_eig[2] * 100, 2), "%)")) +
  ggtitle(sprintf("PERMANOVA: R² = %.3f, p = %.3g", adonis_carrier_base$R2[1], adonis_carrier_base$`Pr(>F)`[1]))

ggsave("pcoa_carrier_status.svg", plot_pcoa_carrier, width = 8, height = 6)
