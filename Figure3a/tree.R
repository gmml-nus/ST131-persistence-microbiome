library(taxize)
library(ape)
library(dplyr)
library(ggtree)
library(ggnewscale)
library(ggtreeExtra)
library(readxl)
library(ggsci)
library(ggplot2)

# import data
ad.meta <- read_excel("padj_carrier.xlsx", sheet = "full")
species <- ad.meta$species_v1 %>% unique()

# validate the species name in ncbi database
classifications <- classification(species, db = "ncbi")
ad.tree <- class2tree(classifications)$phylo
plot(ad.tree, cex = 0.6)

ad.tip_metadata <- read_excel("ad.tip_metadata_carrier.xlsx", sheet = "sheet1")
ad.tip_metadata <- ad.tip_metadata %>%
  rename(species_name = label)

# create the display names
create_display_names <- function(name) {
  words <- strsplit(name, " ")[[1]]
  if (length(words) > 1) {
    first_line <- words[1]
    second_line <- paste(words[2:length(words)], collapse = " ")
    return(paste(first_line, second_line, sep = "\n"))
  } else {
    return(name)
  }
}

ad.tip_metadata <- ad.tip_metadata %>%
  mutate(Effect = factor(Effect, levels = c(1, -1), labels = c("Enriched", "Depleted"))) %>%
  mutate(Phyla = factor(Phyla)) %>%
  mutate(display_name = sapply(species_name, create_display_names))

# assign a color to internal branches when all immediate children share the same phyla
assign_node_colors <- function(tree_data) {
  node_phyla <- as.character(tree_data$Phyla)
  for (iter in 1:10) {
    changed <- FALSE
    for (i in seq_len(nrow(tree_data))) {
      if (is.na(node_phyla[i])) {
        children_idx <- which(tree_data$parent == tree_data$node[i])
        if (length(children_idx) > 0) {
          child_vals <- node_phyla[children_idx]
          child_vals <- child_vals[!is.na(child_vals)]
          if (length(child_vals) > 0) {
            uniq <- unique(child_vals)
            if (length(uniq) == 1) {
              node_phyla[i] <- uniq[1]
              changed <- TRUE
            }
          }
        }
      }
    }
    if (!changed) break
  }
  return(node_phyla)
}

# build base tree and merge metadata
ptree.ad <- ggtree(ad.tree, layout = "circular")
ptree.ad <- ptree.ad %<+% ad.tip_metadata

# compute node phyla for internal branches
node_colors <- assign_node_colors(ptree.ad$data)
ptree.ad$data$Node_Phyla <- factor(node_colors, levels = levels(ad.tip_metadata$Phyla))

# internal branch coloring
ptree.ad <- ptree.ad +
  geom_tree(aes(color = Node_Phyla), size = 0.8, na.rm = FALSE) +
  geom_tippoint(aes(color = Phyla), size = 1, na.rm = FALSE) +
  scale_color_npg(na.value = "grey50") +
  theme_tree2() +
  theme(
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.position = "right",
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
ptree.ad

# show effects only for taxa that have significant results in each specific comparison
ptree.ad$data <- ptree.ad$data %>%
  mutate(across(
    c(int_non_carrier, persistent_non_carrier, persistent_int),
    ~ factor(if_else(!is.na(.x), as.character(Effect), NA_character_), levels = c("Enriched", "Depleted")),
    .names = "{.col}_combined"
  ))

final_plot <- ptree.ad +
  # add heatmap 1(int vs nc)
  new_scale_fill() +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = y, x = x, fill = int_non_carrier_combined),
    width = 5,
    color = "black",
    size = 0.5
  ) +
  scale_fill_manual(
    name = "Effect",
    values = c("Enriched" = "#b03535", "Depleted" = "#2a9df4"),
    na.value = "grey95",
    drop = TRUE
  ) +

  # add heatmap 2(per vs nc)
  new_scale_fill() +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = y, x = x, fill = persistent_non_carrier_combined),
    width = 5,
    color = "black",
    size = 0.5
  ) +
  scale_fill_manual(
    name = "Effect",
    values = c("Enriched" = "#b03535", "Depleted" = "#2a9df4"),
    na.value = "grey95",
    drop = TRUE
  ) +

  # add heatmap 3(per vs int)
  new_scale_fill() +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = y, x = x, fill = persistent_int_combined),
    width = 5,
    color = "black",
    size = 0.5
  ) +
  scale_fill_manual(
    name = "Effect",
    values = c("Enriched" = "#b03535", "Depleted" = "#2a9df4"),
    na.value = "grey95",
    drop = TRUE
  ) +

  geom_tiplab(
    aes(color = Phyla, label = display_name),
    offset = 40,
    size = 6.0,
    parse = FALSE
  )

# save the plot as SVG to modify legend and size of labels
ggsave("circular_tree_plot.svg", final_plot, width = 26, height = 22)
