<<<<<<< HEAD
# ST131 Microbiome and Carriage Analysis Pipeline

This repository contains the metadata and scripts for generating the figures in the article **"Microbiome features associated with persistent intestinal carriages of Escherichia coli ST131 in a Southeast Asian cohort study"**. It provides the data analysis pipeline and visualization scripts used in the article.

## Project Structure

The repository is organized into directories corresponding to the figures and analyses generated for the manuscript:

- **`Figure1a_1b_1c_1d/`**: Scripts for calculating and visualizing alpha and beta diversity (e.g., PCoA based on Aitchison distance) across carriage groups.
- **`Figure2a_2b_3c/`**: Pathway and taxonomic differential abundance plotting.
- **`Figure3a/` & `Figure3c_3d/`**: Taxa and pathway bubble plots depicting relative abundances and prevalence across carriage states.
- **`Figure4b/`**: Correlation analysis between specific EC numbers and *E. coli* abundance.
- **`Figure5a_5b/` & `Figure5c/`**: Analysis of Antibiotic Resistance Gene (ARG) richness and differential abundance.
- **`Figure6a-6f/`**: Machine learning pipeline using Random Forests to predict ST131 carriage status (Persistent vs. Intermittent, and Persistent vs. Non-carrier) based on clinical metadata, pathways, and taxonomic profiles. Includes Recursive Feature Elimination (RFE) and combined ROC curve plotting.
- **`LR_results/`**: Linear regression analysis scripts.

## Metadata Instructions

The master metadata file `meta_2025_v7.xlsx` is central to the entire analysis pipeline and contains several key worksheets that need to be parsed correctly by the R scripts:

- **`R` (Clinical Data):** Contains patient metadata and clinical covariates. Key columns include:
  - `SN`, `rawseqID`: Sample identifiers.
  - `carrier`: Carriage status (`non_carrier`, `per` for Persistent, `int` for Intermittent).
  - `st131_detect`: ST131 detection status (`yes` / `no`).
  - `age`, `sex`, `abx_6months`: Clinical and demographic covariates.
  - `timepoint`: Sampling time points (`Baseline` / `Repeat`).
- **`pwy` (Pathway Abundance):** Contains unnormalized reads per kilobase (RPK) for metabolic pathways. Rows are pathways (e.g. `PWY`), and columns are sample IDs (`SN`).
- **`pwynames` (Pathway Names):** A dictionary mapping `PWY` IDs to their full biological names for interpretation.
- **`taxa` (Taxonomic Abundance):** Contains relative abundance data for microbial taxa. Includes taxonomic lineage columns (Kingdom to Species) and sample columns (`rawseqID`).

## Key Dependencies

To run these scripts, you will need **R** (>= 4.0.0) and the following core packages:

* **Data Manipulation & I/O:** `dplyr`, `tidyr`, `readr`, `readxl`, `tibble`
* **Visualization:** `ggplot2`, `patchwork`, `ggpubr`
* **Statistics & ML:** `randomForest`, `caret`, `pROC`, `rstatix`, `vegan`

## Usage

1. **Clone the repository:**
   ```bash
   git clone https://github.com/your-username/st131.git
   cd st131
   ```

2. **Data Availability:**
   Ensure that the master metadata file (`meta_2025_v7.xlsx`) is present in the root directory. This file is required by all downstream analysis scripts.

3. **Running Scripts:**
   The scripts are designed to be run from the root directory of the repository. For example, to generate the Random Forest ROC curves for Persistent vs Intermittent carriage:
   ```R
   # Inside an R session from the project root
   source("Figure6a-6f/RF_perint_collect.R")
   ```
   
   This will execute the underlying clinical, pathway, and taxa random forest models, compute the 95% Confidence Intervals for AUCs, perform DeLong's tests, and output a combined `combined_ROC_perint.svg` figure.

## License

This project is licensed under the MIT License.
=======
# ST131-persistence-microbiome
Code and analysis for investigating the microbiome signatures of persistent E. coli ST131 carriage in a Southeast Asian cohort.
>>>>>>> e23bd13267c6c1520c600660dda9eae1697946e2
