

############################################################
## 0. Clean environment
############################################################

rm(list = ls())
gc()
graphics.off()

set.seed(1234)


############################################################
## 1. Packages
############################################################

required_packages <- c(
  "Seurat",
  "Matrix",
  "dplyr",
  "tidyr",
  "ggplot2",
  "readr",
  "stringr",
  "ggpubr",
  "patchwork",
  "tibble",
  "scales",
  "DescTools"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Please install the following packages first: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(Seurat)
library(Matrix)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(stringr)
library(ggpubr)
library(patchwork)
library(tibble)
library(scales)
library(DescTools)


############################################################
## 2. Paths
############################################################

tarfile <- "/Users/yinhuang/Downloads/GSE248077_RAW.tar"

raw_dir <- "GSE248077_local"

out_dir <- "results_spatial_final"

dir.create(
  raw_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  out_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

if (!file.exists(tarfile)) {
  stop(
    "Cannot find tar file: ",
    tarfile
  )
}


############################################################
## 3. Extract GEO raw archive
############################################################

untar(
  tarfile,
  exdir = raw_dir
)

cat(
  "\nFirst files in raw directory:\n"
)

print(
  head(
    list.files(raw_dir),
    20
  )
)


############################################################
## 4. Sample information
##
## 6 Chow
## 7 CDAA-HFD
############################################################

sample_info <- data.frame(
  
  GSM = c(
    "GSM7905572",
    "GSM7905573",
    "GSM7905574",
    "GSM7905575",
    "GSM7905576",
    "GSM7905577",
    "GSM7905578",
    "GSM7905579",
    "GSM7905580",
    "GSM7905581",
    "GSM7905582",
    "GSM7905583",
    "GSM7905584"
  ),
  
  prefix = c(
    "GSM7905572_1495592",
    "GSM7905573_1495595",
    "GSM7905574_1495596",
    "GSM7905575_1513084",
    "GSM7905576_1513087",
    "GSM7905577_1513088",
    "GSM7905578_1495593",
    "GSM7905579_1495594",
    "GSM7905580_1495597",
    "GSM7905581_1513085",
    "GSM7905582_1513086",
    "GSM7905583_1513089",
    "GSM7905584_1513090"
  ),
  
  condition = c(
    rep(
      "Chow",
      6
    ),
    rep(
      "CDAA-HFD",
      7
    )
  ),
  
  stringsAsFactors = FALSE
)

print(sample_info)


############################################################
## 5. Helper functions
############################################################

check_file <- function(path) {
  
  if (!file.exists(path)) {
    
    stop(
      "Missing file: ",
      path
    )
  }
  
  invisible(TRUE)
}


############################################################
## Read Visium expression matrix
############################################################

read_visium_expr_only <- function(
    sample_id,
    prefix,
    condition,
    rawdir
) {
  
  message(
    "Reading expression: ",
    sample_id
  )
  
  matrix_file <- file.path(
    rawdir,
    paste0(
      prefix,
      "_matrix.mtx.gz"
    )
  )
  
  features_file <- file.path(
    rawdir,
    paste0(
      prefix,
      "_features.tsv.gz"
    )
  )
  
  barcodes_file <- file.path(
    rawdir,
    paste0(
      prefix,
      "_barcodes.tsv.gz"
    )
  )
  
  check_file(matrix_file)
  check_file(features_file)
  check_file(barcodes_file)
  
  counts <- ReadMtx(
    mtx = matrix_file,
    features = features_file,
    cells = barcodes_file,
    feature.column = 2,
    unique.features = TRUE
  )
  
  obj <- CreateSeuratObject(
    counts = counts,
    assay = "Spatial",
    project = sample_id,
    min.cells = 0,
    min.features = 0
  )
  
  obj$sample <- sample_id
  obj$condition <- condition
  
  DefaultAssay(obj) <- "Spatial"
  
  obj[["percent.mt"]] <- PercentageFeatureSet(
    obj,
    pattern = "^(mt-|Mt-|MT-)"
  )
  
  obj <- RenameCells(
    obj,
    add.cell.id = sample_id
  )
  
  return(obj)
}


############################################################
## Read Visium coordinates
############################################################

read_positions <- function(
    sample_id,
    prefix,
    rawdir
) {
  
  message(
    "Reading coordinates: ",
    sample_id
  )
  
  pos_file <- file.path(
    rawdir,
    paste0(
      prefix,
      "_tissue_positions_list.csv.gz"
    )
  )
  
  check_file(pos_file)
  
  pos <- read_csv(
    pos_file,
    col_names = FALSE,
    show_col_types = FALSE
  )
  
  if (ncol(pos) != 6) {
    
    stop(
      "Unexpected coordinate format for ",
      sample_id
    )
  }
  
  colnames(pos) <- c(
    "barcode",
    "in_tissue",
    "array_row",
    "array_col",
    "pxl_row",
    "pxl_col"
  )
  
  pos <- pos %>%
    mutate(
      sample = sample_id,
      
      cell = paste(
        sample_id,
        barcode,
        sep = "_"
      )
    )
  
  return(pos)
}


############################################################
## Safely intersect gene sets
############################################################

safe_intersect_genes <- function(
    genes,
    obj,
    gene_set_name
) {
  
  genes_found <- intersect(
    genes,
    rownames(obj)
  )
  
  message(
    gene_set_name,
    ": ",
    length(genes_found),
    " / ",
    length(genes),
    " genes found"
  )
  
  if (length(genes_found) == 0) {
    
    warning(
      "No genes found for gene set: ",
      gene_set_name
    )
  }
  
  return(genes_found)
}


############################################################
## P-value formatter
############################################################

format_p <- function(p) {
  
  ifelse(
    is.na(p),
    
    "NA",
    
    ifelse(
      p < 0.001,
      
      formatC(
        p,
        format = "e",
        digits = 2
      ),
      
      signif(
        p,
        3
      )
    )
  )
}


############################################################
## 6. Import all 13 samples
############################################################

obj_list <- lapply(
  seq_len(
    nrow(sample_info)
  ),
  
  function(i) {
    
    read_visium_expr_only(
      
      sample_id =
        sample_info$GSM[i],
      
      prefix =
        sample_info$prefix[i],
      
      condition =
        sample_info$condition[i],
      
      rawdir =
        raw_dir
    )
  }
)

names(obj_list) <- sample_info$GSM

cat(
  "\nImported objects:\n"
)

print(
  names(obj_list)
)


############################################################
## 7. Initial QC visualization
############################################################

qc_merged_raw <- merge(
  obj_list[[1]],
  y = obj_list[-1]
)

p_qc <- VlnPlot(
  qc_merged_raw,
  features = c(
    "nCount_Spatial",
    "nFeature_Spatial",
    "percent.mt"
  ),
  group.by = "sample",
  ncol = 3,
  pt.size = 0
)

ggsave(
  filename = file.path(
    out_dir,
    "QC_violin_all_samples.png"
  ),
  plot = p_qc,
  width = 14,
  height = 5,
  dpi = 300
)


############################################################
## 8. QC filtering
############################################################

obj_list_qc <- lapply(
  obj_list,
  
  function(obj) {
    
    subset(
      obj,
      
      subset =
        nFeature_Spatial > 200 &
        nFeature_Spatial < 7500 &
        nCount_Spatial > 500 &
        percent.mt < 25
    )
  }
)

names(obj_list_qc) <- names(obj_list)

cat(
  "\nSpot numbers after QC:\n"
)

print(
  sapply(
    obj_list_qc,
    ncol
  )
)


############################################################
## 9. Merge QC-filtered objects
############################################################

spatial_all <- merge(
  obj_list_qc[[1]],
  y = obj_list_qc[-1]
)

DefaultAssay(spatial_all) <- "Spatial"

cat(
  "\nCondition table after merge:\n"
)

print(
  table(
    spatial_all$condition
  )
)

print(spatial_all)


############################################################
## 10. Normalize
############################################################

spatial_all <- NormalizeData(
  spatial_all,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

spatial_all <- FindVariableFeatures(
  spatial_all,
  selection.method = "vst",
  nfeatures = 3000,
  verbose = FALSE
)

spatial_all <- ScaleData(
  spatial_all,
  features = VariableFeatures(
    spatial_all
  ),
  verbose = FALSE
)


############################################################
## 11. Define module gene sets
############################################################

kupffer_genes <- c(
  "Clec4f",
  "Marco",
  "Vsig4",
  "Timd4",
  "Cd5l"
)

kupffer_genes <- safe_intersect_genes(
  kupffer_genes,
  spatial_all,
  "Kupffer"
)

if (length(kupffer_genes) == 0) {
  
  stop(
    "No Kupffer markers found."
  )
}


############################################################
## 12. Add Kupffer module score
############################################################

spatial_all <- AddModuleScore(
  spatial_all,
  features = list(
    kupffer_genes
  ),
  name = "Kupffer"
)

## AddModuleScore creates Kupffer1
if (!"Kupffer1" %in% colnames(
  spatial_all@meta.data
)) {
  
  stop(
    "Kupffer1 was not created."
  )
}


############################################################
## 13. Add normalized Tbk1/Casp1 expression
############################################################

required_state_genes <- c(
  "Tbk1",
  "Casp1"
)

missing_state_genes <- setdiff(
  required_state_genes,
  rownames(spatial_all)
)

if (length(missing_state_genes) > 0) {
  
  stop(
    "Missing state genes: ",
    paste(
      missing_state_genes,
      collapse = ", "
    )
  )
}

state_expr <- FetchData(
  spatial_all,
  vars = c(
    "Tbk1",
    "Casp1"
  )
)

spatial_all$Tbk1 <- state_expr$Tbk1
spatial_all$Casp1 <- state_expr$Casp1


############################################################
## 14. Read all spatial coordinates
############################################################

coord_df <- bind_rows(
  
  lapply(
    seq_len(
      nrow(sample_info)
    ),
    
    function(i) {
      
      read_positions(
        
        sample_id =
          sample_info$GSM[i],
        
        prefix =
          sample_info$prefix[i],
        
        rawdir =
          raw_dir
      )
    }
  )
)

cat(
  "\nCoordinate samples:\n"
)

print(
  table(
    coord_df$sample
  )
)


############################################################
## 15. Join metadata and coordinates
############################################################

plot_df <- spatial_all@meta.data %>%
  rownames_to_column(
    "cell"
  ) %>%
  left_join(
    coord_df,
    by = c(
      "cell",
      "sample"
    )
  )

cat(
  "\nMissing coordinates:\n"
)

print(
  sum(
    is.na(
      plot_df$pxl_row
    )
  )
)

plot_df <- plot_df %>%
  filter(
    !is.na(pxl_row),
    !is.na(pxl_col)
  )

if (nrow(plot_df) == 0) {
  
  stop(
    "plot_df is empty after coordinate join."
  )
}


############################################################
## 16. FINAL definitions
############################################################

plot_final <- plot_df %>%
  
  group_by(
    sample
  ) %>%
  
  mutate(
    
    kup_cut = quantile(
      Kupffer1,
      0.90,
      na.rm = TRUE
    ),
    
    Kupffer_region = ifelse(
      Kupffer1 >= kup_cut,
      "Kupffer-high",
      "Other"
    ),
    
    Tbk1_positive =
      Tbk1 > 0,
    
    Casp1_positive =
      Casp1 > 0,
    
    tbk_casp_state = case_when(
      
      !Tbk1_positive &
        !Casp1_positive ~
        "Tbk1-Casp1-",
      
      !Tbk1_positive &
        Casp1_positive ~
        "Tbk1-Casp1+",
      
      Tbk1_positive &
        !Casp1_positive ~
        "Tbk1+Casp1-",
      
      Tbk1_positive &
        Casp1_positive ~
        "Tbk1+Casp1+",
      
      TRUE ~
        NA_character_
    ),
    
    double_pos = ifelse(
      tbk_casp_state ==
        "Tbk1+Casp1+",
      "Tbk1+Casp1+",
      "Other"
    ),
    
    spatial_group = case_when(
      
      double_pos == "Tbk1+Casp1+" &
        Kupffer_region == "Kupffer-high" ~
        "Double+ in Kupffer-high",
      
      double_pos == "Tbk1+Casp1+" ~
        "Double+ only",
      
      TRUE ~
        "Other"
    )
  ) %>%
  
  ungroup() %>%
  
  mutate(
    condition = factor(
      condition,
      levels = c(
        "Chow",
        "CDAA-HFD"
      )
    ),
    
    tbk_casp_state = factor(
      tbk_casp_state,
      levels = c(
        "Tbk1-Casp1-",
        "Tbk1-Casp1+",
        "Tbk1+Casp1-",
        "Tbk1+Casp1+"
      )
    )
  )


cat(
  "\nSpatial group table:\n"
)

print(
  table(
    plot_final$condition,
    plot_final$spatial_group
  )
)


############################################################
############################################################
## PLOT 1
##
## Spatial localization of
## Tbk1+Casp1 co-positive spots
############################################################
############################################################

plot_final2 <- plot_final %>%
  
  mutate(
    
    condition_label = recode(
      as.character(condition),
      
      "Chow" =
        "Control",
      
      "CDAA-HFD" =
        "CDAA-HFD"
    ),
    
    condition_label = factor(
      condition_label,
      levels = c(
        "Control",
        "CDAA-HFD"
      )
    ),
    
    spatial_group2 = case_when(
      
      spatial_group ==
        "Double+ in Kupffer-high" ~
        "Tbk1+Casp1+ in Kupffer-high",
      
      spatial_group ==
        "Double+ only" ~
        "Tbk1+Casp1+",
      
      TRUE ~
        "Other"
    )
  )


plot_final2$spatial_group2 <- factor(
  plot_final2$spatial_group2,
  
  levels = c(
    "Other",
    "Tbk1+Casp1+",
    "Tbk1+Casp1+ in Kupffer-high"
  )
)


############################################################
## Select representative sample
############################################################

representative_samples <- plot_final2 %>%
  
  count(
    condition,
    condition_label,
    sample,
    name = "n_spots"
  ) %>%
  
  group_by(
    condition
  ) %>%
  
  slice_max(
    n_spots,
    n = 1,
    with_ties = FALSE
  ) %>%
  
  ungroup()


cat(
  "\nRepresentative samples:\n"
)

print(
  representative_samples
)


spatial_show <- plot_final2 %>%
  
  semi_join(
    representative_samples %>%
      select(sample),
    
    by = "sample"
  )


p_spatial_representative <- ggplot(
  spatial_show,
  
  aes(
    x = pxl_col,
    y = -pxl_row
  )
) +
  
  geom_point(
    
    data = spatial_show %>%
      filter(
        spatial_group2 == "Other"
      ),
    
    color = "grey82",
    size = 0.28,
    alpha = 0.65
  ) +
  
  geom_point(
    
    data = spatial_show %>%
      filter(
        spatial_group2 ==
          "Tbk1+Casp1+"
      ),
    
    aes(
      color = spatial_group2
    ),
    
    size = 0.70,
    alpha = 0.90
  ) +
  
  geom_point(
    
    data = spatial_show %>%
      filter(
        spatial_group2 ==
          "Tbk1+Casp1+ in Kupffer-high"
      ),
    
    aes(
      color = spatial_group2
    ),
    
    size = 1.00,
    alpha = 1
  ) +
  
  facet_wrap(
    ~condition_label,
    nrow = 1
  ) +
  
  coord_fixed() +
  
  scale_color_manual(
    
    name = NULL,
    
    values = c(
      "Tbk1+Casp1+" =
        "#E69F00",
      
      "Tbk1+Casp1+ in Kupffer-high" =
        "#C62828"
    ),
    
    labels = c(
      expression(
        italic("Tbk1")*"+"*
          italic("Casp1")*"+ spots"
      ),
      
      expression(
        italic("Tbk1")*"+"*
          italic("Casp1")*
          "+ in Kupffer-high spots"
      )
    )
  ) +
  
  labs(
    title = expression(
      paste(
        "Spatial localization of ",
        italic("Tbk1"),
        "+",
        italic("Casp1"),
        " co-positive spots"
      )
    )
  ) +
  
  theme_void(
    base_size = 14
  ) +
  
  theme(
    
    strip.text = element_text(
      size = 13,
      face = "plain",
      color = "black"
    ),
    
    legend.position = "bottom",
    
    legend.text = element_text(
      size = 12,
      face = "plain",
      color = "black"
    ),
    
    legend.key.size = unit(
      0.45,
      "cm"
    ),
    
    plot.title = element_text(
      size = 16,
      face = "plain",
      hjust = 0.5,
      color = "black"
    ),
    
    plot.background = element_rect(
      fill = "white",
      color = NA
    ),
    
    panel.background = element_rect(
      fill = "white",
      color = NA
    ),
    
    plot.margin = margin(
      6,
      6,
      6,
      6
    )
  )


print(
  p_spatial_representative
)


ggsave(
  filename = file.path(
    out_dir,
    "Plot1_Spatial_Tbk1_Casp1_localization.png"
  ),
  plot = p_spatial_representative,
  width = 7.2,
  height = 3.4,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(
    out_dir,
    "Plot1_Spatial_Tbk1_Casp1_localization.pdf"
  ),
  plot = p_spatial_representative,
  width = 7.2,
  height = 3.4,
  bg = "white"
)


############################################################
############################################################
## PLOT 2
##
## Kupffer-associated Tbk1+Casp1 co-positivity
############################################################
############################################################

sample_quant <- plot_final %>%
  
  group_by(
    condition,
    sample
  ) %>%
  
  summarise(
    
    total_spots = n(),
    
    kupffer_high_spots = sum(
      Kupffer_region ==
        "Kupffer-high",
      na.rm = TRUE
    ),
    
    double_in_kupffer = sum(
      
      double_pos ==
        "Tbk1+Casp1+" &
        
        Kupffer_region ==
        "Kupffer-high",
      
      na.rm = TRUE
    ),
    
    frac_double_in_kupffer = ifelse(
      
      kupffer_high_spots > 0,
      
      double_in_kupffer /
        kupffer_high_spots,
      
      NA_real_
    ),
    
    .groups = "drop"
  )


cat(
  "\nSample-level quantification:\n"
)

print(
  sample_quant
)


############################################################
## Wilcoxon test
############################################################

wilcox_kupffer <- wilcox.test(
  frac_double_in_kupffer ~ condition,
  data = sample_quant,
  exact = FALSE
)

cat(
  "\nWilcoxon test:\n"
)

print(
  wilcox_kupffer
)


write.csv(
  sample_quant,
  file.path(
    out_dir,
    "Plot2_sample_level_Tbk1_Casp1_Kupffer.csv"
  ),
  row.names = FALSE
)


############################################################
## Boxplot
############################################################

p_kupffer_copositivity <- ggplot(
  
  sample_quant,
  
  aes(
    x = condition,
    y = frac_double_in_kupffer,
    fill = condition
  )
) +
  
  geom_boxplot(
    width = 0.50,
    outlier.shape = NA,
    alpha = 0.85
  ) +
  
  geom_jitter(
    width = 0.08,
    size = 2
  ) +
  
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    size = 3.5,
    label.y.npc = 0.95
  ) +
  
  scale_fill_manual(
    values = c(
      "Chow" =
        "#BDBDBD",
      
      "CDAA-HFD" =
        "#C41E1E"
    )
  ) +
  
  theme_classic(
    base_size = 11
  ) +
  
  labs(
    title =
      "Kupffer-associated Tbk1+Casp1+ co-positivity",
    
    x = NULL,
    
    y =
      "Fraction within Kupffer-high regions"
  ) +
  
  theme(
    
    legend.position =
      "none",
    
    plot.title = element_text(
      size = 12,
      face = "plain",
      hjust = 0.5
    ),
    
    axis.title.y = element_text(
      size = 10
    ),
    
    axis.text = element_text(
      size = 10
    )
  )


print(
  p_kupffer_copositivity
)


ggsave(
  filename = file.path(
    out_dir,
    "Plot2_Kupffer_associated_Tbk1_Casp1_copositivity.png"
  ),
  plot = p_kupffer_copositivity,
  width = 5.2,
  height = 3.6,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(
    out_dir,
    "Plot2_Kupffer_associated_Tbk1_Casp1_copositivity.pdf"
  ),
  plot = p_kupffer_copositivity,
  width = 5.2,
  height = 3.6,
  bg = "white"
)


############################################################
############################################################
## Prepare sensor expression
##
## Used by Plot 3, Plot 4 and Plot 5
############################################################
############################################################

sensor_genes <- c(
  Aim2 = "Aim2",
  Casp4 = "Casp4",
  Nlrp3 = "Nlrp3"
)

sensor_genes <- sensor_genes[
  sensor_genes %in%
    rownames(spatial_all)
]

cat(
  "\nSensor genes found:\n"
)

print(
  sensor_genes
)

if (length(sensor_genes) != 3) {
  
  warning(
    "Not all Aim2/Casp4/Nlrp3 genes were found."
  )
}


sensor_expr <- FetchData(
  spatial_all,
  vars = unname(
    sensor_genes
  )
)

colnames(sensor_expr) <- names(
  sensor_genes
)


############################################################
############################################################
## PLOT 3
##
## Sensor positivity in
## Tbk1+Casp1+ Kupffer-high CDAA-HFD spots
############################################################
############################################################

target_df <- plot_final %>%
  
  select(
    cell,
    sample,
    condition,
    Kupffer_region,
    Tbk1,
    Casp1
  ) %>%
  
  left_join(
    sensor_expr %>%
      rownames_to_column(
        "cell"
      ),
    
    by = "cell"
  ) %>%
  
  filter(
    condition ==
      "CDAA-HFD",
    
    Kupffer_region ==
      "Kupffer-high",
    
    Tbk1 > 0,
    
    Casp1 > 0
  )


cat(
  "\nNumber of target double-positive spots:\n"
)

print(
  nrow(target_df)
)

cat(
  "\nTarget spots per sample:\n"
)

print(
  table(
    target_df$sample
  )
)


############################################################
## Positivity matrix
##
## Positive = normalized expression > 0
############################################################

pos_mat <- target_df %>%
  
  transmute(
    
    Aim2 =
      Aim2 > 0,
    
    Casp4 =
      Casp4 > 0,
    
    Nlrp3 =
      Nlrp3 > 0
  )


############################################################
## Summary
############################################################

sensor_summary <- pos_mat %>%
  
  summarise(
    
    Aim2_positive =
      sum(Aim2),
    
    Casp4_positive =
      sum(Casp4),
    
    Nlrp3_positive =
      sum(Nlrp3),
    
    total_spots =
      n()
  ) %>%
  
  pivot_longer(
    
    cols = c(
      Aim2_positive,
      Casp4_positive,
      Nlrp3_positive
    ),
    
    names_to =
      "sensor",
    
    values_to =
      "positive_spots"
  ) %>%
  
  mutate(
    
    sensor = recode(
      
      sensor,
      
      "Aim2_positive" =
        "Aim2",
      
      "Casp4_positive" =
        "Casp4\n(Casp11)",
      
      "Nlrp3_positive" =
        "Nlrp3"
    ),
    
    positive_fraction =
      positive_spots /
      total_spots,
    
    sensor = factor(
      sensor,
      levels = c(
        "Aim2",
        "Casp4\n(Casp11)",
        "Nlrp3"
      )
    )
  )


cat(
  "\nSensor positivity summary:\n"
)

print(
  sensor_summary
)


############################################################
## Cochran's Q
############################################################

cochran_res <- CochranQTest(
  as.matrix(
    pos_mat
  )
)

cat(
  "\nCochran's Q test:\n"
)

print(
  cochran_res
)


############################################################
## Pairwise McNemar
############################################################

pairwise_mcnemar <- function(
    x,
    y,
    name_x,
    name_y
) {
  
  tab <- table(
    x,
    y
  )
  
  test <- mcnemar.test(
    tab
  )
  
  data.frame(
    comparison = paste(
      name_x,
      "vs",
      name_y
    ),
    
    p_value =
      test$p.value
  )
}


mcnemar_res <- bind_rows(
  
  pairwise_mcnemar(
    pos_mat$Aim2,
    pos_mat$Casp4,
    "Aim2",
    "Casp4"
  ),
  
  pairwise_mcnemar(
    pos_mat$Aim2,
    pos_mat$Nlrp3,
    "Aim2",
    "Nlrp3"
  ),
  
  pairwise_mcnemar(
    pos_mat$Casp4,
    pos_mat$Nlrp3,
    "Casp4",
    "Nlrp3"
  )
) %>%
  
  mutate(
    p_adj = p.adjust(
      p_value,
      method = "BH"
    )
  )


cat(
  "\nPairwise McNemar tests:\n"
)

print(
  mcnemar_res
)


############################################################
## Bar plot
############################################################

p_sensor_positivity <- ggplot(
  
  sensor_summary,
  
  aes(
    x = sensor,
    y = positive_fraction
  )
) +
  
  geom_col(
    fill = "#C41E1E",
    color = "black",
    width = 0.6,
    alpha = 0.85
  ) +
  
  geom_text(
    
    aes(
      label = paste0(
        positive_spots,
        "/",
        total_spots,
        "\n",
        round(
          positive_fraction * 100,
          1
        ),
        "%"
      )
    ),
    
    vjust = -0.3,
    size = 4
  ) +
  
  scale_y_continuous(
    
    limits = c(
      0,
      max(
        sensor_summary$positive_fraction
      ) * 1.25
    ),
    
    labels =
      percent_format(
        accuracy = 1
      )
  ) +
  
  labs(
    title =
      "Sensor positivity in Tbk1+Casp1+ Kupffer-high CDAA-HFD spots",
    
    x = NULL,
    
    y =
      "Positive fraction"
  ) +
  
  theme_classic(
    base_size = 15
  ) +
  
  theme(
    
    plot.title = element_text(
      face = "plain",
      hjust = 0.5
    ),
    
    axis.text.x = element_text(
      face = "plain"
    ),
    
    legend.position =
      "none"
  )


print(
  p_sensor_positivity
)


ggsave(
  filename = file.path(
    out_dir,
    "Plot3_Sensor_positivity_double_positive_Kupffer_high.png"
  ),
  plot = p_sensor_positivity,
  width = 6,
  height = 4.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(
    out_dir,
    "Plot3_Sensor_positivity_double_positive_Kupffer_high.pdf"
  ),
  plot = p_sensor_positivity,
  width = 6,
  height = 4.5,
  bg = "white"
)


write.csv(
  sensor_summary,
  file.path(
    out_dir,
    "Plot3_sensor_positivity_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  mcnemar_res,
  file.path(
    out_dir,
    "Plot3_pairwise_McNemar_BH.csv"
  ),
  row.names = FALSE
)


############################################################
############################################################
## PLOT 4
##
## Sensor enrichment:
##
## Tbk1+Casp1+
## versus
## Tbk1-Casp1-
##
## CDAA-HFD
## Kupffer-high only
############################################################
############################################################

analysis_df <- plot_final %>%
  
  select(
    cell,
    sample,
    condition,
    Kupffer_region,
    tbk_casp_state
  ) %>%
  
  left_join(
    sensor_expr %>%
      rownames_to_column(
        "cell"
      ),
    
    by = "cell"
  ) %>%
  
  filter(
    
    condition ==
      "CDAA-HFD",
    
    Kupffer_region ==
      "Kupffer-high",
    
    tbk_casp_state %in%
      c(
        "Tbk1-Casp1-",
        "Tbk1+Casp1+"
      )
  ) %>%
  
  mutate(
    
    tbk_casp_group = factor(
      as.character(
        tbk_casp_state
      ),
      
      levels = c(
        "Tbk1-Casp1-",
        "Tbk1+Casp1+"
      )
    )
  )


cat(
  "\nPlot 4 comparison groups:\n"
)

print(
  table(
    analysis_df$tbk_casp_group
  )
)

cat(
  "\nPlot 4 groups by sample:\n"
)

print(
  table(
    analysis_df$sample,
    analysis_df$tbk_casp_group
  )
)


############################################################
## Fisher OR
############################################################

sensor_or <- lapply(
  
  names(sensor_genes),
  
  function(g) {
    
    tab <- table(
      
      analysis_df$tbk_casp_group,
      
      analysis_df[[g]] > 0
    )
    
    cat(
      "\nContingency table for ",
      g,
      ":\n",
      sep = ""
    )
    
    print(tab)
    
    if (
      nrow(tab) < 2 ||
      ncol(tab) < 2
    ) {
      
      return(
        data.frame(
          
          gene = g,
          
          odds_ratio =
            NA_real_,
          
          ci_low =
            NA_real_,
          
          ci_high =
            NA_real_,
          
          p_value =
            NA_real_
        )
      )
    }
    
    ft <- fisher.test(
      tab
    )
    
    data.frame(
      
      gene = g,
      
      odds_ratio =
        as.numeric(
          ft$estimate
        ),
      
      ci_low =
        ft$conf.int[1],
      
      ci_high =
        ft$conf.int[2],
      
      p_value =
        ft$p.value
    )
  }
) %>%
  
  bind_rows()


sensor_or <- sensor_or %>%
  
  mutate(
    
    p_adj = p.adjust(
      p_value,
      method = "BH"
    ),
    
    gene_plot = recode(
      
      gene,
      
      "Aim2" =
        "Aim2",
      
      "Casp4" =
        "Casp4\n(Casp11)",
      
      "Nlrp3" =
        "Nlrp3"
    ),
    
    gene_plot = factor(
      
      gene_plot,
      
      levels = rev(
        c(
          "Aim2",
          "Casp4\n(Casp11)",
          "Nlrp3"
        )
      )
    ),
    
    label = paste0(
      
      "OR = ",
      round(
        odds_ratio,
        2
      ),
      
      "\nP = ",
      format_p(
        p_value
      )
    )
  )


cat(
  "\nPlot 4 OR table:\n"
)

print(
  sensor_or
)


write.csv(
  sensor_or,
  file.path(
    out_dir,
    "Plot4_sensor_enrichment_OR_table.csv"
  ),
  row.names = FALSE
)


############################################################
## Forest plot
############################################################

x_max <- max(
  sensor_or$ci_high,
  na.rm = TRUE
)

if (!is.finite(x_max)) {
  
  x_max <- max(
    sensor_or$odds_ratio,
    na.rm = TRUE
  )
}

if (!is.finite(x_max)) {
  
  x_max <- 3
}


theme_or_matched <- theme_classic(
  base_size = 14
) +
  
  theme(
    
    plot.title = element_text(
      face = "plain",
      size = 16,
      hjust = 0.5,
      color = "black"
    ),
    
    axis.title.x = element_text(
      face = "plain",
      size = 15,
      color = "black"
    ),
    
    axis.title.y =
      element_blank(),
    
    axis.text.x = element_text(
      face = "plain",
      size = 12,
      color = "black"
    ),
    
    axis.text.y = element_text(
      face = "plain",
      size = 12,
      color = "black"
    ),
    
    axis.line = element_line(
      linewidth = 0.65,
      color = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.6,
      color = "black"
    ),
    
    plot.margin = margin(
      6,
      8,
      6,
      8
    )
  )


p_sensor_or <- ggplot(
  
  sensor_or,
  
  aes(
    x = odds_ratio,
    y = gene_plot
  )
) +
  
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    color = "grey55",
    linewidth = 0.6
  ) +
  
  geom_errorbarh(
    
    aes(
      xmin = ci_low,
      xmax = ci_high
    ),
    
    height = 0.15,
    linewidth = 0.65,
    color = "black"
  ) +
  
  geom_point(
    size = 2.8,
    color = "black"
  ) +
  
  geom_text(
    
    aes(
      label = label
    ),
    
    nudge_y = -0.20,
    hjust = 0.5,
    vjust = 1,
    size = 4.0,
    lineheight = 0.9,
    color = "black"
  ) +
  
  scale_y_discrete(
    
    labels = c(
      
      "Aim2" =
        expression(
          italic("Aim2")
        ),
      
      "Casp4\n(Casp11)" =
        expression(
          italic("Casp4")~
            "("*
            italic("Casp11")*
            ")"
        ),
      
      "Nlrp3" =
        expression(
          italic("Nlrp3")
        )
    )
  ) +
  
  scale_x_continuous(
    
    limits = c(
      0,
      x_max * 1.25
    ),
    
    expand = expansion(
      mult = c(
        0.02,
        0.05
      )
    )
  ) +
  
  labs(
    
    title = expression(
      paste(
        "Sensor enrichment in ",
        italic("Tbk1"),
        "+",
        italic("Casp1"),
        "+ vs ",
        italic("Tbk1"),
        "-",
        italic("Casp1"),
        "- Kupffer-high CDAA-HFD spots"
      )
    ),
    
    x =
      "Odds ratio",
    
    y =
      NULL
  ) +
  
  theme_or_matched


print(
  p_sensor_or
)


ggsave(
  filename = file.path(
    out_dir,
    "Plot4_Sensor_enrichment_OR.png"
  ),
  plot = p_sensor_or,
  width = 7.2,
  height = 4.2,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(
    out_dir,
    "Plot4_Sensor_enrichment_OR.pdf"
  ),
  plot = p_sensor_or,
  width = 7.2,
  height = 4.2,
  bg = "white"
)


############################################################
############################################################
## PLOT 5
##
## Casp4 coupling across Tbk1/Casp1 states
############################################################


coupling_df_ref <- plot_final %>%
  
  select(
    cell,
    sample,
    condition,
    Kupffer_region,
    tbk_casp_state
  ) %>%
  
  left_join(
    sensor_expr %>%
      rownames_to_column(
        "cell"
      ) %>%
      select(
        cell,
        Casp4
      ),
    by = "cell"
  ) %>%
  
  filter(
    
    condition ==
      "CDAA-HFD",
    
    Kupffer_region ==
      "Kupffer-high",
    
    !is.na(
      tbk_casp_state
    )
  ) %>%
  
  mutate(
    Casp4_pos =
      Casp4 > 0
  )


cat(
  "\nPlot 5 state counts:\n"
)

print(
  table(
    coupling_df_ref$tbk_casp_state
  )
)

cat(
  "\nPlot 5 Casp4 positivity:\n"
)

print(
  table(
    coupling_df_ref$Casp4_pos
  )
)


############################################################
## Fisher helper
############################################################

calc_fisher_vs_ref <- function(
    df,
    target_state,
    ref_state = "Tbk1-Casp1-"
) {
  
  df_sub <- df %>%
    
    filter(
      tbk_casp_state %in%
        c(
          ref_state,
          target_state
        )
    ) %>%
    
    mutate(
      
      group = factor(
        
        ifelse(
          tbk_casp_state ==
            ref_state,
          
          ref_state,
          
          target_state
        ),
        
        levels = c(
          ref_state,
          target_state
        )
      )
    )
  
  
  cat(
    "\nComparison: ",
    target_state,
    " vs ",
    ref_state,
    "\n",
    sep = ""
  )
  
  
  print(
    table(
      df_sub$group
    )
  )
  
  
  tab <- table(
    df_sub$group,
    df_sub$Casp4_pos
  )
  
  print(tab)
  
  
  if (
    nrow(tab) < 2 ||
    ncol(tab) < 2
  ) {
    
    return(
      data.frame(
        
        comparison = paste0(
          target_state,
          " vs ",
          ref_state
        ),
        
        OR =
          NA_real_,
        
        CI_low =
          NA_real_,
        
        CI_high =
          NA_real_,
        
        p_value =
          NA_real_,
        
        n_ref =
          sum(
            df_sub$group ==
              ref_state
          ),
        
        n_target =
          sum(
            df_sub$group ==
              target_state
          )
      )
    )
  }
  
  
  ft <- fisher.test(
    tab
  )
  
  
  data.frame(
    
    comparison = paste0(
      target_state,
      " vs ",
      ref_state
    ),
    
    OR = as.numeric(
      ft$estimate
    ),
    
    CI_low =
      ft$conf.int[1],
    
    CI_high =
      ft$conf.int[2],
    
    p_value =
      ft$p.value,
    
    n_ref =
      sum(
        df_sub$group ==
          ref_state
      ),
    
    n_target =
      sum(
        df_sub$group ==
          target_state
      )
  )
}


############################################################
## Compare all three states directly
## with double-negative reference
############################################################

coupling_or_ref <- bind_rows(
  
  calc_fisher_vs_ref(
    df =
      coupling_df_ref,
    
    target_state =
      "Tbk1-Casp1+",
    
    ref_state =
      "Tbk1-Casp1-"
  ),
  
  calc_fisher_vs_ref(
    df =
      coupling_df_ref,
    
    target_state =
      "Tbk1+Casp1-",
    
    ref_state =
      "Tbk1-Casp1-"
  ),
  
  calc_fisher_vs_ref(
    df =
      coupling_df_ref,
    
    target_state =
      "Tbk1+Casp1+",
    
    ref_state =
      "Tbk1-Casp1-"
  )
) %>%
  
  mutate(
    
    p_adj = p.adjust(
      p_value,
      method = "BH"
    ),
    
    comparison = factor(
      
      comparison,
      
      levels = rev(
        c(
          "Tbk1-Casp1+ vs Tbk1-Casp1-",
          "Tbk1+Casp1- vs Tbk1-Casp1-",
          "Tbk1+Casp1+ vs Tbk1-Casp1-"
        )
      )
    ),
    
    label = paste0(
      
      "OR = ",
      round(
        OR,
        2
      ),
      
      "\nP = ",
      format_p(
        p_value
      )
    )
  )


cat(
  "\nFinal Casp4 Fisher OR table:\n"
)

print(
  coupling_or_ref
)


write.csv(
  coupling_df_ref,
  file.path(
    out_dir,
    "Plot5_Casp4_coupling_source_spots.csv"
  ),
  row.names = FALSE
)

write.csv(
  coupling_or_ref,
  file.path(
    out_dir,
    "Plot5_Casp4_coupling_OR_table.csv"
  ),
  row.names = FALSE
)


############################################################
## Forest plot
############################################################

x_max_coupling <- max(
  coupling_or_ref$CI_high,
  na.rm = TRUE
)

if (!is.finite(x_max_coupling)) {
  
  x_max_coupling <- max(
    coupling_or_ref$OR,
    na.rm = TRUE
  )
}

if (!is.finite(x_max_coupling)) {
  
  x_max_coupling <- 3
}


theme_ref_or_matched <- theme_classic(
  base_size = 14
) +
  
  theme(
    
    plot.title = element_text(
      face = "plain",
      size = 16,
      hjust = 0.5,
      color = "black"
    ),
    
    axis.title.x = element_text(
      face = "plain",
      size = 15,
      color = "black"
    ),
    
    axis.title.y =
      element_blank(),
    
    axis.text.x = element_text(
      face = "plain",
      size = 12,
      color = "black"
    ),
    
    axis.text.y = element_text(
      face = "plain",
      size = 12,
      color = "black"
    ),
    
    axis.line = element_line(
      linewidth = 0.65,
      color = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.6,
      color = "black"
    ),
    
    plot.margin = margin(
      8,
      40,
      8,
      8
    )
  )


p_casp4_coupling <- ggplot(
  
  coupling_or_ref,
  
  aes(
    x = OR,
    y = comparison
  )
) +
  
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    color = "grey55",
    linewidth = 0.6
  ) +
  
  geom_errorbarh(
    
    aes(
      xmin = CI_low,
      xmax = CI_high
    ),
    
    height = 0.15,
    linewidth = 0.65,
    color = "black"
  ) +
  
  geom_point(
    size = 2.8,
    color = "black"
  ) +
  
  geom_text(
    
    aes(
      label = label
    ),
    
    nudge_x = 0.08,
    nudge_y = -0.22,
    hjust = 0,
    vjust = 1,
    size = 3.8,
    lineheight = 0.9,
    color = "black"
  ) +
  
  scale_y_discrete(
    
    labels = c(
      
      "Tbk1-Casp1+ vs Tbk1-Casp1-" =
        
        expression(
          italic("Tbk1")*"-"*
            italic("Casp1")*"+"~
            "vs"~
            italic("Tbk1")*"-"*
            italic("Casp1")*"-"
        ),
      
      "Tbk1+Casp1- vs Tbk1-Casp1-" =
        
        expression(
          italic("Tbk1")*"+"*
            italic("Casp1")*"-"~
            "vs"~
            italic("Tbk1")*"-"*
            italic("Casp1")*"-"
        ),
      
      "Tbk1+Casp1+ vs Tbk1-Casp1-" =
        
        expression(
          italic("Tbk1")*"+"*
            italic("Casp1")*"+"~
            "vs"~
            italic("Tbk1")*"-"*
            italic("Casp1")*"-"
        )
    )
  ) +
  
  coord_cartesian(
    
    xlim = c(
      0,
      x_max_coupling * 1.35
    ),
    
    clip =
      "off"
  ) +
  
  labs(
    
    title = expression(
      paste(
        italic("Casp4"),
        " coupling across ",
        italic("Tbk1"),
        "/",
        italic("Casp1"),
        " states in CDAA-HFD Kupffer-high spots"
      )
    ),
    
    x =
      "Odds ratio",
    
    y =
      NULL
  ) +
  
  theme_ref_or_matched


print(
  p_casp4_coupling
)


ggsave(
  filename = file.path(
    out_dir,
    "Plot5_Casp4_coupling_Tbk1_Casp1_states.png"
  ),
  plot = p_casp4_coupling,
  width = 7.2,
  height = 4.2,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(
    out_dir,
    "Plot5_Casp4_coupling_Tbk1_Casp1_states.pdf"
  ),
  plot = p_casp4_coupling,
  width = 7.2,
  height = 4.2,
  bg = "white"
)
