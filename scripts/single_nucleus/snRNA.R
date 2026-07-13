rm(list = ls())
gc()
set.seed(1234)

############################################################
## Human MASLD snRNA-seq
## Data import
## Macrophage extraction
## Macrophage re-clustering
## Kupffer-cell identification
############################################################


############################################################
## 0. Load packages
############################################################

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)
library(ggrepel)
library(tibble)
library(rstatix)
library(igraph)
library(lme4)
library(emmeans)
library(grid)


############################################################
## 1. Import data
############################################################

obj <- readRDS(
  "/Users/yinhuang/Downloads/MASLD_snRNA_seq_seurat_v4.rds"
)

############################################################
## 1.1 Global definitions used by all downstream analyses
############################################################

disease_levels <- c(
  "control",
  "MASLD",
  "eMASH",
  "aMASH"
)

state_levels <- c(
  "Other",
  "TBK1+",
  "CASP1+",
  "TBK1+CASP1+"
)

disease_cols <- c(
  control = "#BDBDBD",
  MASLD   = "#E69F00",
  eMASH   = "#56B4E9",
  aMASH   = "#C41E1E"
)

state_cols <- c(
  "Other" = "#BDBDBD",
  "TBK1+" = "#E69F00",
  "CASP1+" = "#56B4E9",
  "TBK1+CASP1+" = "#C41E1E"
)

obj$Disease_group <- as.character(obj$Disease_group)
obj$Disease_group[obj$Disease_group == "MASL"] <- "MASLD"
obj$Disease_group <- factor(
  obj$Disease_group,
  levels = disease_levels
)

DefaultAssay(obj) <- "RNA"


############################################################
## 2. Check metadata
############################################################

colnames(obj@meta.data)

head(obj@meta.data)

dim(obj@meta.data)

grep(
  "group|status|sample|diag|disease|condition",
  colnames(obj@meta.data),
  value = TRUE,
  ignore.case = TRUE
)

table(obj$Disease_group)

unique(obj$Disease_group)

sort(unique(obj$Cell_type_broad))

grep(
  "macro|inflam",
  unique(obj$Cell_type_detailed),
  value = TRUE,
  ignore.case = TRUE
)


############################################################
## 3. Extract all macrophages
############################################################


Idents(obj) <- "Cell_type_broad"

macro <- subset(
  obj,
  idents = "Macrophage"
)

DefaultAssay(macro) <- "RNA"


############################################################
## 4. Re-cluster macrophages
############################################################

macro <- NormalizeData(macro)

macro <- FindVariableFeatures(macro)

macro <- ScaleData(macro)

macro <- RunPCA(macro)

macro <- FindNeighbors(
  macro,
  dims = 1:20
)

macro <- FindClusters(
  macro,
  resolution = 0.4
)

macro <- RunUMAP(
  macro,
  dims = 1:20
)


############################################################
## 5. Check macrophage re-clustering
############################################################

DimPlot(
  macro,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
)

table(macro$seurat_clusters)


############################################################
## 6. Inspect canonical Kupffer-cell markers
############################################################

FeaturePlot(
  macro,
  features = c(
    "MARCO",
    "VSIG4",
    "C1QA",
    "C1QB",
    "C1QC"
  ),
  reduction = "umap",
  ncol = 3
)


DotPlot(
  macro,
  features = c(
    "MARCO",
    "VSIG4",
    "C1QA",
    "C1QB",
    "C1QC"
  )
) +
  RotatedAxis()


VlnPlot(
  macro,
  features = c(
    "MARCO",
    "VSIG4",
    "C1QA",
    "C1QB",
    "C1QC"
  ),
  group.by = "seurat_clusters",
  pt.size = 0
)


############################################################
## 7. Define Kupffer-cell clusters
############################################################

kupffer_clusters <- c(
  "0",
  "1",
  "7",
  "8"
)


############################################################
## 7.1 Define Kupffer and non-Kupffer macrophages once
############################################################

kupffer_cells <- colnames(macro)[
  macro$seurat_clusters %in% kupffer_clusters
]

macrophage_non_kupffer_cells <- colnames(macro)[
  !(macro$seurat_clusters %in% kupffer_clusters)
]

macro$macrophage_subtype <- ifelse(
  colnames(macro) %in% kupffer_cells,
  "Kupffer",
  "Macrophage_non_Kupffer"
)

macro$macrophage_subtype <- factor(
  macro$macrophage_subtype,
  levels = c("Kupffer", "Macrophage_non_Kupffer")
)

obj$cell_type_with_kupffer_updated <- as.character(
  obj$Cell_type_broad
)

obj$cell_type_with_kupffer_updated[
  colnames(obj) %in% kupffer_cells
] <- "Kupffer"

obj$cell_type_with_kupffer_updated[
  colnames(obj) %in% macrophage_non_kupffer_cells
] <- "Macrophage_non_Kupffer"

obj$cell_type_with_kupffer_updated <- factor(
  obj$cell_type_with_kupffer_updated
)

table(macro$macrophage_subtype)
table(obj$cell_type_with_kupffer_updated, useNA = "ifany")

############################################################
## 8. Extract Kupffer cells
############################################################

kup <- subset(
  macro,
  subset = seurat_clusters %in% kupffer_clusters
)

DefaultAssay(kup) <- "RNA"

macro$Disease_group <- factor(macro$Disease_group, levels = disease_levels)
kup$Disease_group <- factor(kup$Disease_group, levels = disease_levels)

############################################################
## 9. Check extracted Kupffer cells
############################################################

kup

ncol(kup)

table(kup$Disease_group)

table(kup$seurat_clusters)


############################################################
## 10. Confirm that the macrophage UMAP is retained in kup
############################################################

Reductions(kup)

if (!"umap" %in% Reductions(kup)) {
  stop(
    "UMAP reduction was not retained in the Kupffer object."
  )
}


############################################################
## Step 1 
############################################################


DefaultAssay(kup) <- "RNA"


############################################################
## 1.1 Extract UMAP and expression
############################################################

umap_df <- as.data.frame(
  Embeddings(kup, "umap")
)

colnames(umap_df) <- c(
  "UMAP_1",
  "UMAP_2"
)

umap_df$cell <- rownames(umap_df)

expr_df <- FetchData(
  kup,
  vars = c(
    "TBK1",
    "CASP1"
  )
)

expr_df$cell <- rownames(expr_df)

umap_df <- umap_df %>%
  left_join(
    expr_df,
    by = "cell"
  )


############################################################
## 1.2 Re-define TBK1/CASP1 states
############################################################

umap_df <- umap_df %>%
  mutate(
    TBK1_pos = TBK1 > 0,
    CASP1_pos = CASP1 > 0,
    
    state_new = case_when(
      TBK1_pos & CASP1_pos ~
        "TBK1+CASP1+",
      
      TBK1_pos & !CASP1_pos ~
        "TBK1+",
      
      !TBK1_pos & CASP1_pos ~
        "CASP1+",
      
      TRUE ~
        "Other"
    )
  )

umap_df$state_new <- factor(
  umap_df$state_new,
  levels = state_levels
)

print(
  umap_df %>%
    count(state_new) %>%
    mutate(
      percent = 100 * n / sum(n)
    )
)


############################################################
## 1.3 Write TBK1/CASP1 states back to kup
############################################################

state_match <- match(
  colnames(kup),
  umap_df$cell
)

kup$TBK1_expr <- umap_df$TBK1[state_match]
kup$CASP1_expr <- umap_df$CASP1[state_match]

kup$TBK1_pos <- umap_df$TBK1_pos[state_match]
kup$CASP1_pos <- umap_df$CASP1_pos[state_match]

kup$TBK1_positive <- kup$TBK1_pos
kup$CASP1_positive <- kup$CASP1_pos

kup$state_plot <- factor(
  umap_df$state_new[state_match],
  levels = state_levels
)

table(kup$state_plot)


############################################################
## 1.4 Cap expression for display only
############################################################

expr_cap <- 2

umap_df <- umap_df %>%
  mutate(
    TBK1_plot = pmin(
      TBK1,
      expr_cap
    ),
    
    CASP1_plot = pmin(
      CASP1,
      expr_cap
    )
  )


############################################################
## 1.5 Randomize state plotting order
############################################################


umap_df_state <- umap_df %>%
  slice_sample(
    prop = 1
  )


############################################################
## 1.6
############################################################

theme_umap_consistent <- theme_classic(
  base_size = 14
) +
  theme(
    plot.title = element_text(
      face = "plain",
      size = 18,
      hjust = 0,
      color = "black"
    ),
    
    axis.title.x = element_text(
      face = "bold",
      size = 16,
      color = "black"
    ),
    
    axis.title.y = element_text(
      face = "bold",
      size = 16,
      color = "black"
    ),
    
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
      linewidth = 0.8,
      color = "black"
    ),
    
    legend.title = element_blank(),
    
    legend.text = element_text(
      face = "plain",
      size = 12,
      color = "black"
    ),
    
    legend.key.size = unit(
      0.45,
      "cm"
    ),
    
    plot.margin = margin(
      6,
      6,
      6,
      6
    )
  )


############################################################
## 1.7 TBK1 expression UMAP
############################################################

p_umap_TBK1 <- ggplot(
  umap_df,
  aes(
    x = UMAP_1,
    y = UMAP_2
  )
) +
  geom_point(
    aes(
      color = TBK1_plot
    ),
    size = 0.28,
    alpha = 0.85
  ) +
  scale_color_gradient(
    low = "grey90",
    high = "#FF3B30",
    limits = c(
      0,
      expr_cap
    ),
    breaks = c(
      0,
      0.5,
      1.0,
      1.5,
      2.0
    ),
    oob = scales::squish,
    name = NULL
  ) +
  labs(
    title = expression(
      italic("TBK1")
    ),
    x = "UMAP_1",
    y = "UMAP_2"
  ) +
  theme_umap_consistent +
  theme(
    legend.text = element_text(
      size = 11
    ),
    
    legend.key.height = unit(
      0.45,
      "cm"
    )
  )


############################################################
## 1.8 CASP1 expression UMAP
############################################################

p_umap_CASP1 <- ggplot(
  umap_df,
  aes(
    x = UMAP_1,
    y = UMAP_2
  )
) +
  geom_point(
    aes(
      color = CASP1_plot
    ),
    size = 0.28,
    alpha = 0.85
  ) +
  scale_color_gradient(
    low = "grey90",
    high = "#FF3B30",
    limits = c(
      0,
      expr_cap
    ),
    breaks = c(
      0,
      0.5,
      1.0,
      1.5,
      2.0
    ),
    oob = scales::squish,
    name = NULL
  ) +
  labs(
    title = expression(
      italic("CASP1")
    ),
    x = "UMAP_1",
    y = "UMAP_2"
  ) +
  theme_umap_consistent +
  theme(
    legend.text = element_text(
      size = 11
    ),
    
    legend.key.height = unit(
      0.45,
      "cm"
    )
  )


############################################################
## 1.9 TBK1/CASP1 state UMAP
############################################################

p_umap_state <- ggplot(
  umap_df_state,
  aes(
    x = UMAP_1,
    y = UMAP_2,
    color = state_new
  )
) +
  geom_point(
    size = 0.25,
    alpha = 0.85
  ) +
  scale_color_manual(
    values = state_cols,
    
    breaks = c(
      "Other",
      "TBK1+",
      "CASP1+",
      "TBK1+CASP1+"
    ),
    
    labels = c(
      "Other",
      
      expression(
        italic("TBK1") * "+"
      ),
      
      expression(
        italic("CASP1") * "+"
      ),
      
      expression(
        italic("TBK1") * "+" *
          italic("CASP1") * "+"
      )
    ),
    
    name = NULL
  ) +
  labs(
    title = expression(
      paste(
        italic("TBK1"),
        "/",
        italic("CASP1"),
        " states"
      )
    ),
    x = "UMAP_1",
    y = "UMAP_2"
  ) +
  theme_umap_consistent +
  theme(
    legend.position = "right",
    
    legend.text = element_text(
      size = 12,
      color = "black"
    )
  )


############################################################
## 1.10 Combine
############################################################

p_umap_TBK1_CASP1_state_consistent <-
  p_umap_TBK1 +
  p_umap_CASP1 +
  p_umap_state +
  plot_layout(
    ncol = 3,
    widths = c(
      1,
      1,
      1.25
    )
  )

p_umap_TBK1_CASP1_state_consistent


############################################################
## 1.11 Save
############################################################

ggsave(
  "Step1_Kupffer_TBK1_CASP1_state_UMAP_consistent_style.pdf",
  p_umap_TBK1_CASP1_state_consistent,
  width = 10.5,
  height = 3.8
)

ggsave(
  "Step1_Kupffer_TBK1_CASP1_state_UMAP_consistent_style.png",
  p_umap_TBK1_CASP1_state_consistent,
  width = 10.5,
  height = 3.8,
  dpi = 600
)

write.csv(
  umap_df,
  "Step1_Kupffer_TBK1_CASP1_state_UMAP_consistent_style_source_data.csv",
  row.names = FALSE
)


############################################################
## 1.12 Save final Kupffer object
############################################################

saveRDS(
  kup,
  "Kupffer_from_macrophage_clusters_with_TBK1_CASP1_states.rds"
)


############################################################
## Step2
## Shared Kupffer-cell inflammatory pseudotime
## Then extract MASLD-specific gene dynamics
############################################################


############################################################
## 2.1 Check input object
############################################################

if (!exists("kup")) {
  stop(
    paste0(
      "Object 'kup' is not found. ",
      "Please run the Kupffer extraction and Step1 first."
    )
  )
}

DefaultAssay(kup) <- "RNA"

cat("Current Kupffer object:\n")
print(kup)

cat("Number of Kupffer cells:\n")
print(ncol(kup))


############################################################
## 2.2 Confirm inherited disease labels and TBK1/CASP1 states
############################################################

cat("Disease groups in current Kupffer object:\n")
print(table(kup$Disease_group, useNA = "ifany"))

cat("TBK1/CASP1 states by disease group:\n")
print(table(kup$Disease_group, kup$state_plot))

if (sum(kup$Disease_group == "MASLD", na.rm = TRUE) == 0) {
  stop("No MASLD cells were found in kup.")
}

############################################################
## 2.3 Confirm PCA is available
############################################################

print(Reductions(kup))

if (!"pca" %in% Reductions(kup)) {
  
  message(
    paste0(
      "PCA was not found in kup. ",
      "PCA will now be calculated using the existing Kupffer object."
    )
  )
  
  kup <- FindVariableFeatures(
    kup
  )
  
  kup <- ScaleData(
    kup
  )
  
  kup <- RunPCA(
    kup,
    npcs = 30
  )
}

if (!"umap" %in% Reductions(kup)) {
  
  message(
    paste0(
      "UMAP was not found in kup. ",
      "UMAP will be calculated from the existing PCA."
    )
  )
  
  kup <- RunUMAP(
    kup,
    dims = 1:20
  )
}


############################################################
## 2.4 Define genes for inflammatory trajectory score
############################################################

trajectory_score_genes <- c(
  "TBK1",
  "IRF1",
  "STAT1",
  "IFIH1",
  "CASP1",
  "CASP4",
  "GSDMD",
  "IL1B"
)

trajectory_score_genes <- trajectory_score_genes[
  trajectory_score_genes %in%
    rownames(kup)
]

cat("Genes used for trajectory score:\n")
print(trajectory_score_genes)

if (
  length(trajectory_score_genes) == 0
) {
  stop(
    "None of the trajectory-score genes were found in kup."
  )
}


############################################################
## 2.5 Calculate inflammatory trajectory score
############################################################


kup <- AddModuleScore(
  object = kup,
  features = list(
    trajectory_score_genes
  ),
  name = "Trajectory_score_"
)

score_col <- grep(
  "^Trajectory_score_",
  colnames(kup@meta.data),
  value = TRUE
)[1]

cat("Trajectory score column:\n")
print(score_col)

if (
  is.na(score_col) ||
  length(score_col) == 0
) {
  stop(
    "Trajectory score column was not generated."
  )
}


############################################################
## 2.6 Extract PCA coordinates from ALL Kupffer cells
############################################################

pca_embeddings <- Embeddings(
  kup,
  reduction = "pca"
)

n_pcs_available <- ncol(
  pca_embeddings
)

n_pcs_use <- min(
  20,
  n_pcs_available
)

cat("Number of PCs used:\n")
print(n_pcs_use)

pca_mat <- pca_embeddings[
  ,
  seq_len(n_pcs_use),
  drop = FALSE
]

cell_names <- rownames(
  pca_mat
)

if (
  !identical(
    sort(cell_names),
    sort(colnames(kup))
  )
) {
  stop(
    "PCA cell names do not match the current Kupffer object."
  )
}


############################################################
## 2.7 Calculate distances in PCA space
############################################################

dist_mat <- as.matrix(
  dist(pca_mat)
)

rownames(dist_mat) <- cell_names
colnames(dist_mat) <- cell_names

dim(dist_mat)


############################################################
## 2.8 Build k-nearest-neighbor graph
############################################################

k <- 15

if (
  nrow(dist_mat) <= k
) {
  stop(
    "The number of Kupffer cells is too small for k = 15."
  )
}

knn_edges <- lapply(
  seq_len(nrow(dist_mat)),
  function(i) {
    
    nearest_cells <- order(
      dist_mat[i, ]
    )[2:(k + 1)]
    
    data.frame(
      from = cell_names[i],
      to = cell_names[nearest_cells],
      weight = dist_mat[
        i,
        nearest_cells
      ],
      stringsAsFactors = FALSE
    )
  }
)

edge_df <- bind_rows(
  knn_edges
)

g <- graph_from_data_frame(
  d = edge_df,
  directed = FALSE,
  vertices = data.frame(
    name = cell_names,
    stringsAsFactors = FALSE
  )
)

cat("Graph summary:\n")
print(g)


############################################################
## 2.9 Check graph connectivity
############################################################

graph_components <- components(
  g
)

cat("Number of graph components:\n")
print(graph_components$no)

cat("Size of graph components:\n")
print(graph_components$csize)


############################################################
## 2.10 Choose root cells
############################################################

meta_tmp <- kup@meta.data

meta_tmp <- meta_tmp[
  colnames(kup),
  ,
  drop = FALSE
]

root_candidates <- rownames(meta_tmp)[
  meta_tmp$Disease_group == "control" &
    meta_tmp$state_plot == "Other"
]

cat("Initial control + Other root candidates:\n")
print(length(root_candidates))

if (
  length(root_candidates) < 5
) {
  
  root_candidates <- rownames(meta_tmp)[
    meta_tmp$Disease_group == "control"
  ]
  
  message(
    paste0(
      "Fewer than 5 control-Other cells were found. ",
      "All control Kupffer cells will be considered as root candidates."
    )
  )
}

if (
  length(root_candidates) < 5
) {
  
  root_candidates <- rownames(
    meta_tmp
  )
  
  message(
    paste0(
      "Fewer than 5 control Kupffer cells were found. ",
      "All Kupffer cells will be considered as root candidates."
    )
  )
}

root_candidates <- root_candidates[
  !is.na(
    meta_tmp[
      root_candidates,
      score_col
    ]
  )
]

root_candidates <- root_candidates[
  order(
    meta_tmp[
      root_candidates,
      score_col
    ]
  )
]

root_cells <- root_candidates[
  seq_len(
    min(
      20,
      length(root_candidates)
    )
  )
]

cat("Number of selected root cells:\n")
print(length(root_cells))

cat("Root cells overlapping graph vertices:\n")

print(
  sum(
    root_cells %in%
      V(g)$name
  )
)

if (
  sum(
    root_cells %in%
    V(g)$name
  ) == 0
) {
  stop(
    paste0(
      "No root cells overlap with graph vertices. ",
      "Please check cell names."
    )
  )
}

root_cells <- root_cells[
  root_cells %in%
    V(g)$name
]


############################################################
## 2.11 Calculate graph-distance pseudotime
############################################################

dist_from_roots <- distances(
  graph = g,
  v = root_cells,
  to = V(g),
  weights = E(g)$weight
)

pseudo_raw <- apply(
  dist_from_roots,
  2,
  min,
  na.rm = TRUE
)

pseudo_raw[
  !is.finite(pseudo_raw)
] <- NA_real_

pseudo_raw <- pseudo_raw[
  colnames(kup)
]

names(pseudo_raw) <- colnames(
  kup
)

if (
  all(is.na(pseudo_raw))
) {
  stop(
    "All pseudotime distances are NA."
  )
}


############################################################
## 2.12 Scale pseudotime to 0-1
############################################################

pseudo_min <- min(
  pseudo_raw,
  na.rm = TRUE
)

pseudo_max <- max(
  pseudo_raw,
  na.rm = TRUE
)

if (
  pseudo_max == pseudo_min
) {
  stop(
    "Pseudotime has no variation."
  )
}

pseudo_scaled <- (
  pseudo_raw - pseudo_min
) / (
  pseudo_max - pseudo_min
)

names(pseudo_scaled) <- names(
  pseudo_raw
)


############################################################
## 2.13 Orient pseudotime
############################################################

score_vec <- kup@meta.data[
  colnames(kup),
  score_col
]

cor_pt_score <- suppressWarnings(
  cor(
    pseudo_scaled,
    score_vec,
    method = "spearman",
    use = "complete.obs"
  )
)

cat(
  paste0(
    "Spearman correlation between shared pseudotime ",
    "and inflammatory score:\n"
  )
)

print(cor_pt_score)

if (
  !is.na(cor_pt_score) &&
  cor_pt_score < 0
) {
  
  pseudo_scaled <- 1 -
    pseudo_scaled
  
  message(
    paste0(
      "Pseudotime was reversed so that the inflammatory ",
      "trajectory score increases along pseudotime."
    )
  )
}


############################################################
## 2.14 Add shared pseudotime to Kupffer object
############################################################

kup <- AddMetaData(
  object = kup,
  metadata = data.frame(
    shared_pseudotime = as.numeric(
      pseudo_scaled
    ),
    row.names = names(
      pseudo_scaled
    )
  )
)

cat("Shared pseudotime summary:\n")

print(
  summary(
    kup$shared_pseudotime
  )
)

cat("Pseudotime-related metadata columns:\n")

print(
  grep(
    "pseudo|time",
    colnames(kup@meta.data),
    value = TRUE,
    ignore.case = TRUE
  )
)


############################################################
## 2.15 Define genes shown in the final figure
############################################################

genes_use <- c(
  "IRF1",
  "TBK1",
  "IL18",
  "STAT1",
  "CASP4",
  "CASP1",
  "IL1B",
  "GSDMD"
)

genes_use <- genes_use[
  genes_use %in%
    rownames(kup)
]

missing_trajectory_genes <- setdiff(
  c(
    "IRF1",
    "TBK1",
    "IL18",
    "STAT1",
    "CASP4",
    "CASP1",
    "IL1B",
    "GSDMD"
  ),
  genes_use
)

cat("Genes used for plotting:\n")
print(genes_use)

cat("Missing genes:\n")
print(missing_trajectory_genes)

if (
  length(genes_use) == 0
) {
  stop(
    "None of the requested trajectory genes were found."
  )
}


############################################################
## 2.16 Prepare ALL-disease shared pseudotime data
############################################################

traj_df <- FetchData(
  kup,
  vars = c(
    genes_use,
    "shared_pseudotime",
    "Disease_group",
    "state_plot",
    "DonorID"
  )
) %>%
  filter(
    !is.na(shared_pseudotime),
    !is.na(Disease_group)
  )

cat("Cells used in all-disease pseudotime analysis:\n")
print(nrow(traj_df))

cat("Cells used by disease group:\n")
print(table(traj_df$Disease_group))


############################################################
## 2.17 Convert all-disease data to long format
############################################################

traj_long <- traj_df %>%
  pivot_longer(
    cols = all_of(
      genes_use
    ),
    names_to = "gene",
    values_to = "expr"
  )


############################################################
## 2.18 Summarize gene expression across all disease stages
############################################################

traj_bin_stage <- traj_long %>%
  mutate(
    pt_bin = cut(
      shared_pseudotime,
      breaks = seq(
        0,
        1,
        length.out = 26
      ),
      include.lowest = TRUE
    )
  ) %>%
  group_by(
    Disease_group,
    gene,
    pt_bin
  ) %>%
  summarise(
    mean_expr = mean(
      expr,
      na.rm = TRUE
    ),
    mean_pt = mean(
      shared_pseudotime,
      na.rm = TRUE
    ),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  filter(
    n_cells >= 10
  )

write.csv(
  traj_bin_stage,
  "Step4_all_disease_stage_stratified_pseudotime_gene_dynamics.csv",
  row.names = FALSE
)


############################################################
## 2.19 Plot all-disease pseudotime dynamics
############################################################

p_stage_gene <- ggplot(
  traj_bin_stage,
  aes(
    x = mean_pt,
    y = mean_expr,
    color = Disease_group,
    group = Disease_group
  )
) +
  geom_line(
    linewidth = 1.1
  ) +
  geom_point(
    size = 1.4,
    alpha = 0.85
  ) +
  facet_wrap(
    ~ gene,
    scales = "free_y",
    ncol = 4
  ) +
  scale_color_manual(
    values = disease_cols
  ) +
  labs(
    title =
      "Stage-stratified gene dynamics along shared Kupffer-cell pseudotime",
    x = "Shared inflammatory pseudotime",
    y = "Mean normalized expression"
  ) +
  theme_classic(
    base_size = 14
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    strip.text = element_text(
      face = "bold"
    ),
    axis.text.x = element_text(
      face = "bold"
    )
  )

p_stage_gene


############################################################
## 2.20 Save all-disease check figure
############################################################

ggsave(
  "Step4_all_disease_gene_dynamics_shared_pseudotime.pdf",
  p_stage_gene,
  width = 12,
  height = 7
)

ggsave(
  "Step4_all_disease_gene_dynamics_shared_pseudotime.png",
  p_stage_gene,
  width = 12,
  height = 7,
  dpi = 600
)


############################################################
## 2.21 Extract MASLD cells
############################################################

masld_cells <- rownames(
  kup@meta.data
)[
  kup$Disease_group == "MASLD" &
    !is.na(
      kup$shared_pseudotime
    )
]

cat(
  "Number of MASLD cells with shared pseudotime:\n"
)

print(
  length(
    masld_cells
  )
)

if (
  length(masld_cells) == 0
) {
  stop(
    "No MASLD cells with shared pseudotime were found."
  )
}


############################################################
## 2.22 Extract MASLD gene expression and pseudotime
############################################################

expr_masld <- FetchData(
  kup,
  vars = c(
    genes_use,
    "shared_pseudotime"
  ),
  cells = masld_cells
)

expr_masld <- expr_masld %>%
  rename(
    pseudotime =
      shared_pseudotime
  ) %>%
  filter(
    !is.na(pseudotime)
  ) %>%
  arrange(
    pseudotime
  )

cat("MASLD shared-pseudotime range:\n")

print(
  range(
    expr_masld$pseudotime,
    na.rm = TRUE
  )
)


############################################################
## 2.23 Bin MASLD shared pseudotime
############################################################

n_bins <- 20

expr_masld$bin <- cut(
  expr_masld$pseudotime,
  breaks = n_bins,
  labels = FALSE
)

plot_df_masld <- expr_masld %>%
  group_by(
    bin
  ) %>%
  summarise(
    pseudotime = mean(
      pseudotime,
      na.rm = TRUE
    ),
    across(
      all_of(genes_use),
      mean,
      na.rm = TRUE
    ),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  filter(
    n_cells >= 5
  ) %>%
  pivot_longer(
    cols = all_of(
      genes_use
    ),
    names_to = "gene",
    values_to = "expression"
  )

if (
  nrow(plot_df_masld) == 0
) {
  stop(
    paste0(
      "No MASLD pseudotime bins remained after requiring ",
      "at least 5 cells per bin."
    )
  )
}


############################################################
## 2.24 Set facet order
############################################################

gene_order <- c(
  "CASP1",
  "CASP4",
  "GSDMD",
  "IL18",
  "IL1B",
  "IRF1",
  "STAT1",
  "TBK1"
)

gene_order <- gene_order[
  gene_order %in%
    unique(
      plot_df_masld$gene
    )
]

plot_df_masld$gene <- factor(
  plot_df_masld$gene,
  levels = gene_order
)

cat("Final gene order:\n")
print(gene_order)

cat("Number of retained MASLD pseudotime bins:\n")

print(
  length(
    unique(
      plot_df_masld$bin
    )
  )
)

cat("Retained MASLD bins:\n")

print(
  plot_df_masld %>%
    distinct(
      bin,
      pseudotime,
      n_cells
    ) %>%
    arrange(
      bin
    )
)


############################################################
## 2.25 Export MASLD source table
############################################################

write.csv(
  plot_df_masld,
  "Step4_MASLD_gene_dynamics_after_all_disease_shared_pseudotime.csv",
  row.names = FALSE
)


############################################################
## 2.26 Prepare italic and bold gene labels
############################################################

gene_label_map <- setNames(
  paste0(
    "bolditalic('",
    gene_order,
    "')"
  ),
  gene_order
)


############################################################
## 2.27 Plot MASLD-specific gene dynamics
############################################################

p_masld_traj <- ggplot(
  plot_df_masld,
  aes(
    x = pseudotime,
    y = expression
  )
) +
  geom_line(
    color = "#E69F00",
    linewidth = 1.4
  ) +
  geom_point(
    color = "#E69F00",
    size = 2
  ) +
  facet_wrap(
    ~ gene,
    scales = "free_y",
    ncol = 4,
    labeller = as_labeller(
      gene_label_map,
      label_parsed
    )
  ) +
  labs(
    title =
      "MASLD-specific gene dynamics along shared inflammatory pseudotime",
    x =
      "Shared inflammatory pseudotime",
    y =
      "Mean normalized expression"
  ) +
  theme_classic(
    base_size = 14
  ) +
  theme(
    plot.title = element_text(
      face = "plain",
      size = 18,
      hjust = 0.5,
      color = "black"
    ),
    
    axis.title.x = element_text(
      face = "plain",
      size = 16,
      color = "black"
    ),
    
    axis.title.y = element_text(
      face = "plain",
      size = 16,
      color = "black"
    ),
    
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
    
    strip.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.8
    ),
    
    strip.text = element_text(
      face = "bold.italic",
      size = 13,
      color = "black"
    ),
    
    panel.spacing = unit(
      0.75,
      "lines"
    ),
    
    plot.margin = margin(
      8,
      8,
      8,
      8
    )
  )

p_masld_traj


############################################################
## 2.28 Save final MASLD figure
############################################################

ggsave(
  "Step4_MASLD_gene_dynamics_shared_inflammatory_pseudotime.pdf",
  p_masld_traj,
  width = 10,
  height = 6
)

ggsave(
  "Step4_MASLD_gene_dynamics_shared_inflammatory_pseudotime.png",
  p_masld_traj,
  width = 10,
  height = 6,
  dpi = 600
)


############################################################
## 2.29 Export pseudotime metadata
############################################################

pseudotime_output <- kup@meta.data %>%
  select(
    DonorID,
    Disease_group,
    state_plot,
    all_of(score_col),
    shared_pseudotime
  )

write.csv(
  pseudotime_output,
  "Step4_all_Kupffer_shared_inflammatory_pseudotime_values.csv",
  row.names = FALSE
)


############################################################
## 2.30 Save updated Kupffer object
############################################################

saveRDS(
  kup,
  "Step4_Kupffer_with_shared_inflammatory_pseudotime.rds"
)


############################################################
## Figure 1A and Figure 1B
############################################################

DefaultAssay(obj) <- "RNA"
DefaultAssay(macro) <- "RNA"
DefaultAssay(kup) <- "RNA"

############################################################
## Figure A
############################################################

DefaultAssay(obj) <- "RNA"


############################################################
## A.2 Check cell-type column
############################################################

if (
  !"cell_type_with_kupffer_updated" %in%
  colnames(obj@meta.data)
) {
  stop(
    paste0(
      "cell_type_with_kupffer_updated is not found. ",
      "Please run the Kupffer annotation step first."
    )
  )
}

table(
  obj$cell_type_with_kupffer_updated
)


############################################################
## A.3 Remove useless cell groups
############################################################

remove_celltypes <- c(
  "Unidentified",
  "Suspected_doublet"
)

cat("Cell types to remove:\n")
print(remove_celltypes)

cat("Before filtering:\n")

print(
  table(
    obj$cell_type_with_kupffer_updated
  )
)


############################################################
## A.4 Extract TBK1 expression and metadata
############################################################

tbk1_df <- FetchData(
  obj,
  vars = c(
    "TBK1",
    "cell_type_with_kupffer_updated",
    "Disease_group",
    "DonorID"
  )
)

colnames(tbk1_df) <- c(
  "TBK1_expr",
  "cell_type",
  "Disease_group",
  "DonorID"
)

tbk1_df <- tbk1_df %>%
  filter(
    !is.na(TBK1_expr),
    !is.na(cell_type),
    !is.na(Disease_group),
    !is.na(DonorID),
    !(cell_type %in%
        remove_celltypes)
  ) %>%
  mutate(
    Disease_group = factor(
      Disease_group,
      levels = disease_levels
    ),
    cell_type = as.character(
      cell_type
    ),
    TBK1_positive =
      TBK1_expr > 0
  )

cat("After filtering:\n")
print(table(tbk1_df$cell_type))

write.csv(
  tbk1_df,
  "FigureA_TBK1_single_cell_expression_cleaned_with_Kupffer.csv",
  row.names = FALSE
)


############################################################
## A.5 Donor-level summary
############################################################

tbk1_donor_celltype <- tbk1_df %>%
  group_by(
    DonorID,
    Disease_group,
    cell_type
  ) %>%
  summarise(
    mean_expr = mean(
      TBK1_expr,
      na.rm = TRUE
    ),
    median_expr = median(
      TBK1_expr,
      na.rm = TRUE
    ),
    frac_TBK1_positive = mean(
      TBK1_positive,
      na.rm = TRUE
    ),
    n_cells = n(),
    .groups = "drop"
  )

write.csv(
  tbk1_donor_celltype,
  "FigureA_TBK1_donor_level_expression_cleaned_with_Kupffer.csv",
  row.names = FALSE
)


############################################################
## A.6 Heatmap source table
############################################################

tbk1_heat_base <- tbk1_donor_celltype %>%
  group_by(
    Disease_group,
    cell_type
  ) %>%
  summarise(
    mean_expr_group = mean(
      mean_expr,
      na.rm = TRUE
    ),
    median_expr_group = median(
      mean_expr,
      na.rm = TRUE
    ),
    n_donors = n_distinct(
      DonorID
    ),
    .groups = "drop"
  )


############################################################
## A.7 Cell-type order
## Ordered by overall raw mean TBK1 expression
############################################################

celltype_order <- tbk1_heat_base %>%
  group_by(cell_type) %>%
  summarise(
    overall_mean_expr = mean(
      mean_expr_group,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  arrange(
    overall_mean_expr
  ) %>%
  pull(
    cell_type
  )

tbk1_heat_base$cell_type <- factor(
  tbk1_heat_base$cell_type,
  levels = celltype_order
)

tbk1_heat_base$Disease_group <- factor(
  tbk1_heat_base$Disease_group,
  levels = disease_levels
)

write.csv(
  tbk1_heat_base,
  "FigureA_TBK1_raw_mean_expression_heatmap_source_table.csv",
  row.names = FALSE
)


############################################################
## A.8 Define enhanced color limits
############################################################

tbk1_color_limits <- quantile(
  tbk1_heat_base$mean_expr_group,
  probs = c(
    0.05,
    0.95
  ),
  na.rm = TRUE
)

tbk1_color_limits


############################################################
## A.9 Unified theme for heatmap
############################################################

final_theme_big <- theme_classic(
  base_size = 16
) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 22,
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      size = 16,
      hjust = 0.5
    ),
    axis.title = element_blank(),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold",
      size = 20,
      color = "black"
    ),
    
    axis.text.y = element_text(
      face = "plain",
      size = 20,
      color = "black"
    ),
    
    legend.title = element_text(
      face = "bold",
      size = 15
    ),
    
    legend.text = element_text(
      size = 14
    ),
    
    plot.margin = margin(
      8,
      8,
      8,
      8
    )
  )


############################################################
## A.10 Plot TBK1 heatmap
############################################################

p_figure_A <- ggplot(
  tbk1_heat_base,
  aes(
    x = Disease_group,
    y = cell_type,
    fill = mean_expr_group
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.6
  ) +
  scale_fill_gradientn(
    colors = c(
      "white",
      "#FEE5D9",
      "#FCAE91",
      "#FB6A4A",
      "#CB181D"
    ),
    limits = c(
      tbk1_color_limits[1],
      tbk1_color_limits[2]
    ),
    oob = scales::squish,
    name = expression(
      atop(
        paste(
          "Mean ",
          italic("TBK1")
        ),
        "expression"
      )
    )
  ) +
  labs(
    title = expression(
      paste(
        "Disease-stage dynamics of ",
        italic("TBK1"),
        " expression"
      )
    ),
    subtitle =
      "Raw donor-level mean expression",
    x = "",
    y = ""
  ) +
  final_theme_big

p_figure_A


############################################################
## A.11 Save Figure A
############################################################

ggsave(
  "FigureA_TBK1_raw_mean_expression_heatmap.pdf",
  p_figure_A,
  width = 7.8,
  height = 6.4
)

ggsave(
  "FigureA_TBK1_raw_mean_expression_heatmap.png",
  p_figure_A,
  width = 7.8,
  height = 6.4,
  dpi = 600
)


############################################################
## Figure B
############################################################

DefaultAssay(kup) <- "RNA"


############################################################
## B.2 Use the existing TBK1-positive definition
## Includes TBK1+CASP1- and TBK1+CASP1+ cells
############################################################

kup$TBK1_status <- factor(
  ifelse(kup$TBK1_positive, "TBK1+", "TBK1-"),
  levels = c("TBK1-", "TBK1+")
)

table(kup$Disease_group, kup$TBK1_status)

############################################################
## B.3 Subset TBK1-positive Kupffer cells
############################################################

kup_tbk1_pos <- subset(
  kup,
  subset =
    TBK1_positive == TRUE
)

DefaultAssay(kup_tbk1_pos) <- "RNA"

kup_tbk1_pos$Disease_group <- factor(
  kup_tbk1_pos$Disease_group,
  levels = disease_levels
)

cat("Number of TBK1-positive Kupffer cells:\n")
print(ncol(kup_tbk1_pos))

cat("TBK1-positive Kupffer cells by disease group:\n")

print(
  table(
    kup_tbk1_pos$Disease_group
  )
)


############################################################
## B.4 Define inflammasome-related genes
############################################################

inflam_genes <- c(
  "AIM2",
  "NLRP3",
  "NLRC4",
  "NLRP1",
  "NLRP12",
  "PYCARD",
  "CASP1",
  "CASP4",
  "CASP5",
  "GSDMD",
  "IL18",
  "IL1B"
)

inflam_genes_use <- inflam_genes[
  inflam_genes %in%
    rownames(kup_tbk1_pos)
]

missing_inflam_genes <- setdiff(
  inflam_genes,
  inflam_genes_use
)

cat("Genes used:\n")
print(inflam_genes_use)

cat("Missing genes:\n")
print(missing_inflam_genes)

if (
  length(inflam_genes_use) == 0
) {
  stop(
    paste0(
      "None of the inflammasome-related genes ",
      "were found in kup_tbk1_pos."
    )
  )
}


############################################################
## B.5 Extract expression
## in TBK1-positive Kupffer cells
############################################################

inflam_expr_tbk1_kup <- FetchData(
  kup_tbk1_pos,
  vars = c(
    inflam_genes_use,
    "Disease_group",
    "DonorID"
  )
) %>%
  filter(
    !is.na(Disease_group),
    !is.na(DonorID)
  ) %>%
  pivot_longer(
    cols = all_of(
      inflam_genes_use
    ),
    names_to = "gene",
    values_to = "expr"
  ) %>%
  mutate(
    Disease_group = factor(
      Disease_group,
      levels = disease_levels
    ),
    gene = factor(
      gene,
      levels = inflam_genes
    ),
    gene_positive =
      expr > 0
  )

write.csv(
  inflam_expr_tbk1_kup,
  "FigureB_TBK1_positive_Kupffer_inflammasome_single_cell_expression.csv",
  row.names = FALSE
)


############################################################
## B.6 Donor-level summary
############################################################

inflam_tbk1_kup_donor <-
  inflam_expr_tbk1_kup %>%
  group_by(
    DonorID,
    Disease_group,
    gene
  ) %>%
  summarise(
    mean_expr = mean(
      expr,
      na.rm = TRUE
    ),
    median_expr = median(
      expr,
      na.rm = TRUE
    ),
    frac_positive = mean(
      gene_positive,
      na.rm = TRUE
    ),
    n_cells = n(),
    .groups = "drop"
  )

write.csv(
  inflam_tbk1_kup_donor,
  "FigureB_TBK1_positive_Kupffer_inflammasome_donor_level.csv",
  row.names = FALSE
)


############################################################
## B.7 Heatmap source table
############################################################

heat_expr_raw_tbk1_kup <-
  inflam_tbk1_kup_donor %>%
  group_by(
    Disease_group,
    gene
  ) %>%
  summarise(
    mean_expr_group = mean(
      mean_expr,
      na.rm = TRUE
    ),
    median_expr_group = median(
      mean_expr,
      na.rm = TRUE
    ),
    n_donors = n_distinct(
      DonorID
    ),
    .groups = "drop"
  ) %>%
  complete(
    Disease_group = factor(
      disease_levels,
      levels = disease_levels
    ),
    gene = factor(
      inflam_genes,
      levels = inflam_genes
    ),
    fill = list(
      mean_expr_group = NA_real_,
      median_expr_group = NA_real_,
      n_donors = 0
    )
  ) %>%
  mutate(
    Disease_group = factor(
      Disease_group,
      levels = disease_levels
    ),
    gene = factor(
      gene,
      levels = rev(
        inflam_genes
      )
    )
  )

write.csv(
  heat_expr_raw_tbk1_kup,
  "FigureB_TBK1_positive_Kupffer_inflammasome_heatmap_source_table.csv",
  row.names = FALSE
)


############################################################
## B.8 Color limits
############################################################

inflam_color_limits <- c(
  0,
  quantile(
    heat_expr_raw_tbk1_kup$mean_expr_group,
    probs = 0.95,
    na.rm = TRUE
  )
)

inflam_color_limits


############################################################
## B.9 Theme
############################################################

final_theme_gene_heatmap <- theme_classic(
  base_size = 16
) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 22,
      hjust = 0.5
    ),
    
    plot.subtitle = element_text(
      size = 16,
      hjust = 0.5
    ),
    
    axis.title = element_blank(),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold",
      size = 20,
      color = "black"
    ),
    
    axis.text.y = element_text(
      face = "italic",
      size = 20,
      color = "black"
    ),
    
    legend.title = element_text(
      face = "bold",
      size = 15
    ),
    
    legend.text = element_text(
      size = 14
    ),
    
    plot.margin = margin(
      8,
      8,
      8,
      8
    )
  )


############################################################
## B.10 Plot inflammasome-gene heatmap
############################################################

p_figure_B <- ggplot(
  heat_expr_raw_tbk1_kup,
  aes(
    x = Disease_group,
    y = gene,
    fill = mean_expr_group
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.6
  ) +
  scale_fill_gradientn(
    colors = c(
      "white",
      "#FEE5D9",
      "#FCAE91",
      "#FB6A4A",
      "#CB181D"
    ),
    limits = inflam_color_limits,
    oob = scales::squish,
    na.value = "white",
    name = "Mean\nexpression"
  ) +
  labs(
    title = expression(
      paste(
        "Inflammasome-gene expression in ",
        italic("TBK1"),
        "+ Kupffer cells"
      )
    ),
    subtitle =
      "Raw donor-level mean expression",
    x = "",
    y = ""
  ) +
  final_theme_gene_heatmap

p_figure_B


############################################################
## B.11 Save Figure B
############################################################

ggsave(
  "FigureB_TBK1_positive_Kupffer_inflammasome_heatmap.pdf",
  p_figure_B,
  width = 7.8,
  height = 6.4
)

ggsave(
  "FigureB_TBK1_positive_Kupffer_inflammasome_heatmap.png",
  p_figure_B,
  width = 7.8,
  height = 6.4,
  dpi = 600
)


############################################################
## Save updated objects
############################################################

saveRDS(
  obj,
  "obj_with_Kupffer_and_Macrophage_non_Kupffer_annotation.rds"
)

saveRDS(
  macro,
  "macro_with_Kupffer_and_Macrophage_non_Kupffer_annotation.rds"
)

saveRDS(
  kup,
  "Kupffer_with_TBK1_positive_annotation.rds"
)

saveRDS(
  kup_tbk1_pos,
  "TBK1_positive_Kupffer_cells.rds"
)


############################################################
## Step5
############################################################


############################################################
## 5.1 Check input
############################################################

if (!exists("kup")) {
  stop(
    paste0(
      "Object 'kup' is not found. ",
      "Please run the Kupffer extraction first."
    )
  )
}

DefaultAssay(kup) <- "RNA"

if (
  !"state_plot" %in%
  colnames(kup@meta.data)
) {
  stop(
    paste0(
      "state_plot is not found in kup. ",
      "Please run the TBK1/CASP1 state-definition step first."
    )
  )
}

table(
  kup$state_plot,
  useNA = "ifany"
)


############################################################
## 5.2 Confirm state order
############################################################

kup$state_plot <- factor(kup$state_plot, levels = state_levels)
table(kup$state_plot)

############################################################
## 5.3 Build state summary table
############################################################

state_summary <- data.frame(
  state = kup$state_plot
) %>%
  filter(
    !is.na(state)
  ) %>%
  mutate(
    state = factor(
      state,
      levels = state_levels
    )
  ) %>%
  count(
    state,
    name = "n_cells"
  ) %>%
  complete(
    state = factor(
      state_levels,
      levels = state_levels
    ),
    fill = list(
      n_cells = 0
    )
  ) %>%
  mutate(
    total_cells = sum(
      n_cells
    ),
    fraction =
      n_cells / total_cells,
    percentage =
      fraction * 100
  )

state_summary

write.csv(
  state_summary,
  "Step5_Kupffer_TBK1_CASP1_state_overall_distribution.csv",
  row.names = FALSE
)


############################################################
## 5.4 Plot labels
############################################################

state_labels <- c(
  "Other",
  expression(italic("TBK1") * "+"),
  expression(italic("CASP1") * "+"),
  expression(italic("TBK1") * "+" * italic("CASP1") * "+")
)

## Figure-specific palette retained to reproduce the original panel
state_colors_overall <- c(
  "Other" = "grey70",
  "TBK1+" = "#E69F00",
  "CASP1+" = "#56B4E9",
  "TBK1+CASP1+" = "#D55E00"
)

############################################################
## 5.5 Original plot theme
############################################################

state_bar_theme <- theme_classic(
  base_size = 16
) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 22,
      hjust = 0.5
    ),
    
    axis.title.x = element_blank(),
    
    axis.title.y = element_text(
      face = "bold",
      size = 20,
      color = "black"
    ),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold",
      size = 18,
      color = "black"
    ),
    
    axis.text.y = element_text(
      size = 18,
      color = "black"
    ),
    
    legend.position = "none",
    
    plot.margin = margin(
      8,
      8,
      8,
      8
    )
  )


############################################################
## 5.6 Plot overall state distribution
############################################################

p_state_distribution <- ggplot(
  state_summary,
  aes(
    x = state,
    y = percentage,
    fill = state
  )
) +
  geom_col(
    width = 0.65,
    color = "black",
    linewidth = 0.7
  ) +
  scale_fill_manual(
    values = state_colors_overall,
    breaks = state_levels
  ) +
  scale_x_discrete(
    breaks = state_levels,
    labels = state_labels
  ) +
  scale_y_continuous(
    limits = c(
      0,
      max(
        state_summary$percentage,
        na.rm = TRUE
      ) * 1.18
    ),
    expand = expansion(
      mult = c(
        0,
        0.03
      )
    )
  ) +
  labs(
    title = "State Distribution",
    y = "Fraction of Kupffer cells (%)"
  ) +
  state_bar_theme

p_state_distribution


############################################################
## 5.7 Save Step5
############################################################

ggsave(
  "Step5_Kupffer_TBK1_CASP1_state_overall_distribution_barplot.pdf",
  p_state_distribution,
  width = 5.4,
  height = 5.8
)

ggsave(
  "Step5_Kupffer_TBK1_CASP1_state_overall_distribution_barplot.png",
  p_state_distribution,
  width = 5.4,
  height = 5.8,
  dpi = 600
)


############################################################
## Step6
############################################################


############################################################
## 6.1 Check input object
############################################################

if (!exists("kup")) {
  stop(
    paste0(
      "Object 'kup' is not found. ",
      "Please run the Kupffer extraction first."
    )
  )
}

DefaultAssay(kup) <- "RNA"


############################################################
## 6.3 Extract TBK1 and CASP1 expression
############################################################

coupling_df <- FetchData(
  kup,
  vars = c(
    "TBK1",
    "CASP1",
    "Disease_group",
    "DonorID"
  )
) %>%
  filter(
    !is.na(TBK1),
    !is.na(CASP1),
    !is.na(Disease_group),
    !is.na(DonorID)
  ) %>%
  mutate(
    Disease_group =
      as.character(
        Disease_group
      ),
    
    Disease_group = ifelse(
      Disease_group == "MASL",
      "MASLD",
      Disease_group
    ),
    
    Disease_group = factor(
      Disease_group,
      levels = disease_levels
    ),
    
    TBK1_positive =
      TBK1 > 0,
    
    CASP1_positive =
      CASP1 > 0
  )

table(coupling_df$Disease_group)

table(
  coupling_df$TBK1_positive,
  coupling_df$CASP1_positive
)

table(
  coupling_df$Disease_group,
  coupling_df$TBK1_positive,
  coupling_df$CASP1_positive
)

write.csv(
  coupling_df,
  "Step6_Kupffer_TBK1_CASP1_coupling_single_cell_source_table.csv",
  row.names = FALSE
)


############################################################
## 6.4 Fisher's exact test function
############################################################

get_fisher_or <- function(
    df_one_disease
) {
  
  tab <- table(
    TBK1_positive =
      df_one_disease$TBK1_positive,
    
    CASP1_positive =
      df_one_disease$CASP1_positive
  )
  
  ##########################################################
  ## Ensure a complete 2 x 2 table
  ##########################################################
  
  full_tab <- matrix(
    0,
    nrow = 2,
    ncol = 2,
    dimnames = list(
      TBK1_positive = c(
        "FALSE",
        "TRUE"
      ),
      CASP1_positive = c(
        "FALSE",
        "TRUE"
      )
    )
  )
  
  full_tab[
    rownames(tab),
    colnames(tab)
  ] <- tab
  
  fisher_res <- fisher.test(
    full_tab
  )
  
  out <- data.frame(
    n_TBK1_neg_CASP1_neg =
      full_tab[
        "FALSE",
        "FALSE"
      ],
    
    n_TBK1_neg_CASP1_pos =
      full_tab[
        "FALSE",
        "TRUE"
      ],
    
    n_TBK1_pos_CASP1_neg =
      full_tab[
        "TRUE",
        "FALSE"
      ],
    
    n_TBK1_pos_CASP1_pos =
      full_tab[
        "TRUE",
        "TRUE"
      ],
    
    odds_ratio =
      as.numeric(
        fisher_res$estimate
      ),
    
    conf_low =
      fisher_res$conf.int[1],
    
    conf_high =
      fisher_res$conf.int[2],
    
    p_value =
      fisher_res$p.value
  )
  
  return(out)
}


############################################################
## 6.5 Run Fisher's exact test
## within each disease group
############################################################

coupling_or_df <- coupling_df %>%
  group_by(
    Disease_group
  ) %>%
  group_modify(
    ~ get_fisher_or(.x)
  ) %>%
  ungroup() %>%
  mutate(
    Disease_group = factor(
      Disease_group,
      levels = disease_levels
    ),
    
    label = paste0(
      "OR=",
      sprintf(
        "%.2f",
        odds_ratio
      ),
      "\nP=",
      signif(
        p_value,
        2
      )
    )
  )

coupling_or_df

write.csv(
  coupling_or_df,
  "Step6_Kupffer_TBK1_CASP1_coupling_OR_by_disease.csv",
  row.names = FALSE
)


############################################################
## 6.6 Interaction test
############################################################

glm_input <- coupling_df %>%
  mutate(
    CASP1_positive_num =
      as.integer(
        CASP1_positive
      ),
    
    TBK1_positive_factor = factor(
      TBK1_positive,
      levels = c(
        FALSE,
        TRUE
      )
    ),
    
    Disease_group = factor(
      Disease_group,
      levels = disease_levels
    )
  )


############################################################
## Model without interaction
############################################################

glm_no_interaction <- glm(
  CASP1_positive_num ~
    TBK1_positive_factor +
    Disease_group,
  data = glm_input,
  family = binomial()
)


############################################################
## Model with interaction
############################################################

glm_interaction <- glm(
  CASP1_positive_num ~
    TBK1_positive_factor *
    Disease_group,
  data = glm_input,
  family = binomial()
)


############################################################
## Likelihood-ratio test
############################################################

interaction_test <- anova(
  glm_no_interaction,
  glm_interaction,
  test = "Chisq"
)

interaction_test

interaction_p <-
  interaction_test$`Pr(>Chi)`[2]

interaction_label <- paste0(
  "Interaction P = ",
  signif(
    interaction_p,
    2
  )
)

interaction_label

write.csv(
  data.frame(
    interaction_p =
      interaction_p
  ),
  "Step6_Kupffer_TBK1_CASP1_coupling_interaction_test.csv",
  row.names = FALSE
)


############################################################
## 6.7 Prepare plot range
############################################################

finite_ci_high <-
  coupling_or_df$conf_high[
    is.finite(
      coupling_or_df$conf_high
    )
  ]

finite_ci_low <-
  coupling_or_df$conf_low[
    is.finite(
      coupling_or_df$conf_low
    )
  ]


############################################################
## Define upper Y-axis boundary
############################################################

if (
  length(finite_ci_high) == 0
) {
  
  y_upper <- max(
    coupling_or_df$odds_ratio,
    na.rm = TRUE
  ) * 2
  
} else {
  
  y_upper <- max(
    coupling_or_df$odds_ratio,
    finite_ci_high,
    na.rm = TRUE
  ) * 1.25
}

y_upper <- max(
  y_upper,
  2.8
)


############################################################
## Define lower Y-axis boundary
############################################################

if (
  length(finite_ci_low) == 0
) {
  
  y_lower <- 0.5
  
} else {
  
  y_lower <- min(
    coupling_or_df$conf_low,
    coupling_or_df$odds_ratio,
    na.rm = TRUE
  ) * 0.9
}

y_lower <- max(
  0.4,
  min(
    y_lower,
    0.7
  )
)


############################################################
## Plotting values
############################################################

coupling_plot_df <- coupling_or_df %>%
  mutate(
    conf_low_plot = ifelse(
      is.finite(
        conf_low
      ),
      conf_low,
      y_lower
    ),
    
    conf_high_plot = ifelse(
      is.finite(
        conf_high
      ),
      conf_high,
      y_upper * 0.95
    ),
    
    label = paste0(
      "OR=",
      sprintf(
        "%.2f",
        odds_ratio
      ),
      "\nP=",
      signif(
        p_value,
        2
      )
    ),
    
    label_y = pmin(
      conf_high_plot * 1.08,
      y_upper * 0.88
    )
  )


############################################################
## Position of the interaction P-value
############################################################

interaction_label_df <- data.frame(
  x = 3.0,
  y = y_upper * 0.9,
  label = paste0(
    "Interaction P = ",
    signif(
      interaction_p,
      2
    )
  )
)


############################################################
## 6.8 Original coupling plot theme
############################################################

coupling_theme <- theme_classic(
  base_size = 16
) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 22,
      hjust = 0.5
    ),
    
    axis.title.x = element_blank(),
    
    axis.title.y = element_text(
      face = "bold",
      size = 20,
      color = "black"
    ),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold",
      size = 18,
      color = "black"
    ),
    
    axis.text.y = element_text(
      size = 18,
      color = "black"
    ),
    
    plot.margin = margin(
      8,
      8,
      8,
      8
    )
  )


############################################################
## 6.9 Plot coupling odds ratios
############################################################

p_coupling_or <- ggplot(
  coupling_plot_df,
  aes(
    x = Disease_group,
    y = odds_ratio
  )
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.8,
    color = "grey45"
  ) +
  
  ##########################################################
## Vertical confidence-interval line
##########################################################

geom_linerange(
  aes(
    ymin = conf_low_plot,
    ymax = conf_high_plot
  ),
  linewidth = 1.0,
  color = "black"
) +
  
  ##########################################################
## Horizontal confidence-interval caps
##########################################################

geom_errorbar(
  aes(
    ymin = conf_low_plot,
    ymax = conf_high_plot
  ),
  width = 0.12,
  linewidth = 1.0,
  color = "black"
) +
  
  geom_point(
    size = 4.6,
    color = "black"
  ) +
  
  geom_text(
    aes(
      y = label_y,
      label = label
    ),
    hjust = -0.05,
    vjust = 0,
    size = 4.6,
    color = "black"
  ) +
  
  geom_text(
    data =
      interaction_label_df,
    aes(
      x = x,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    size = 5.2,
    fontface = "bold",
    color = "black"
  ) +
  
  scale_y_continuous(
    limits = c(
      y_lower,
      y_upper
    ),
    expand = expansion(
      mult = c(
        0.02,
        0.08
      )
    )
  ) +
  
  labs(
    title = expression(
      paste(
        italic("TBK1"),
        "\u2013",
        italic("CASP1"),
        " coupling is preserved across disease stages"
      )
    ),
    y = "Odds ratio"
  ) +
  
  coupling_theme

p_coupling_or


############################################################
## 6.10 Save Step6
############################################################

ggsave(
  "Step6_Kupffer_TBK1_CASP1_coupling_OR_by_disease_fixed_CI.pdf",
  p_coupling_or,
  width = 7.8,
  height = 5.8
)

ggsave(
  "Step6_Kupffer_TBK1_CASP1_coupling_OR_by_disease_fixed_CI.png",
  p_coupling_or,
  width = 7.8,
  height = 5.8,
  dpi = 600
)


############################################################
## 6.11 Save current Kupffer object
############################################################

saveRDS(
  kup,
  "Step6_Kupffer_with_TBK1_CASP1_states_and_coupling.rds"
)


############################################################
## Functional rewiring heatmap
############################################################


############################################################
## Confirm input and reuse existing states
############################################################

DefaultAssay(kup) <- "RNA"
cat("Kupffer cells by disease group and state:\n")
print(table(kup$Disease_group, kup$state_plot))

############################################################
## Define pathway gene sets
############################################################

module_gene_list_v2 <- list(
  
  TBK1_associated_IFN = c(
    "IRF1",
    "IRF3",
    "IRF7",
    "STAT1",
    "STAT2",
    "IFIH1",
    "ISG15"
  ),
  
  ER_stress = c(
    "DDIT3",
    "HSPA5",
    "ATF4",
    "XBP1",
    "ERN1"
  ),
  
  Mito_stress = c(
    "HSPD1",
    "HSPE1",
    "DNAJA3",
    "LONP1"
  ),
  
  Mitophagy = c(
    "PINK1",
    "PRKN",
    "BNIP3",
    "BNIP3L",
    "OPTN",
    "SQSTM1"
  ),
  
  Autophagy_adaptor = c(
    "TAX1BP1",
    "CALCOCO2",
    "SQSTM1",
    "OPTN",
    "ATG5",
    "ATG7"
  ),
  
  Ubiquitin_scaffold = c(
    "TRAF6",
    "TAB1",
    "TAB2",
    "TAB3",
    "MAP3K7",
    "CYLD",
    "OTULIN",
    "RNF31",
    "SHARPIN",
    "RBCK1"
  ),
  
  Canonical_NFkB = c(
    "NFKB1",
    "RELB",
    "NFKBIA",
    "NFKBIZ",
    "TNFAIP3",
    "BIRC3"
  ),
  
  Inflammasome_program = c(
    "CASP4",
    "CASP5",
    "GSDMD",
    "IL18",
    "IL1B",
    "PYCARD",
    "NLRP3",
    "AIM2"
  )
)


############################################################
##Keep only genes present in the dataset
############################################################

module_gene_list_v2 <- lapply(
  module_gene_list_v2,
  function(x) {
    intersect(
      x,
      rownames(kup)
    )
  }
)

module_gene_list_v2 <- module_gene_list_v2[
  lengths(module_gene_list_v2) > 0
]

cat("Genes used in each pathway module:\n")
print(module_gene_list_v2)

write.csv(
  data.frame(
    pathway = rep(
      names(module_gene_list_v2),
      lengths(module_gene_list_v2)
    ),
    gene = unlist(
      module_gene_list_v2
    )
  ),
  "Functional_rewiring_pathway_genes_used_TBK1_CASP1_removed.csv",
  row.names = FALSE
)


############################################################
##Calculate module scores
## Same prefix and logic as the original non-double figure
############################################################


for (nm in names(module_gene_list_v2)) {
  
  kup <- AddModuleScore(
    object = kup,
    features = list(
      module_gene_list_v2[[nm]]
    ),
    name = paste0(
      nm,
      "_v2_"
    )
  )
}


############################################################
##Define module-score columns
############################################################

module_cols_v2 <- c(
  TBK1_associated_IFN =
    "TBK1_associated_IFN_v2_1",
  
  ER_stress =
    "ER_stress_v2_1",
  
  Mito_stress =
    "Mito_stress_v2_1",
  
  Mitophagy =
    "Mitophagy_v2_1",
  
  Autophagy_adaptor =
    "Autophagy_adaptor_v2_1",
  
  Ubiquitin_scaffold =
    "Ubiquitin_scaffold_v2_1",
  
  Canonical_NFkB =
    "Canonical_NFkB_v2_1",
  
  Inflammasome_program =
    "Inflammasome_program_v2_1"
)

module_cols_v2 <- module_cols_v2[
  module_cols_v2 %in%
    colnames(kup@meta.data)
]

cat("Module-score columns used:\n")
print(module_cols_v2)

if (length(module_cols_v2) == 0) {
  stop(
    "No module-score columns were generated."
  )
}


############################################################
##Pathway order and labels
############################################################

pathway_order_v2 <- c(
  "Ubiquitin_scaffold",
  "Autophagy_adaptor",
  "Mitophagy",
  "Mito_stress",
  "ER_stress",
  "Canonical_NFkB",
  "TBK1_associated_IFN",
  "Inflammasome_program"
)

pathway_order_v2 <- pathway_order_v2[
  pathway_order_v2 %in%
    names(module_cols_v2)
]

pathway_label_parse <- c(
  Ubiquitin_scaffold =
    "'Ubiquitin scaffold'",
  
  Autophagy_adaptor =
    "'Autophagy adaptor'",
  
  Mitophagy =
    "'Mitophagy'",
  
  Mito_stress =
    "'Mito stress'",
  
  ER_stress =
    "'ER stress'",
  
  Canonical_NFkB =
    "'Canonical NF-'*kappa*'B'",
  
  TBK1_associated_IFN =
    "italic('TBK1')*'-associated IFN'",
  
  Inflammasome_program =
    "'Inflammasome-associated'"
)


############################################################
##Theme for upper heatmap
############################################################

pathway_heatmap_theme <- theme_classic(
  base_size = 16
) +
  theme(
    plot.title = element_text(
      face = "plain",
      size = 22,
      hjust = 0.5,
      color = "black"
    ),
    
    plot.subtitle = element_text(
      face = "plain",
      size = 16,
      hjust = 0.5,
      color = "black"
    ),
    
    axis.title = element_blank(),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold",
      size = 20,
      color = "black"
    ),
    
    axis.text.y = element_text(
      face = "plain",
      size = 20,
      color = "black"
    ),
    
    legend.title = element_text(
      face = "bold",
      size = 15,
      color = "black"
    ),
    
    legend.text = element_text(
      size = 14,
      color = "black"
    ),
    
    legend.key.height = unit(
      0.65,
      "cm"
    ),
    
    legend.key.width = unit(
      0.35,
      "cm"
    ),
    
    plot.margin = margin(
      8,
      8,
      8,
      8
    )
  )


############################################################
## TBK1+CASP1+ versus all non-double Kupffer cells
############################################################


############################################################
##Define Double and Non_double
############################################################

kup$double_status <- ifelse(
  kup$state_plot == "TBK1+CASP1+",
  "Double",
  "Non_double"
)

kup$double_status <- factor(
  kup$double_status,
  levels = c(
    "Non_double",
    "Double"
  )
)

cat("Double versus non-double counts:\n")

print(
  table(
    kup$Disease_group,
    kup$double_status
  )
)


############################################################
##Generate long dataframe
## Single-cell-level module scores
############################################################

heat_long_v2 <- kup@meta.data %>%
  mutate(
    Disease_group = factor(
      Disease_group,
      levels = disease_levels
    ),
    
    double_status = ifelse(
      state_plot == "TBK1+CASP1+",
      "Double",
      "Non_double"
    )
  ) %>%
  filter(
    !is.na(Disease_group),
    !is.na(double_status)
  ) %>%
  select(
    Disease_group,
    DonorID,
    double_status,
    all_of(
      unname(module_cols_v2)
    )
  ) %>%
  pivot_longer(
    cols = all_of(
      unname(module_cols_v2)
    ),
    names_to = "module_col",
    values_to = "score"
  ) %>%
  mutate(
    pathway = names(module_cols_v2)[
      match(
        module_col,
        unname(module_cols_v2)
      )
    ]
  )

write.csv(
  heat_long_v2,
  "Functional_rewiring_DP_vs_non_double_single_cell_module_scores.csv",
  row.names = FALSE
)


############################################################
##Calculate Double minus Non_double
############################################################

heat_stats_v2 <- heat_long_v2 %>%
  group_by(
    Disease_group,
    pathway
  ) %>%
  summarise(
    delta =
      mean(
        score[
          double_status == "Double"
        ],
        na.rm = TRUE
      ) -
      mean(
        score[
          double_status == "Non_double"
        ],
        na.rm = TRUE
      ),
    
    n_cell_double = sum(
      double_status == "Double"
    ),
    
    n_cell_non_double = sum(
      double_status == "Non_double"
    ),
    
    p = tryCatch(
      wilcox.test(
        score[
          double_status == "Double"
        ],
        score[
          double_status == "Non_double"
        ]
      )$p.value,
      error = function(e) NA_real_
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(
      p,
      method = "BH"
    ),
    
    sig = case_when(
      is.na(p_adj) ~ "",
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ ""
    )
  )


############################################################
##Complete missing combinations
############################################################

heat_stats_v2 <- heat_stats_v2 %>%
  complete(
    Disease_group = factor(
      disease_levels,
      levels = disease_levels
    ),
    
    pathway = pathway_order_v2,
    
    fill = list(
      delta = NA_real_,
      n_cell_double = 0,
      n_cell_non_double = 0,
      p = NA_real_,
      p_adj = NA_real_,
      sig = ""
    )
  ) %>%
  mutate(
    Disease_group = factor(
      Disease_group,
      levels = disease_levels
    ),
    
    pathway = factor(
      pathway,
      levels = rev(
        pathway_order_v2
      )
    )
  )

print(heat_stats_v2)

write.csv(
  heat_stats_v2,
  "Functional_rewiring_DP_vs_non_double_heatmap_stats.csv",
  row.names = FALSE
)


############################################################
##Plot upper heatmap
############################################################

p_functional_rewiring_non_double <- ggplot(
  heat_stats_v2,
  aes(
    x = Disease_group,
    y = pathway,
    fill = delta
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.6
  ) +
  scale_y_discrete(
    labels = function(x) {
      parse(
        text = pathway_label_parse[x]
      )
    }
  ) +
  scale_fill_gradient2(
    low = "#3B6FB6",
    mid = "white",
    high = "#C41E1E",
    midpoint = 0,
    limits = c(
      -0.15,
      0.15
    ),
    oob = scales::squish,
    name = "Double minus\nnon-double"
  ) +
  labs(
    title = "Functional rewiring",
    
    subtitle = expression(
      paste(
        italic("TBK1"),
        "+",
        italic("CASP1"),
        "+ Kupffer cells versus non-double Kupffer cells"
      )
    ),
    
    x = "",
    y = ""
  ) +
  pathway_heatmap_theme

p_functional_rewiring_non_double


############################################################
##Save upper heatmap
############################################################

ggsave(
  "Functional_rewiring_DP_vs_non_double_heatmap.pdf",
  p_functional_rewiring_non_double,
  width = 8.8,
  height = 6.4
)

ggsave(
  "Functional_rewiring_DP_vs_non_double_heatmap.png",
  p_functional_rewiring_non_double,
  width = 8.8,
  height = 6.4,
  dpi = 600
)


############################################################
## DotPlot:
## MASLD Kupffer TBK1/CASP1 states
############################################################


DefaultAssay(kup) <- "RNA"


############################################################
## 1. Extract MASLD Kupffer cells
############################################################

kup_masld <- subset(
  kup,
  subset = Disease_group == "MASLD"
)

DefaultAssay(kup_masld) <- "RNA"

kup_masld$state_plot <- factor(
  kup_masld$state_plot,
  levels = c(
    "Other",
    "TBK1+",
    "CASP1+",
    "TBK1+CASP1+"
  )
)

table(kup_masld$state_plot)


############################################################
## 2. Gene panels
############################################################

features_tbk_masld <- list(
  
  "TBK1-associated IFN" = c(
    "IFIH1",
    "IRF1",
    "IRF3",
    "IRF7",
    "STAT1",
    "STAT2",
    "ISG15"
  ),
  
  "Inflammasome-associated" = c(
    "CASP4",
    "CASP5",
    "GSDMD",
    "IL18",
    "IL1B",
    "PYCARD",
    "NLRP3",
    "AIM2"
  )
)


############################################################
## 3. Keep genes present in the object
############################################################

features_tbk_masld <- lapply(
  features_tbk_masld,
  function(x) {
    x[x %in% rownames(kup_masld)]
  }
)

print(features_tbk_masld)


############################################################
## 4. Original theme
############################################################

dot_theme_masld <- theme_classic(
  base_size = 16
) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 22,
      hjust = 0.5
    ),
    
    axis.title = element_blank(),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 17,
      color = "black",
      face = "italic"
    ),
    
    axis.text.y = element_text(
      size = 18,
      color = "black"
    ),
    
    legend.title = element_text(
      face = "bold",
      size = 15
    ),
    
    legend.text = element_text(
      size = 14
    ),
    
    plot.margin = margin(
      10,
      15,
      35,
      10
    )
  )


############################################################
## 5. DotPlot B
## TBK1-associated IFN
############################################################

ifn_genes_use <- features_tbk_masld[[
  "TBK1-associated IFN"
]]

p_ifn_dot_masld <- DotPlot(
  kup_masld,
  features = ifn_genes_use,
  group.by = "state_plot",
  cols = c(
    "#D9D9D9",
    "#3B00FF"
  ),
  dot.scale = 6
) +
  scale_y_discrete(
    limits = rev(
      c(
        "Other",
        "TBK1+",
        "CASP1+",
        "TBK1+CASP1+"
      )
    ),
    
    labels = c(
      expression(
        italic("TBK1") * "+" *
          italic("CASP1") * "+"
      ),
      
      expression(
        italic("CASP1") * "+"
      ),
      
      expression(
        italic("TBK1") * "+"
      ),
      
      "Other"
    )
  ) +
  labs(
    title = expression(
      bolditalic("TBK1") *
        bold("-associated IFN")
    ),
    x = "",
    y = ""
  ) +
  dot_theme_masld

p_ifn_dot_masld


############################################################
## 6. DotPlot C
## Inflammasome-associated
############################################################

inflam_genes_use <- features_tbk_masld[[
  "Inflammasome-associated"
]]

p_inflam_dot_masld <- DotPlot(
  kup_masld,
  features = inflam_genes_use,
  group.by = "state_plot",
  cols = c(
    "#D9D9D9",
    "#3B00FF"
  ),
  dot.scale = 6
) +
  scale_y_discrete(
    limits = rev(
      c(
        "Other",
        "TBK1+",
        "CASP1+",
        "TBK1+CASP1+"
      )
    ),
    
    labels = c(
      expression(
        italic("TBK1") * "+" *
          italic("CASP1") * "+"
      ),
      
      expression(
        italic("CASP1") * "+"
      ),
      
      expression(
        italic("TBK1") * "+"
      ),
      
      "Other"
    )
  ) +
  labs(
    title = "Inflammasome-associated",
    x = "",
    y = ""
  ) +
  dot_theme_masld

p_inflam_dot_masld


############################################################
## 7. Combine
############################################################

p_dot_masld_combined <-
  p_ifn_dot_masld +
  p_inflam_dot_masld +
  plot_layout(
    ncol = 2,
    widths = c(
      1,
      1
    )
  )

p_dot_masld_combined


############################################################
## 8. Save
############################################################

ggsave(
  "MASLD_TBK1_IFN_and_inflammasome_dotplots.pdf",
  p_dot_masld_combined,
  width = 14,
  height = 7.2
)

ggsave(
  "MASLD_TBK1_IFN_and_inflammasome_dotplots.png",
  p_dot_masld_combined,
  width = 14,
  height = 7.2,
  dpi = 600
)


############################################################
## Heatmap
############################################################


############################################################
## 1. Define CASP4 positivity only
############################################################

if (!"CASP4" %in% rownames(kup)) {
  stop("CASP4 was not found in the Seurat object.")
}

casp4_expr <- FetchData(kup, vars = "CASP4")
kup$CASP4_expr <- casp4_expr$CASP4
kup$CASP4_pos <- kup$CASP4_expr > 0

############################################################
## 2. Subset TBK1+CASP1+ Kupffer cells
############################################################

kup_double <- subset(
  kup,
  subset = state_plot == "TBK1+CASP1+"
)

DefaultAssay(kup_double) <- "RNA"

kup_double$Disease_group <- factor(
  kup_double$Disease_group,
  levels = disease_levels
)


############################################################
## 3. Define CASP4 status within double-positive cells
############################################################

kup_double$CASP4_status <- ifelse(
  kup_double$CASP4_pos,
  "CASP4+",
  "CASP4-"
)

kup_double$CASP4_status <- factor(
  kup_double$CASP4_status,
  levels = c(
    "CASP4-",
    "CASP4+"
  )
)

table(
  kup_double$Disease_group,
  kup_double$CASP4_status
)


############################################################
## 4. Define IFN genes
## Keep original genes and order
############################################################

ifn_genes <- c(
  "IRF1",
  "IRF3",
  "IRF7",
  "STAT1",
  "STAT2",
  "IFIH1",
  "ISG15"
)

ifn_genes_use <- ifn_genes[
  ifn_genes %in%
    rownames(kup_double)
]

if (length(ifn_genes_use) == 0) {
  stop(
    "None of the IFN genes were found in the Seurat object."
  )
}

cat("IFN genes used:\n")
print(ifn_genes_use)

cat("Missing IFN genes:\n")
print(
  setdiff(
    ifn_genes,
    ifn_genes_use
  )
)


############################################################
## 5. Extract IFN-gene expression
############################################################

expr_ifn <- FetchData(
  kup_double,
  vars = c(
    ifn_genes_use,
    "Disease_group",
    "CASP4_status"
  )
)

expr_ifn_long <- expr_ifn %>%
  pivot_longer(
    cols = all_of(
      ifn_genes_use
    ),
    names_to = "gene",
    values_to = "expression"
  ) %>%
  filter(
    !is.na(Disease_group),
    !is.na(CASP4_status)
  )

write.csv(
  expr_ifn_long,
  "TBK1_IFN_genes_CASP4_single_cell_expression_source_table.csv",
  row.names = FALSE
)


############################################################
## 6. Calculate mean expression difference
##
## delta = mean(CASP4+) - mean(CASP4-)
############################################################

heat_df <- expr_ifn_long %>%
  group_by(
    Disease_group,
    gene
  ) %>%
  summarise(
    mean_CASP4_pos = mean(
      expression[
        CASP4_status == "CASP4+"
      ],
      na.rm = TRUE
    ),
    
    mean_CASP4_neg = mean(
      expression[
        CASP4_status == "CASP4-"
      ],
      na.rm = TRUE
    ),
    
    n_CASP4_pos = sum(
      CASP4_status == "CASP4+"
    ),
    
    n_CASP4_neg = sum(
      CASP4_status == "CASP4-"
    ),
    
    delta =
      mean_CASP4_pos -
      mean_CASP4_neg,
    
    p = tryCatch(
      wilcox.test(
        expression[
          CASP4_status == "CASP4+"
        ],
        expression[
          CASP4_status == "CASP4-"
        ]
      )$p.value,
      error = function(e) {
        NA_real_
      }
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(
      p,
      method = "BH"
    ),
    
    sig = case_when(
      is.na(p_adj) ~ "",
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ ""
    )
  )


############################################################
## 7. Complete missing disease/gene combinations
############################################################

heat_df <- heat_df %>%
  complete(
    Disease_group = factor(
      disease_levels,
      levels = disease_levels
    ),
    
    gene = ifn_genes_use,
    
    fill = list(
      mean_CASP4_pos = NA_real_,
      mean_CASP4_neg = NA_real_,
      n_CASP4_pos = 0,
      n_CASP4_neg = 0,
      delta = NA_real_,
      p = NA_real_,
      p_adj = NA_real_,
      sig = ""
    )
  )

heat_df$Disease_group <- factor(
  heat_df$Disease_group,
  levels = disease_levels
)

heat_df$gene <- factor(
  heat_df$gene,
  levels = rev(
    ifn_genes_use
  )
)

print(heat_df)


############################################################
## 8. Original theme
############################################################

final_theme <- theme_classic(
  base_size = 12
) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 13,
      hjust = 0
    ),
    
    plot.subtitle = element_text(
      size = 9,
      hjust = 0.5
    ),
    
    axis.title = element_text(
      face = "bold",
      size = 11
    ),
    
    axis.text = element_text(
      color = "black",
      size = 9
    ),
    
    legend.title = element_text(
      face = "bold",
      size = 9
    ),
    
    legend.text = element_text(
      size = 8
    ),
    
    strip.text = element_text(
      face = "bold",
      size = 10
    ),
    
    plot.margin = margin(
      4,
      4,
      4,
      4
    )
  )


############################################################
## 9. Plot heatmap
############################################################

p_ifn_casp4_heat <- ggplot(
  heat_df,
  aes(
    x = Disease_group,
    y = gene,
    fill = delta
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.55
  ) +
  scale_fill_gradient2(
    low = "#3B6FB6",
    mid = "white",
    high = "#C41E1E",
    midpoint = 0,
    
    limits = c(
      -0.25,
      0.25
    ),
    
    oob = scales::squish,
    
    name = "CASP4+ minus\nCASP4-"
  ) +
  labs(
    title =
      expression(
        paste(
          italic("TBK1"),
          "-associated IFN genes in ",
          italic("TBK1"),
          "+",
          italic("CASP1"),
          "+ Kupffer cells"
        )
      ),
    
    subtitle =
      expression(
        paste(
          italic("TBK1"),
          "+",
          italic("CASP1"),
          "+",
          italic("CASP4"),
          "+ versus ",
          italic("TBK1"),
          "+",
          italic("CASP1"),
          "+",
          italic("CASP4"),
          "-"
        )
      ),
    
    x = "",
    y = ""
  ) +
  final_theme +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 13,
      hjust = 0
    ),
    
    plot.subtitle = element_text(
      face = "bold",
      size = 8.5,
      hjust = 0.5
    ),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold"
    ),
    
    axis.text.y = element_text(
      face = "italic"
    ),
    
    legend.position = "right"
  )

p_ifn_casp4_heat


############################################################
## 10. Save figure and source table
############################################################

ggsave(
  "TBK1_IFN_genes_CASP4_delta_heatmap_TBK1_CASP1_Kupffer.pdf",
  p_ifn_casp4_heat,
  width = 5.4,
  height = 4.1
)

ggsave(
  "TBK1_IFN_genes_CASP4_delta_heatmap_TBK1_CASP1_Kupffer.png",
  p_ifn_casp4_heat,
  width = 5.4,
  height = 4.1,
  dpi = 600
)

write.csv(
  heat_df,
  "TBK1_IFN_genes_CASP4_delta_heatmap_source_table.csv",
  row.names = FALSE
)


############################################################
## 11. Save double-positive Kupffer object
############################################################

saveRDS(
  kup_double,
  "TBK1_CASP1_double_positive_Kupffer_with_CASP4_status.rds"
)


############################################################
## AIM2, NLRP3 and CASP4 positivity
############################################################


DefaultAssay(kup) <- "RNA"


############################################################
## 1. Check the existing columns
############################################################

if (
  !"state_plot" %in%
  colnames(kup@meta.data)
) {
  stop(
    paste0(
      "state_plot is not found in kup. ",
      "Please run the previous TBK1/CASP1 state-definition code first."
    )
  )
}

if (
  !"Disease_group" %in%
  colnames(kup@meta.data)
) {
  stop(
    "Disease_group is not found in kup."
  )
}

if (
  !"DonorID" %in%
  colnames(kup@meta.data)
) {
  stop(
    "DonorID is not found in kup."
  )
}


############################################################
## 2. Confirm the existing disease and state order
############################################################

kup$Disease_group <- factor(kup$Disease_group, levels = disease_levels)
kup$state_plot <- factor(kup$state_plot, levels = state_levels)
print(table(kup$Disease_group, kup$state_plot))

############################################################
## 3. Genes to analyse
############################################################

genes_use <- c(
  "AIM2",
  "NLRP3",
  "CASP4"
)

genes_use <- genes_use[
  genes_use %in%
    rownames(kup)
]

cat("Genes used:\n")
print(genes_use)

cat("Missing genes:\n")
print(
  setdiff(
    c(
      "AIM2",
      "NLRP3",
      "CASP4"
    ),
    genes_use
  )
)

if (
  length(genes_use) == 0
) {
  stop(
    "None of AIM2, NLRP3 or CASP4 were found in kup."
  )
}


############################################################
## 4. Extract single-cell expression
############################################################

expr_state_gene <- FetchData(
  kup,
  vars = c(
    genes_use,
    "Disease_group",
    "state_plot",
    "DonorID"
  )
) %>%
  filter(
    !is.na(Disease_group),
    !is.na(state_plot),
    !is.na(DonorID)
  ) %>%
  pivot_longer(
    cols = all_of(
      genes_use
    ),
    names_to = "gene",
    values_to = "expr"
  ) %>%
  mutate(
    gene_pos = expr > 0,
    
    Disease_group = factor(
      Disease_group,
      levels = disease_levels
    ),
    
    state_plot = factor(
      state_plot,
      levels = state_levels
    )
  )

write.csv(
  expr_state_gene,
  "AIM2_NLRP3_CASP4_by_TBK1_CASP1_state_single_cell_source_data.csv",
  row.names = FALSE
)


############################################################
## 5. Donor-level summary
############################################################

expr_state_gene_donor <- expr_state_gene %>%
  group_by(
    DonorID,
    Disease_group,
    state_plot,
    gene
  ) %>%
  summarise(
    mean_expr = mean(
      expr,
      na.rm = TRUE
    ),
    
    frac_positive = mean(
      gene_pos,
      na.rm = TRUE
    ),
    
    n_cells = n(),
    
    .groups = "drop"
  )

write.csv(
  expr_state_gene_donor,
  "AIM2_NLRP3_CASP4_by_TBK1_CASP1_state_donor_values.csv",
  row.names = FALSE
)

print(
  expr_state_gene_donor %>%
    count(
      gene,
      state_plot,
      Disease_group,
      name = "n_donors"
    )
)


############################################################
## 6. Kruskal-Wallis tests
############################################################

stats_frac <- expr_state_gene_donor %>%
  group_by(
    gene,
    state_plot
  ) %>%
  summarise(
    n_disease_groups = n_distinct(
      Disease_group[
        !is.na(frac_positive)
      ]
    ),
    
    n_donors = sum(
      !is.na(frac_positive)
    ),
    
    p_kruskal = tryCatch(
      {
        dat_test <- data.frame(
          frac_positive = frac_positive,
          Disease_group = Disease_group
        ) %>%
          filter(
            !is.na(frac_positive),
            !is.na(Disease_group)
          )
        
        if (
          n_distinct(
            dat_test$Disease_group
          ) >= 2
        ) {
          kruskal.test(
            frac_positive ~ Disease_group,
            data = dat_test
          )$p.value
        } else {
          NA_real_
        }
      },
      error = function(e) {
        NA_real_
      }
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(
      p_kruskal,
      method = "BH"
    ),
    
    sig = case_when(
      is.na(p_adj) ~ "",
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

print(stats_frac)

write.csv(
  stats_frac,
  "AIM2_NLRP3_CASP4_positive_fraction_Kruskal_Wallis_BH_stats.csv",
  row.names = FALSE
)


############################################################
## 7.
############################################################

significant_stats_frac <- stats_frac %>%
  filter(
    !is.na(p_adj),
    p_adj < 0.05
  )

cat("Comparisons significant after BH correction:\n")
print(significant_stats_frac)


############################################################
## 8. Labels for italic gene names and states
############################################################

gene_labeller <- as_labeller(
  c(
    "AIM2" =
      "italic('AIM2')",
    
    "CASP4" =
      "italic('CASP4')",
    
    "NLRP3" =
      "italic('NLRP3')"
  ),
  default = label_parsed
)

state_labeller <- as_labeller(
  c(
    "Other" =
      "'Other'",
    
    "TBK1+" =
      "italic('TBK1')*'+'",
    
    "CASP1+" =
      "italic('CASP1')*'+'",
    
    "TBK1+CASP1+" =
      "italic('TBK1')*'+'*italic('CASP1')*'+'"
  ),
  default = label_parsed
)


############################################################
## 9. Original theme
############################################################

theme_step19 <- theme_classic(
  base_size = 16
) +
  theme(
    legend.position = "none",
    
    plot.title = element_text(
      face = "plain",
      size = 18,
      hjust = 0.5,
      color = "black"
    ),
    
    axis.title.x = element_blank(),
    
    axis.title.y = element_text(
      face = "plain",
      size = 16,
      color = "black"
    ),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold",
      size = 12,
      color = "black"
    ),
    
    axis.text.y = element_text(
      size = 12,
      color = "black"
    ),
    
    strip.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.8
    ),
    
    strip.text.x = element_text(
      face = "bold",
      size = 11,
      color = "black"
    ),
    
    strip.text.y = element_text(
      face = "bold",
      size = 11,
      color = "black"
    ),
    
    panel.spacing.x = unit(
      0.35,
      "lines"
    ),
    
    panel.spacing.y = unit(
      0.55,
      "lines"
    ),
    
    plot.margin = margin(
      8,
      8,
      8,
      8
    )
  )


############################################################
## 10. Plot positive fraction per donor
############################################################

p_frac_state_gene <- ggplot(
  expr_state_gene_donor,
  aes(
    x = Disease_group,
    y = frac_positive,
    fill = Disease_group
  )
) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA,
    alpha = 0.85,
    color = "black",
    linewidth = 0.55
  ) +
  geom_jitter(
    width = 0.12,
    size = 1.8,
    alpha = 0.85,
    color = "black"
  ) +
  facet_grid(
    rows = vars(
      gene
    ),
    cols = vars(
      state_plot
    ),
    scales = "free_y",
    labeller = labeller(
      gene = gene_labeller,
      state_plot = state_labeller
    )
  ) +
  scale_y_continuous(
    labels = scales::percent_format(
      accuracy = 1
    )
  ) +
  scale_fill_manual(
    values = disease_cols
  ) +
  labs(
    title = expression(
      paste(
        italic("AIM2"),
        ", ",
        italic("NLRP3"),
        " and ",
        italic("CASP4"),
        " positivity across ",
        italic("TBK1"),
        "/",
        italic("CASP1"),
        " Kupffer-cell states"
      )
    ),
    x = NULL,
    y = "Positive fraction per donor"
  ) +
  theme_step19

p_frac_state_gene


############################################################
## 11. Save figure
############################################################

ggsave(
  "AIM2_NLRP3_CASP4_positive_fraction_by_TBK1_CASP1_state_italic.pdf",
  p_frac_state_gene,
  width = 12,
  height = 8
)

ggsave(
  "AIM2_NLRP3_CASP4_positive_fraction_by_TBK1_CASP1_state_italic.png",
  p_frac_state_gene,
  width = 12,
  height = 8,
  dpi = 600
)


############################################################
## TBK1/CASP1 state composition across disease stages
############################################################


DefaultAssay(kup) <- "RNA"


############################################################
## 1. Check required existing metadata
############################################################

required_meta <- c(
  "state_plot",
  "Disease_group",
  "DonorID"
)

missing_meta <- setdiff(
  required_meta,
  colnames(kup@meta.data)
)

if (length(missing_meta) > 0) {
  stop(
    paste(
      "The following metadata columns are missing:",
      paste(
        missing_meta,
        collapse = ", "
      )
    )
  )
}


############################################################
## 2. Confirm existing disease/state order and use global colors
############################################################

kup$Disease_group <- factor(kup$Disease_group, levels = disease_levels)
kup$state_plot <- factor(kup$state_plot, levels = state_levels)

state_colors <- c(
  "Other" = "grey70",
  "TBK1+" = "#D9A042",
  "CASP1+" = "#66A3D2",
  "TBK1+CASP1+" = "#C5332D"
)

############################################################
## 3. Build pooled cell-level state-composition table
############################################################

state_disease_df <- data.frame(
  Disease_group = kup$Disease_group,
  state = kup$state_plot
) %>%
  filter(
    !is.na(Disease_group),
    !is.na(state)
  ) %>%
  mutate(
    Disease_group = factor(
      Disease_group,
      levels = disease_levels
    ),
    
    state = factor(
      state,
      levels = state_levels
    )
  ) %>%
  count(
    Disease_group,
    state,
    name = "n_cells"
  ) %>%
  complete(
    Disease_group = factor(
      disease_levels,
      levels = disease_levels
    ),
    
    state = factor(
      state_levels,
      levels = state_levels
    ),
    
    fill = list(
      n_cells = 0
    )
  ) %>%
  group_by(
    Disease_group
  ) %>%
  mutate(
    total_cells = sum(
      n_cells
    ),
    
    fraction =
      n_cells / total_cells,
    
    percentage =
      fraction * 100,
    
    percentage_label = sprintf(
      "%.1f%%",
      percentage
    )
  ) %>%
  ungroup()

print(state_disease_df)

write.csv(
  state_disease_df,
  "TBK1_CASP1_state_composition_by_disease_with_percentages.csv",
  row.names = FALSE
)


############################################################
## 4. Prepare cell-level binary data for statistics
############################################################

double_positive_stats_df <- kup@meta.data %>%
  transmute(
    DonorID = factor(
      DonorID
    ),
    
    Disease_group = factor(
      as.character(Disease_group),
      levels = disease_levels
    ),
    
    state = factor(
      as.character(state_plot),
      levels = state_levels
    ),
    
    double_positive = as.integer(
      state == "TBK1+CASP1+"
    )
  ) %>%
  filter(
    !is.na(DonorID),
    !is.na(Disease_group),
    !is.na(state)
  )

table(
  double_positive_stats_df$Disease_group,
  double_positive_stats_df$double_positive
)

write.csv(
  double_positive_stats_df,
  "TBK1_CASP1_double_positive_cell_level_statistics_input.csv",
  row.names = FALSE
)


############################################################
## 5. Four-stage mixed-effects logistic regression
############################################################

model_stage_null <- glmer(
  double_positive ~
    1 +
    (1 | DonorID),
  
  data = double_positive_stats_df,
  
  family = binomial(),
  
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(
      maxfun = 2e5
    )
  )
)

model_stage_full <- glmer(
  double_positive ~
    Disease_group +
    (1 | DonorID),
  
  data = double_positive_stats_df,
  
  family = binomial(),
  
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(
      maxfun = 2e5
    )
  )
)


############################################################
## 6 Overall likelihood-ratio test
############################################################

overall_stage_test <- anova(
  model_stage_null,
  model_stage_full,
  test = "Chisq"
)

print(overall_stage_test)

overall_stage_p <-
  overall_stage_test$`Pr(>Chisq)`[2]

overall_stage_label <- paste0(
  "Disease-stage effect P = ",
  format.pval(
    overall_stage_p,
    digits = 2,
    eps = 0.001
  )
)

print(overall_stage_label)

write.csv(
  data.frame(
    test =
      "Overall four-stage disease effect",
    
    p_value =
      overall_stage_p
  ),
  "TBK1_CASP1_double_positive_overall_disease_stage_test.csv",
  row.names = FALSE
)


############################################################
## 7. Pairwise comparisons among four disease stages
############################################################

stage_emmeans <- emmeans(
  model_stage_full,
  specs = ~ Disease_group,
  type = "response"
)

stage_pairwise <- contrast(
  stage_emmeans,
  method = "pairwise",
  adjust = "BH"
)

stage_pairwise_df <- as.data.frame(
  stage_pairwise
)

print(stage_pairwise_df)

write.csv(
  stage_pairwise_df,
  "TBK1_CASP1_double_positive_pairwise_disease_comparisons_BH.csv",
  row.names = FALSE
)


############################################################
## 8. Estimated probability of double positivity
############################################################

stage_probability_df <- as.data.frame(
  stage_emmeans
)

print(stage_probability_df)

write.csv(
  stage_probability_df,
  "TBK1_CASP1_double_positive_estimated_probabilities_by_disease.csv",
  row.names = FALSE
)


############################################################
## 9. Exploratory control versus all disease groups combined
############################################################

double_positive_binary_df <-
  double_positive_stats_df %>%
  mutate(
    Control_vs_disease = ifelse(
      Disease_group == "control",
      "control",
      "disease"
    ),
    
    Control_vs_disease = factor(
      Control_vs_disease,
      levels = c(
        "control",
        "disease"
      )
    )
  )

table(
  double_positive_binary_df$Control_vs_disease,
  double_positive_binary_df$double_positive
)


############################################################
## 9.1 Null and full models
############################################################

model_control_disease_null <- glmer(
  double_positive ~
    1 +
    (1 | DonorID),
  
  data = double_positive_binary_df,
  
  family = binomial(),
  
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(
      maxfun = 2e5
    )
  )
)

model_control_disease_full <- glmer(
  double_positive ~
    Control_vs_disease +
    (1 | DonorID),
  
  data = double_positive_binary_df,
  
  family = binomial(),
  
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(
      maxfun = 2e5
    )
  )
)


############################################################
## 9.2 Likelihood-ratio test
############################################################

control_vs_disease_test <- anova(
  model_control_disease_null,
  model_control_disease_full,
  test = "Chisq"
)

print(control_vs_disease_test)

control_vs_disease_p <-
  control_vs_disease_test$`Pr(>Chisq)`[2]


############################################################
## 9.3 Extract odds ratio and confidence interval
############################################################

control_vs_disease_coef <- summary(
  model_control_disease_full
)$coefficients

control_vs_disease_log_or <-
  control_vs_disease_coef[
    "Control_vs_diseasedisease",
    "Estimate"
  ]

control_vs_disease_se <-
  control_vs_disease_coef[
    "Control_vs_diseasedisease",
    "Std. Error"
  ]

control_vs_disease_or <- exp(
  control_vs_disease_log_or
)

control_vs_disease_ci_low <- exp(
  control_vs_disease_log_or -
    1.96 * control_vs_disease_se
)

control_vs_disease_ci_high <- exp(
  control_vs_disease_log_or +
    1.96 * control_vs_disease_se
)

control_vs_disease_result <- data.frame(
  comparison =
    "All disease groups versus control",
  
  odds_ratio =
    control_vs_disease_or,
  
  conf_low =
    control_vs_disease_ci_low,
  
  conf_high =
    control_vs_disease_ci_high,
  
  p_value =
    control_vs_disease_p
)

print(control_vs_disease_result)

write.csv(
  control_vs_disease_result,
  "TBK1_CASP1_double_positive_control_vs_all_disease_exploratory.csv",
  row.names = FALSE
)


############################################################
## 9.4 Estimated probabilities
############################################################

control_vs_disease_emmeans <- emmeans(
  model_control_disease_full,
  specs = ~ Control_vs_disease,
  type = "response"
)

control_vs_disease_probability_df <-
  as.data.frame(
    control_vs_disease_emmeans
  )

print(control_vs_disease_probability_df)

write.csv(
  control_vs_disease_probability_df,
  "TBK1_CASP1_double_positive_control_vs_disease_probabilities.csv",
  row.names = FALSE
)


############################################################
## 10. Original stacked-bar theme
############################################################

stacked_theme <- theme_classic(
  base_size = 16
) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 22,
      hjust = 0.5
    ),
    
    axis.title.x = element_blank(),
    
    axis.title.y = element_text(
      face = "bold",
      size = 20,
      color = "black"
    ),
    
    axis.text.x = element_text(
      face = "bold",
      size = 18,
      color = "black"
    ),
    
    axis.text.y = element_text(
      size = 18,
      color = "black"
    ),
    
    legend.title = element_blank(),
    
    legend.text = element_text(
      size = 16,
      color = "black"
    ),
    
    plot.margin = margin(
      8,
      8,
      8,
      8
    )
  )


############################################################
## 11. Stacked bar plot with percentage labels
############################################################

p_state_composition_disease <- ggplot(
  state_disease_df,
  aes(
    x = Disease_group,
    y = percentage,
    fill = state
  )
) +
  geom_col(
    width = 0.75,
    color = "black",
    linewidth = 0.4
  ) +
  
  ##########################################################
## Percentage labels inside each segment
##########################################################

geom_text(
  aes(
    label = percentage_label
  ),
  position = position_stack(
    vjust = 0.5
  ),
  size = 4.5,
  color = "black"
) +
  
  scale_fill_manual(
    values = state_colors,
    breaks = state_levels,
    labels = c(
      "Other",
      
      expression(
        italic("TBK1") * "+"
      ),
      
      expression(
        italic("CASP1") * "+"
      ),
      
      expression(
        italic("TBK1") * "+" *
          italic("CASP1") * "+"
      )
    )
  ) +
  
  scale_y_continuous(
    breaks = c(
      0,
      25,
      50,
      75,
      100
    ),
    
    labels = function(x) {
      paste0(
        x,
        "%"
      )
    },
    
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +
  
  coord_cartesian(
    ylim = c(
      0,
      115
    ),
    clip = "off"
  ) +
  
##########################################################
## Overall disease-stage P value
##########################################################

annotate(
  geom = "text",
  x = 2.5,
  y = 109,
  label = overall_stage_label,
  size = 4.5,
  color = "black"
) +
  
  labs(
    title = expression(
      paste(
        italic("TBK1"),
        "/",
        italic("CASP1"),
        " state composition across disease stages"
      )
    ),
    
    y =
      "Fraction of Kupffer cells"
  ) +
  
  stacked_theme

p_state_composition_disease


############################################################
## 12. Save figure
############################################################

ggsave(
  "TBK1_CASP1_state_composition_across_disease_stages_with_percentages_and_statistics.pdf",
  p_state_composition_disease,
  width = 8.2,
  height = 6.1
)

ggsave(
  "TBK1_CASP1_state_composition_across_disease_stages_with_percentages_and_statistics.png",
  p_state_composition_disease,
  width = 8.2,
  height = 6.1,
  dpi = 600
)


############################################################
## Final save
############################################################

saveRDS(
  obj,
  "MASLD_snRNA_obj_final_annotated.rds"
)

saveRDS(
  macro,
  "MASLD_snRNA_macrophages_final.rds"
)

saveRDS(
  kup,
  "MASLD_snRNA_Kupffer_final_with_all_metadata.rds"
)
