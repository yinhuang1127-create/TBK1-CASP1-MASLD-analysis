
############################################################
## 1. Check objects inherited from the main script
############################################################

required_objects <- c(
  "spatial_all",
  "plot_final",
  "sensor_expr",
  "out_dir"
)

missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    logical(1),
    inherits = TRUE
  )
]

if (length(missing_objects) > 0) {
  stop(
    "Missing objects from the main script: ",
    paste(
      missing_objects,
      collapse = ", "
    ),
    ". Run the main spatial-transcriptomics script first, ",
    "then run this sensitivity-analysis script in the same R session."
  )
}


############################################################
## 2. Output directory
############################################################

sensitivity_out_dir <- file.path(
  out_dir,
  "Kupffer_threshold_sensitivity"
)

dir.create(
  sensitivity_out_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


############################################################
## 3. Check columns already defined in the main script
############################################################

required_meta_cols <- c(
  "cell",
  "sample",
  "condition",
  "Kupffer1",
  "tbk_casp_state"
)

missing_meta_cols <- setdiff(
  required_meta_cols,
  colnames(plot_final)
)

if (length(missing_meta_cols) > 0) {
  stop(
    "Missing required columns in plot_final: ",
    paste(
      missing_meta_cols,
      collapse = ", "
    )
  )
}

required_sensor_cols <- c(
  "Aim2",
  "Casp4",
  "Nlrp3"
)

missing_sensor_cols <- setdiff(
  required_sensor_cols,
  colnames(sensor_expr)
)

if (length(missing_sensor_cols) > 0) {
  stop(
    "Missing required columns in sensor_expr: ",
    paste(
      missing_sensor_cols,
      collapse = ", "
    )
  )
}


############################################################
## 4. Build common dataframe
##
## Reuse the Tbk1/Casp1 state defined in plot_final.
## Sensor positivity is defined once here.
############################################################

base_df <- plot_final %>%
  select(
    all_of(
      required_meta_cols
    )
  ) %>%
  left_join(
    sensor_expr %>%
      rownames_to_column(
        "cell"
      ),
    by = "cell"
  ) %>%
  mutate(
    condition = factor(
      condition,
      levels = c(
        "Chow",
        "CDAA-HFD"
      )
    ),

    tbk_casp_state = factor(
      as.character(
        tbk_casp_state
      ),
      levels = c(
        "Tbk1-Casp1-",
        "Tbk1-Casp1+",
        "Tbk1+Casp1-",
        "Tbk1+Casp1+"
      )
    ),

    Aim2_pos =
      Aim2 > 0,

    Casp4_pos =
      Casp4 > 0,

    Nlrp3_pos =
      Nlrp3 > 0
  )

cat(
  "\nTbk1/Casp1 states inherited from the main script:\n"
)

print(
  table(
    base_df$tbk_casp_state,
    useNA = "ifany"
  )
)


############################################################
## 5. Threshold settings
############################################################

cutoff_levels <- c(
  "Top 5%",
  "Top 10%",
  "Top 15%"
)

threshold_settings <- tibble(
  cutoff = cutoff_levels,

  quantile_prob = c(
    0.95,
    0.90,
    0.85
  )
)


############################################################
## 6. Safe Fisher exact test
############################################################

safe_fisher_or <- function(
    df,
    group_col,
    outcome_col,
    ref_group,
    target_group
) {

  df_sub <- df %>%
    filter(
      .data[[group_col]] %in%
        c(
          ref_group,
          target_group
        )
    )

  n_ref <- sum(
    df_sub[[group_col]] ==
      ref_group,
    na.rm = TRUE
  )

  n_target <- sum(
    df_sub[[group_col]] ==
      target_group,
    na.rm = TRUE
  )

  empty_result <- function() {
    tibble(
      OR = NA_real_,
      CI_low = NA_real_,
      CI_high = NA_real_,
      p_value = NA_real_,
      n_ref = n_ref,
      n_target = n_target
    )
  }

  if (
    n_ref == 0 ||
      n_target == 0
  ) {
    return(
      empty_result()
    )
  }

  tab <- table(
    factor(
      df_sub[[group_col]],
      levels = c(
        ref_group,
        target_group
      )
    ),

    factor(
      df_sub[[outcome_col]],
      levels = c(
        FALSE,
        TRUE
      )
    )
  )

  if (!all(dim(tab) == c(2, 2))) {
    return(
      empty_result()
    )
  }

  if (
    sum(tab[, "FALSE"]) == 0 ||
      sum(tab[, "TRUE"]) == 0
  ) {
    return(
      empty_result()
    )
  }

  ft <- tryCatch(
    fisher.test(
      tab
    ),
    error = function(e) {
      NULL
    }
  )

  if (is.null(ft)) {
    return(
      empty_result()
    )
  }

  tibble(
    OR = as.numeric(
      ft$estimate
    ),

    CI_low = as.numeric(
      ft$conf.int[1]
    ),

    CI_high = as.numeric(
      ft$conf.int[2]
    ),

    p_value =
      ft$p.value,

    n_ref =
      n_ref,

    n_target =
      n_target
  )
}


############################################################
## 7. Safe Wilcoxon rank-sum test
############################################################

safe_wilcox <- function(
    df,
    value_col,
    group_col = "condition"
) {

  df_test <- df %>%
    filter(
      !is.na(
        .data[[value_col]]
      ),
      !is.na(
        .data[[group_col]]
      )
    )

  if (
    nrow(df_test) == 0 ||
      n_distinct(
        df_test[[group_col]]
      ) < 2
  ) {
    return(
      tibble(
        statistic = NA_real_,
        p_value = NA_real_
      )
    )
  }

  wt <- tryCatch(
    suppressWarnings(
      wilcox.test(
        df_test[[value_col]] ~
          df_test[[group_col]],
        exact = TRUE
      )
    ),
    error = function(e) {
      NULL
    }
  )

  if (is.null(wt)) {
    wt <- tryCatch(
      suppressWarnings(
        wilcox.test(
          df_test[[value_col]] ~
            df_test[[group_col]],
          exact = FALSE
        )
      ),
      error = function(e) {
        NULL
      }
    )
  }

  if (is.null(wt)) {
    return(
      tibble(
        statistic = NA_real_,
        p_value = NA_real_
      )
    )
  }

  tibble(
    statistic = as.numeric(
      wt$statistic
    ),

    p_value =
      wt$p.value
  )
}


############################################################
## 8. Run one threshold
############################################################

run_threshold_analysis <- function(
    base_df,
    cutoff_label,
    quantile_prob
) {

  ##########################################################
  ## A. Recalculate only the Kupffer-high threshold
  ## independently within each sample
  ##########################################################

  df_thr <- base_df %>%
    group_by(
      sample
    ) %>%
    mutate(
      kup_cut = quantile(
        Kupffer1,
        probs = quantile_prob,
        na.rm = TRUE,
        names = FALSE
      ),

      Kupffer_region_sens = ifelse(
        Kupffer1 >= kup_cut,
        "Kupffer-high",
        "Other"
      )
    ) %>%
    ungroup()


  ##########################################################
  ## B. Supplementary Table 2 source
  ## Sample-level double-positive representation
  ##########################################################

  double_positive_per_sample <- df_thr %>%
    filter(
      Kupffer_region_sens ==
        "Kupffer-high"
    ) %>%
    group_by(
      condition,
      sample
    ) %>%
    summarise(
      n_double_positive = sum(
        tbk_casp_state ==
          "Tbk1+Casp1+",
        na.rm = TRUE
      ),

      n_kupffer_high =
        n(),

      frac_double_in_kupffer =
        n_double_positive /
        n_kupffer_high,

      .groups =
        "drop"
    ) %>%
    mutate(
      cutoff =
        cutoff_label
    )


  ##########################################################
  ## Chow versus CDAA-HFD test
  ##########################################################

  wilcox_result <- safe_wilcox(
    double_positive_per_sample,

    value_col =
      "frac_double_in_kupffer",

    group_col =
      "condition"
  ) %>%
    mutate(
      cutoff =
        cutoff_label
    )


  ##########################################################
  ## C. Supplementary Table 3 source
  ## Sensor enrichment:
  ## double-positive versus double-negative
  ##########################################################

  sensor_compare_df <- df_thr %>%
    filter(
      condition ==
        "CDAA-HFD",

      Kupffer_region_sens ==
        "Kupffer-high",

      tbk_casp_state %in%
        c(
          "Tbk1-Casp1-",
          "Tbk1+Casp1+"
        )
    )


  sensor_results <- bind_rows(

    safe_fisher_or(
      df =
        sensor_compare_df,

      group_col =
        "tbk_casp_state",

      outcome_col =
        "Aim2_pos",

      ref_group =
        "Tbk1-Casp1-",

      target_group =
        "Tbk1+Casp1+"
    ) %>%
      mutate(
        gene =
          "Aim2"
      ),

    safe_fisher_or(
      df =
        sensor_compare_df,

      group_col =
        "tbk_casp_state",

      outcome_col =
        "Casp4_pos",

      ref_group =
        "Tbk1-Casp1-",

      target_group =
        "Tbk1+Casp1+"
    ) %>%
      mutate(
        gene =
          "Casp4"
      ),

    safe_fisher_or(
      df =
        sensor_compare_df,

      group_col =
        "tbk_casp_state",

      outcome_col =
        "Nlrp3_pos",

      ref_group =
        "Tbk1-Casp1-",

      target_group =
        "Tbk1+Casp1+"
    ) %>%
      mutate(
        gene =
          "Nlrp3"
      )
  ) %>%
    mutate(
      cutoff =
        cutoff_label,

      p_adj_BH = p.adjust(
        p_value,
        method = "BH"
      )
    )


  ##########################################################
  ## D. Supplementary Table 4 source
  ## Casp4 coupling across Tbk1/Casp1 states
  ##########################################################

  casp4_source_df <- df_thr %>%
    filter(
      condition ==
        "CDAA-HFD",

      Kupffer_region_sens ==
        "Kupffer-high"
    )


  casp4_coupling_results <- bind_rows(

    safe_fisher_or(
      df =
        casp4_source_df,

      group_col =
        "tbk_casp_state",

      outcome_col =
        "Casp4_pos",

      ref_group =
        "Tbk1-Casp1-",

      target_group =
        "Tbk1-Casp1+"
    ) %>%
      mutate(
        comparison =
          "Tbk1-Casp1+ vs Tbk1-Casp1-"
      ),

    safe_fisher_or(
      df =
        casp4_source_df,

      group_col =
        "tbk_casp_state",

      outcome_col =
        "Casp4_pos",

      ref_group =
        "Tbk1-Casp1-",

      target_group =
        "Tbk1+Casp1-"
    ) %>%
      mutate(
        comparison =
          "Tbk1+Casp1- vs Tbk1-Casp1-"
      ),

    safe_fisher_or(
      df =
        casp4_source_df,

      group_col =
        "tbk_casp_state",

      outcome_col =
        "Casp4_pos",

      ref_group =
        "Tbk1-Casp1-",

      target_group =
        "Tbk1+Casp1+"
    ) %>%
      mutate(
        comparison =
          "Tbk1+Casp1+ vs Tbk1-Casp1-"
      )
  ) %>%
    mutate(
      cutoff =
        cutoff_label,

      p_adj_BH = p.adjust(
        p_value,
        method = "BH"
      )
    )


  ##########################################################
  ## Return objects needed for Supplementary Tables 2–4
  ##########################################################

  list(
    double_positive_per_sample =
      double_positive_per_sample,

    wilcox_result =
      wilcox_result,

    sensor_results =
      sensor_results,

    casp4_coupling_results =
      casp4_coupling_results
  )
}


############################################################
## 9. Run all thresholds
############################################################

threshold_results <- lapply(
  seq_len(
    nrow(
      threshold_settings
    )
  ),

  function(i) {
    run_threshold_analysis(
      base_df =
        base_df,

      cutoff_label =
        threshold_settings$cutoff[i],

      quantile_prob =
        threshold_settings$quantile_prob[i]
    )
  }
)

names(threshold_results) <-
  threshold_settings$cutoff


############################################################
## 10. Combine threshold results
############################################################

all_double_positive_per_sample <- bind_rows(
  lapply(
    threshold_results,
    `[[`,
    "double_positive_per_sample"
  )
)

all_wilcox <- bind_rows(
  lapply(
    threshold_results,
    `[[`,
    "wilcox_result"
  )
)

all_sensor_results <- bind_rows(
  lapply(
    threshold_results,
    `[[`,
    "sensor_results"
  )
)

all_casp4_coupling <- bind_rows(
  lapply(
    threshold_results,
    `[[`,
    "casp4_coupling_results"
  )
)


############################################################
## 11. Supplementary Table 2
## Representation of Tbk1+Casp1+ spots
############################################################

double_positive_summary <- all_double_positive_per_sample %>%
  group_by(
    cutoff,
    condition
  ) %>%
  summarise(
    n_samples =
      n(),

    total_double_positive =
      sum(
        n_double_positive
      ),

    mean_per_sample =
      mean(
        n_double_positive
      ),

    median_per_sample =
      median(
        n_double_positive
      ),

    Q1 = quantile(
      n_double_positive,
      0.25
    ),

    Q3 = quantile(
      n_double_positive,
      0.75
    ),

    min_per_sample =
      min(
        n_double_positive
      ),

    max_per_sample =
      max(
        n_double_positive
      ),

    n_samples_with_zero =
      sum(
        n_double_positive == 0
      ),

    .groups =
      "drop"
  )


supp_table2 <- double_positive_summary %>%
  left_join(
    all_wilcox %>%
      select(
        cutoff,

        `Between-condition P` =
          p_value
      ),

    by =
      "cutoff"
  ) %>%
  mutate(
    cutoff = factor(
      cutoff,
      levels = cutoff_levels
    ),

    IQR = paste0(
      signif(
        Q1,
        3
      ),
      "–",
      signif(
        Q3,
        3
      )
    ),

    Range = paste0(
      min_per_sample,
      "–",
      max_per_sample
    ),

    samples_with_zero = paste0(
      n_samples_with_zero,
      "/",
      n_samples
    )
  ) %>%
  arrange(
    cutoff,
    condition
  ) %>%
  transmute(
    Cutoff =
      as.character(
        cutoff
      ),

    Condition =
      as.character(
        condition
      ),

    `Total double-positive spots` =
      total_double_positive,

    `Mean spots/sample` =
      round(
        mean_per_sample,
        2
      ),

    `Median spots/sample` =
      round(
        median_per_sample,
        2
      ),

    IQR =
      IQR,

    Range =
      Range,

    `Samples with zero DP spots` =
      samples_with_zero,

    `Between-condition P` =
      signif(
        `Between-condition P`,
        4
      )
  )


############################################################
## 12. Supplementary Table 3
## Sensor-enrichment sensitivity
############################################################

supp_table3 <- all_sensor_results %>%
  mutate(
    cutoff = factor(
      cutoff,
      levels = cutoff_levels
    ),

    gene = factor(
      gene,
      levels = c(
        "Aim2",
        "Casp4",
        "Nlrp3"
      )
    ),

    Sensor = recode(
      as.character(
        gene
      ),

      "Aim2" =
        "Aim2",

      "Casp4" =
        "Casp4 (Casp11)",

      "Nlrp3" =
        "Nlrp3"
    ),

    `95% CI` = paste0(
      signif(
        CI_low,
        3
      ),
      "–",
      signif(
        CI_high,
        3
      )
    )
  ) %>%
  arrange(
    cutoff,
    gene
  ) %>%
  transmute(
    Cutoff =
      as.character(
        cutoff
      ),

    Sensor =
      Sensor,

    OR =
      round(
        OR,
        3
      ),

    `95% CI` =
      `95% CI`,

    `Raw P` =
      signif(
        p_value,
        4
      ),

    `BH-adjusted P` =
      signif(
        p_adj_BH,
        4
      ),

    `Reference spots` =
      n_ref,

    `Target spots` =
      n_target
  )


############################################################
## 13. Supplementary Table 4
## Casp4-coupling sensitivity
############################################################

supp_table4 <- all_casp4_coupling %>%
  mutate(
    cutoff = factor(
      cutoff,
      levels = cutoff_levels
    ),

    comparison_order = case_when(
      comparison ==
        "Tbk1-Casp1+ vs Tbk1-Casp1-" ~ 1,

      comparison ==
        "Tbk1+Casp1- vs Tbk1-Casp1-" ~ 2,

      comparison ==
        "Tbk1+Casp1+ vs Tbk1-Casp1-" ~ 3,

      TRUE ~
        99
    ),

    Comparison = recode(
      comparison,

      "Tbk1-Casp1+ vs Tbk1-Casp1-" =
        "Tbk1−Casp1+ vs Tbk1−Casp1−",

      "Tbk1+Casp1- vs Tbk1-Casp1-" =
        "Tbk1+Casp1− vs Tbk1−Casp1−",

      "Tbk1+Casp1+ vs Tbk1-Casp1-" =
        "Tbk1+Casp1+ vs Tbk1−Casp1−"
    ),

    `95% CI` = paste0(
      signif(
        CI_low,
        3
      ),
      "–",
      signif(
        CI_high,
        3
      )
    )
  ) %>%
  arrange(
    cutoff,
    comparison_order
  ) %>%
  transmute(
    Cutoff =
      as.character(
        cutoff
      ),

    Comparison =
      Comparison,

    OR =
      round(
        OR,
        3
      ),

    `95% CI` =
      `95% CI`,

    `Raw P` =
      signif(
        p_value,
        4
      ),

    `BH-adjusted P` =
      signif(
        p_adj_BH,
        4
      ),

    `Reference spots` =
      n_ref,

    `Target spots` =
      n_target
  )


############################################################
## 14. Print final tables
############################################################

cat(
  "\nSupplementary Table 2:\n"
)

print(
  supp_table2,
  n = Inf,
  width = Inf
)

cat(
  "\nSupplementary Table 3:\n"
)

print(
  supp_table3,
  n = Inf,
  width = Inf
)

cat(
  "\nSupplementary Table 4:\n"
)

print(
  supp_table4,
  n = Inf,
  width = Inf
)


############################################################
## 15. Save Supplementary Tables 2–4
############################################################

write.csv(
  supp_table2,
  file.path(
    sensitivity_out_dir,
    "Supplementary_Table_2_threshold_representation.csv"
  ),
  row.names = FALSE
)

write.csv(
  supp_table3,
  file.path(
    sensitivity_out_dir,
    "Supplementary_Table_3_sensor_enrichment_sensitivity.csv"
  ),
  row.names = FALSE
)

write.csv(
  supp_table4,
  file.path(
    sensitivity_out_dir,
    "Supplementary_Table_4_Casp4_coupling_sensitivity.csv"
  ),
  row.names = FALSE
)


############################################################
## 16. Save session information
############################################################

writeLines(
  capture.output(
    sessionInfo()
  ),
  con = file.path(
    sensitivity_out_dir,
    "sessionInfo_threshold_sensitivity.txt"
  )
)


############################################################
## 17. Final message
############################################################

cat(
  "\nThreshold sensitivity analysis completed.\n"
)

cat(
  "Outputs saved to:\n",
  sensitivity_out_dir,
  "\n"
)
