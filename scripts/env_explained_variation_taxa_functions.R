# -*- coding: UTF-8 -*-

# ============================================================
# Fig_env_explained_taxa_functions.R
#
# Environmental effects on microbial taxa and functions
# Domains: Bacteria, Archaea, Fungi
#
# This script generates a figure similar to:
#   top bar: random forest explained variation
#   heatmap: Spearman correlation between environmental variables and features
#   circle size: random forest variable importance
#
# All comments and plot labels are ASCII-only to avoid encoding problems.
# ============================================================

# ============================================================
# Input file 1: environmental metadata table
# ============================================================
#
# Default file:
#   env_metadata.tsv
#
# Required columns:
#   sample_id
#
# Other numeric columns are treated as environmental variables.
#
# Example:
# sample_id        TN     SOC    AP     pH     CEC    Fe     S
# 201704_MF1_0002 1.25   12.3   8.5    6.82   10.1   2.3    0.18
# 201704_MF1_0608 1.10   10.8   7.2    6.75   9.8    2.1    0.16
#
#
# ============================================================
# Input file 2: taxon abundance table
# ============================================================
#
# Default file:
#   all.marker_rpm.genus.long.tsv
#
# Required columns:
#   sample_id
#   domain
#   rank
#   lineage
#   marker_rpm
#
# Example:
# sample_id        domain    rank   lineage                                      marker_rpm
# 201704_MF1_0002 Bacteria  genus  d__Bacteria;p__Proteobacteria;g__TaxonA       12.3
# 201704_MF1_0002 Archaea   genus  d__Archaea;p__Asgardarchaeota;g__TaxonB        0.35
# 201704_MF1_0002 Fungi     genus  k__Fungi;p__Ascomycota;g__Trichophyton         0.14
#
# Notes:
#   lineage is used as the unique feature ID.
#   marker_rpm is used as taxon abundance.
#
#
# ============================================================
# Optional input file 3: functional abundance table
# ============================================================
#
# Default file:
#   function_abundance.tsv
#
# Required columns:
#   sample_id
#   domain
#   function
#   abundance
#
# Example:
# sample_id        domain    function                 abundance
# 201704_MF1_0002 Bacteria  Carbon_metabolism          12.3
# 201704_MF1_0002 Archaea   Methane_metabolism          1.5
# 201704_MF1_0002 Fungi     CAZymes                     5.2
#
# If this file does not exist, only taxon panels will be plotted.
#
# ============================================================

options(stringsAsFactors = FALSE)

# ----------------------------
# 0. Packages
# ----------------------------

packages <- c(
  "tidyverse",
  "data.table",
  "randomForest",
  "patchwork",
  "scales"
)

for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

library(tidyverse)
library(data.table)
library(randomForest)
library(patchwork)
library(scales)

# ----------------------------
# 1. Parameters
# ----------------------------

env_file <- "env_metadata.tsv"
taxa_file <- "all.marker_rpm.genus.long.tsv"
function_file <- "function_abundance.tsv"

output_dir <- "env_explained_taxa_functions"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

domain_keep <- c("Bacteria", "Archaea", "Fungi")

target_rank <- "genus"

# Use NULL to automatically use all numeric environmental variables.
env_vars <- NULL

# Example:
# env_vars <- c("TN", "SOC", "AP", "pH", "CEC", "Fe", "S")

top_n_taxa_per_domain <- 25
top_n_functions_per_domain <- 25

min_detected_samples <- 5
min_samples_for_model <- 15

cor_method <- "spearman"
p_adjust_method <- "BH"

rf_ntree <- 1000
rf_seed <- 123

log_transform_abundance <- TRUE
pseudo_count <- 1e-8

show_correlation_text <- FALSE
show_significance_stars <- TRUE

figure_width_mm <- 240
figure_height_taxa_only_mm <- 120
figure_height_taxa_function_mm <- 220

# ----------------------------
# 2. Theme and helper functions
# ----------------------------

theme_nature <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(color = "black"),
      axis.text = element_text(color = "black", size = base_size),
      axis.title = element_text(color = "black", size = base_size + 1),
      axis.line = element_line(linewidth = 0.35, color = "black"),
      axis.ticks = element_line(linewidth = 0.35, color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size + 1),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 1),
      legend.key.size = unit(0.35, "cm"),
      plot.title = element_text(face = "bold", size = base_size + 1),
      plot.margin = margin(4, 4, 4, 4)
    )
}

p_to_star <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
}

safe_cor_test <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  n <- length(x)
  
  if (n < 4 || length(unique(x)) < 2 || length(unique(y)) < 2) {
    return(tibble(r = NA_real_, p = NA_real_, n = n))
  }
  
  ct <- suppressWarnings(cor.test(x, y, method = method, exact = FALSE))
  
  tibble(
    r = unname(ct$estimate),
    p = ct$p.value,
    n = n
  )
}

safe_log_transform <- function(x) {
  if (log_transform_abundance) {
    log10(x + pseudo_count)
  } else {
    x
  }
}

rank_prefix <- c(
  domain = "d__",
  phylum = "p__",
  class = "c__",
  order = "o__",
  family = "f__",
  genus = "g__",
  species = "s__"
)

extract_rank_label <- function(lineage, rank) {
  prefix <- rank_prefix[[rank]]
  
  if (is.null(prefix)) {
    return(lineage)
  }
  
  parts <- stringr::str_split(lineage, ";", simplify = FALSE)[[1]]
  hit <- parts[stringr::str_starts(parts, fixed(prefix))]
  
  if (length(hit) == 0) {
    out <- tail(parts, 1)
  } else {
    out <- hit[1]
  }
  
  out <- stringr::str_remove(out, paste0("^", prefix))
  out <- stringr::str_replace_all(out, "_", " ")
  out
}

make_unique_labels <- function(df, feature_col = "feature_id", label_col = "feature_label") {
  df %>%
    group_by(domain, .data[[label_col]]) %>%
    mutate(
      label_n = n(),
      sh_id = stringr::str_extract(.data[[feature_col]], "SH[0-9]+\\.[A-Za-z0-9]+"),
      feature_plot_label = case_when(
        label_n == 1 ~ .data[[label_col]],
        !is.na(sh_id) ~ paste0(.data[[label_col]], " | ", sh_id),
        TRUE ~ paste0(.data[[label_col]], " | ", row_number())
      )
    ) %>%
    ungroup() %>%
    select(-label_n, -sh_id)
}

make_safe_feature_names <- function(x) {
  make.unique(gsub("[^A-Za-z0-9_.]+", "_", x))
}

# ----------------------------
# 3. Read environmental table
# ----------------------------

if (!file.exists(env_file)) {
  stop("Environmental metadata table not found: ", env_file)
}

env_raw <- data.table::fread(
  env_file,
  sep = "\t",
  header = TRUE,
  data.table = FALSE,
  check.names = FALSE
)

names(env_raw) <- trimws(names(env_raw))

if (!"sample_id" %in% names(env_raw)) {
  stop("Environmental table must contain sample_id column.")
}

env_raw <- env_raw %>%
  mutate(sample_id = as.character(sample_id))

if (is.null(env_vars)) {
  env_vars <- env_raw %>%
    select(-sample_id) %>%
    select(where(is.numeric)) %>%
    names()
}

if (length(env_vars) == 0) {
  stop("No numeric environmental variables found in env_file.")
}

missing_env_vars <- setdiff(env_vars, names(env_raw))
if (length(missing_env_vars) > 0) {
  stop("Missing environmental variables: ", paste(missing_env_vars, collapse = ", "))
}

env_dat <- env_raw %>%
  select(sample_id, all_of(env_vars)) %>%
  mutate(across(all_of(env_vars), as.numeric))

# ----------------------------
# 4. Read taxon abundance table
# ----------------------------

if (!file.exists(taxa_file)) {
  stop("Taxon abundance table not found: ", taxa_file)
}

taxa_raw <- data.table::fread(
  taxa_file,
  sep = "\t",
  header = TRUE,
  data.table = FALSE,
  check.names = FALSE
)

names(taxa_raw) <- trimws(names(taxa_raw))

required_taxa_cols <- c("sample_id", "domain", "rank", "lineage", "marker_rpm")
missing_taxa_cols <- setdiff(required_taxa_cols, names(taxa_raw))

if (length(missing_taxa_cols) > 0) {
  stop("Missing columns in taxa_file: ", paste(missing_taxa_cols, collapse = ", "))
}

taxa_long <- taxa_raw %>%
  mutate(
    sample_id = as.character(sample_id),
    domain = as.character(domain),
    rank = as.character(rank),
    lineage = as.character(lineage),
    marker_rpm = as.numeric(marker_rpm)
  ) %>%
  filter(
    domain %in% domain_keep,
    rank == target_rank,
    !is.na(marker_rpm),
    marker_rpm > 0
  ) %>%
  mutate(
    feature_type = "Taxa",
    feature_id = lineage,
    feature_label = map_chr(lineage, extract_rank_label, rank = target_rank),
    abundance = marker_rpm
  ) %>%
  select(sample_id, domain, feature_type, feature_id, feature_label, abundance)

# ----------------------------
# 5. Read optional function abundance table
# ----------------------------

function_long <- tibble()

if (file.exists(function_file)) {
  func_raw <- data.table::fread(
    function_file,
    sep = "\t",
    header = TRUE,
    data.table = FALSE,
    check.names = FALSE
  )
  
  names(func_raw) <- trimws(names(func_raw))
  
  required_func_cols <- c("sample_id", "domain", "function", "abundance")
  missing_func_cols <- setdiff(required_func_cols, names(func_raw))
  
  if (length(missing_func_cols) > 0) {
    stop("Missing columns in function_file: ", paste(missing_func_cols, collapse = ", "))
  }
  
  function_long <- func_raw %>%
    mutate(
      sample_id = as.character(sample_id),
      domain = as.character(domain),
      `function` = as.character(`function`),
      abundance = as.numeric(abundance)
    ) %>%
    filter(
      domain %in% domain_keep,
      !is.na(abundance),
      abundance > 0
    ) %>%
    transmute(
      sample_id,
      domain,
      feature_type = "Function",
      feature_id = `function`,
      feature_label = `function`,
      abundance
    )
} else {
  message("[INFO] function_file not found. Only taxon panels will be plotted: ", function_file)
}

feature_long <- bind_rows(taxa_long, function_long)

common_samples <- intersect(env_dat$sample_id, feature_long$sample_id)

if (length(common_samples) < min_samples_for_model) {
  stop("Too few shared samples between env_file and abundance tables.")
}

env_dat <- env_dat %>%
  filter(sample_id %in% common_samples) %>%
  arrange(match(sample_id, common_samples))

feature_long <- feature_long %>%
  filter(sample_id %in% common_samples)

# ----------------------------
# 6. Select abundant features
# ----------------------------

feature_info <- feature_long %>%
  group_by(feature_type, domain, feature_id, feature_label) %>%
  summarise(
    total_abundance = sum(abundance, na.rm = TRUE),
    mean_abundance = mean(abundance, na.rm = TRUE),
    detected_samples = n_distinct(sample_id[abundance > 0]),
    .groups = "drop"
  ) %>%
  filter(detected_samples >= min_detected_samples) %>%
  group_by(feature_type, domain) %>%
  arrange(desc(total_abundance), .by_group = TRUE) %>%
  mutate(
    keep_n = ifelse(feature_type == "Taxa", top_n_taxa_per_domain, top_n_functions_per_domain)
  ) %>%
  filter(row_number() <= keep_n) %>%
  ungroup() %>%
  make_unique_labels(feature_col = "feature_id", label_col = "feature_label")

feature_info <- feature_info %>%
  mutate(feature_safe_id = make_safe_feature_names(feature_id))

feature_long <- feature_long %>%
  inner_join(
    feature_info %>%
      select(feature_type, domain, feature_id, feature_label, feature_safe_id, feature_plot_label),
    by = c("feature_type", "domain", "feature_id", "feature_label")
  )

write_tsv(
  feature_info,
  file.path(output_dir, "selected_features.tsv")
)

# ----------------------------
# 7. Build abundance matrix
# ----------------------------

make_feature_matrix <- function(feature_type_name, domain_name) {
  info_sub <- feature_info %>%
    filter(feature_type == feature_type_name, domain == domain_name)
  
  if (nrow(info_sub) == 0) {
    return(NULL)
  }
  
  features <- info_sub$feature_safe_id
  
  mat <- feature_long %>%
    filter(feature_type == feature_type_name, domain == domain_name) %>%
    group_by(sample_id, feature_safe_id) %>%
    summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
    complete(
      sample_id = common_samples,
      feature_safe_id = features,
      fill = list(abundance = 0)
    ) %>%
    pivot_wider(
      names_from = feature_safe_id,
      values_from = abundance,
      values_fill = 0
    ) %>%
    arrange(match(sample_id, common_samples))
  
  mat
}

# ----------------------------
# 8. Correlation and random forest analysis
# ----------------------------

analyze_one_set <- function(feature_type_name, domain_name) {
  mat_df <- make_feature_matrix(feature_type_name, domain_name)
  
  if (is.null(mat_df)) {
    return(tibble())
  }
  
  sample_ids <- mat_df$sample_id
  
  env_sub <- env_dat %>%
    filter(sample_id %in% sample_ids) %>%
    arrange(match(sample_id, sample_ids))
  
  feature_cols <- setdiff(names(mat_df), "sample_id")
  
  results <- list()
  idx <- 1
  
  for (feature_safe in feature_cols) {
    y_raw <- mat_df[[feature_safe]]
    y <- safe_log_transform(y_raw)
    
    ok_y <- is.finite(y)
    
    if (sum(ok_y) < min_samples_for_model || length(unique(y[ok_y])) < 2) {
      next
    }
    
    rf_dat <- env_sub %>%
      select(all_of(env_vars)) %>%
      mutate(response = y) %>%
      drop_na()
    
    if (nrow(rf_dat) < min_samples_for_model || length(unique(rf_dat$response)) < 2) {
      next
    }
    
    set.seed(rf_seed)
    
    rf_fit <- randomForest(
      x = rf_dat %>% select(all_of(env_vars)),
      y = rf_dat$response,
      ntree = rf_ntree,
      importance = TRUE
    )
    
    rf_importance <- randomForest::importance(rf_fit, type = 1) %>%
      as.data.frame() %>%
      rownames_to_column("env_var")
    
    if (!"%IncMSE" %in% names(rf_importance)) {
      value_col <- setdiff(names(rf_importance), "env_var")[1]
      rf_importance <- rf_importance %>%
        rename(`%IncMSE` = all_of(value_col))
    }
    
    explained_variation <- tail(rf_fit$rsq, 1) * 100
    explained_variation <- max(0, explained_variation, na.rm = TRUE)
    
    for (ev in env_vars) {
      cor_one <- safe_cor_test(
        x = env_sub[[ev]],
        y = y,
        method = cor_method
      )
      
      imp_one <- rf_importance %>%
        filter(env_var == ev) %>%
        pull(`%IncMSE`)
      
      if (length(imp_one) == 0) {
        imp_one <- NA_real_
      }
      
      results[[idx]] <- tibble(
        feature_type = feature_type_name,
        domain = domain_name,
        feature_safe_id = feature_safe,
        env_var = ev,
        correlation = cor_one$r,
        p_value = cor_one$p,
        n = cor_one$n,
        importance_inc_mse = imp_one,
        explained_variation = explained_variation
      )
      
      idx <- idx + 1
    }
  }
  
  bind_rows(results)
}

analysis_results <- expand_grid(
  feature_type = unique(feature_long$feature_type),
  domain = domain_keep
) %>%
  mutate(result = map2(feature_type, domain, analyze_one_set)) %>%
  unnest(result)

if (nrow(analysis_results) == 0) {
  stop("No valid feature-environment model was generated.")
}

analysis_results <- analysis_results %>%
  group_by(feature_type, domain) %>%
  mutate(q_value = p.adjust(p_value, method = p_adjust_method)) %>%
  ungroup() %>%
  left_join(
    feature_info %>%
      select(
        feature_type,
        domain,
        feature_id,
        feature_safe_id,
        feature_label,
        feature_plot_label,
        total_abundance
      ),
    by = c("feature_type", "domain", "feature_safe_id")
  ) %>%
  mutate(
    star = p_to_star(p_value),
    tile_label = case_when(
      show_correlation_text & show_significance_stars ~ paste0(sprintf("%.2f", correlation), "\n", star),
      show_correlation_text ~ sprintf("%.2f", correlation),
      show_significance_stars ~ star,
      TRUE ~ ""
    ),
    importance_plot = ifelse(
      is.na(importance_inc_mse) | importance_inc_mse < 0,
      0,
      importance_inc_mse
    )
  )

write_tsv(
  analysis_results,
  file.path(output_dir, "env_feature_correlation_random_forest.tsv")
)

# ----------------------------
# 9. Plot one panel
# ----------------------------

plot_one_panel <- function(feature_type_name, domain_name) {
  df <- analysis_results %>%
    filter(feature_type == feature_type_name, domain == domain_name)
  
  if (nrow(df) == 0) {
    return(NULL)
  }
  
  feature_order <- df %>%
    group_by(feature_safe_id, feature_plot_label) %>%
    summarise(
      explained_variation = first(explained_variation),
      max_abs_r = max(abs(correlation), na.rm = TRUE),
      total_abundance = first(total_abundance),
      .groups = "drop"
    ) %>%
    arrange(desc(explained_variation), desc(max_abs_r), desc(total_abundance)) %>%
    mutate(feature_plot_label = factor(feature_plot_label, levels = feature_plot_label))
  
  df <- df %>%
    left_join(
      feature_order %>% select(feature_safe_id, feature_plot_label_order = feature_plot_label),
      by = "feature_safe_id"
    ) %>%
    mutate(
      feature_plot_label = factor(feature_plot_label, levels = levels(feature_order$feature_plot_label)),
      env_var = factor(env_var, levels = rev(env_vars))
    )
  
  bar_df <- feature_order %>%
    mutate(feature_plot_label = factor(feature_plot_label, levels = levels(feature_order$feature_plot_label)))
  
  y_max <- max(5, max(bar_df$explained_variation, na.rm = TRUE) * 1.15)
  
  p_bar <- ggplot(
    bar_df,
    aes(x = feature_plot_label, y = explained_variation)
  ) +
    geom_col(
      width = 0.65,
      fill = "#4DBBD5",
      color = "black",
      linewidth = 0.25
    ) +
    scale_y_continuous(
      limits = c(0, y_max),
      expand = expansion(mult = c(0, 0.03))
    ) +
    labs(
      y = "Percent explained\nvariation",
      x = NULL,
      title = domain_name
    ) +
    theme_nature(base_size = 7) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  p_heat <- ggplot(
    df,
    aes(x = feature_plot_label, y = env_var, fill = correlation)
  ) +
    geom_tile(color = "grey88", linewidth = 0.2) +
    geom_point(
      aes(size = importance_plot),
      shape = 21,
      stroke = 0.45,
      color = "black",
      fill = NA
    ) +
    geom_text(
      aes(label = tile_label),
      size = 2.0,
      lineheight = 0.75
    ) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-1, 1),
      name = "Correlation"
    ) +
    scale_size_continuous(
      range = c(0.5, 4.5),
      name = "Importance\n(%IncMSE)"
    ) +
    labs(
      x = NULL,
      y = NULL
    ) +
    theme_nature(base_size = 7) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
      legend.position = "right",
      panel.border = element_rect(color = "grey60", fill = NA, linewidth = 0.25)
    )
  
  p_bar / p_heat + plot_layout(heights = c(0.45, 1.6))
}

# ----------------------------
# 10. Build combined figure
# ----------------------------

panel_list <- list()

feature_types <- c("Taxa", "Function")
feature_types <- feature_types[feature_types %in% unique(feature_long$feature_type)]

for (ft in feature_types) {
  for (dm in domain_keep) {
    p <- plot_one_panel(ft, dm)
    if (!is.null(p)) {
      panel_list[[paste(ft, dm, sep = "_")]] <- p
    }
  }
}

if (length(panel_list) == 0) {
  stop("No panel was generated.")
}

if (length(feature_types) == 2) {
  taxa_panels <- panel_list[stringr::str_starts(names(panel_list), "Taxa_")]
  func_panels <- panel_list[stringr::str_starts(names(panel_list), "Function_")]
  
  fig_taxa <- wrap_plots(taxa_panels, nrow = 1, guides = "collect") +
    plot_annotation(title = "Taxonomic features")
  
  fig_func <- wrap_plots(func_panels, nrow = 1, guides = "collect") +
    plot_annotation(title = "Functional features")
  
  final_fig <- fig_taxa / fig_func +
    plot_layout(heights = c(1, 1), guides = "collect") +
    plot_annotation(
      tag_levels = "a",
      theme = theme(
        plot.tag = element_text(face = "bold", size = 12),
        legend.position = "right"
      )
    )
  
  fig_height_mm <- figure_height_taxa_function_mm
} else {
  final_fig <- wrap_plots(panel_list, nrow = 1, guides = "collect") +
    plot_annotation(
      tag_levels = "a",
      theme = theme(
        plot.tag = element_text(face = "bold", size = 12),
        legend.position = "right"
      )
    )
  
  fig_height_mm <- figure_height_taxa_only_mm
}

final_fig

# ----------------------------
# 11. Save figure
# ----------------------------

ggsave(
  filename = file.path(output_dir, "env_explained_taxa_functions.pdf"),
  plot = final_fig,
  width = figure_width_mm,
  height = fig_height_mm,
  units = "mm",
  device = "pdf"
)

ggsave(
  filename = file.path(output_dir, "env_explained_taxa_functions.png"),
  plot = final_fig,
  width = figure_width_mm,
  height = fig_height_mm,
  units = "mm",
  dpi = 600
)

message("[DONE] Results written to: ", output_dir)
