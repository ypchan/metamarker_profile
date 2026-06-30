
# Mantel test + Spearman correlation heatmap ----

# Author: Yanpeng Chen

# 0. Packages
library(tidyverse)
library(vegan)
library(linkET)
library(ggsci)
library(RColorBrewer)
library(grid)

devtools::install_github("Hy4m/linkET", force = TRUE)

# 1. User parameters


env_file  <- "env.tsv"
comm_file <- "community_long.tsv"

out_dir <- "mantel_correlation_heatmap"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 群落矩阵距离方法
# 常用：
# "bray"       适合群落丰度
# "jaccard"    适合 presence/absence
comm_dist_method <- "bray"

# 环境因子距离方法
# 对单个环境变量通常用 euclidean
env_dist_method <- "euclidean"

# Mantel permutation
n_perm <- 999

# 是否对群落丰度做 Hellinger transformation
# 推荐 TRUE，尤其是生态群落丰度数据
use_hellinger <- TRUE

# 是否对环境因子做 scale
# 推荐 TRUE，特别是 pH、TOC、Fe、Mn 量纲不同
scale_env <- TRUE

# 显示哪些 domain
target_domains <- c("Bacteria", "Archaea", "Fungi")


# 2. Read input data


env_raw <- readr::read_tsv(env_file, show_col_types = FALSE)

comm_raw <- readr::read_tsv(comm_file, show_col_types = FALSE)


# 3. Check input format


required_env_cols <- c("sample_id")
required_comm_cols <- c("sample_id", "domain", "taxon", "abundance")

if (!all(required_env_cols %in% colnames(env_raw))) {
  stop("env.tsv must contain column: sample_id")
}

if (!all(required_comm_cols %in% colnames(comm_raw))) {
  stop("community_long.tsv must contain columns: sample_id, domain, taxon, abundance")
}

if (anyDuplicated(env_raw$sample_id) > 0) {
  stop("env.tsv contains duplicated sample_id")
}

# 检查环境因子是否都是 numeric
env_factor_cols <- setdiff(colnames(env_raw), "sample_id")

non_numeric_env <- env_raw %>%
  select(all_of(env_factor_cols)) %>%
  summarise(across(everything(), ~ !is.numeric(.x))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "non_numeric") %>%
  filter(non_numeric) %>%
  pull(variable)

if (length(non_numeric_env) > 0) {
  stop(
    "These environmental columns are not numeric: ",
    paste(non_numeric_env, collapse = ", ")
  )
}


# 4. Prepare shared samples


shared_samples <- intersect(env_raw$sample_id, comm_raw$sample_id)

if (length(shared_samples) < 5) {
  stop("Too few shared samples between env.tsv and community_long.tsv")
}

env_tbl <- env_raw %>%
  filter(sample_id %in% shared_samples) %>%
  arrange(sample_id)

comm_tbl <- comm_raw %>%
  filter(sample_id %in% shared_samples) %>%
  mutate(
    domain = as.character(domain),
    taxon = as.character(taxon),
    abundance = as.numeric(abundance)
  ) %>%
  filter(!is.na(abundance))

# 只保留目标 domain
comm_tbl <- comm_tbl %>%
  filter(domain %in% target_domains)


# 5. Prepare environmental matrix


env_mat <- env_tbl %>%
  column_to_rownames("sample_id") %>%
  as.data.frame()

# 删除全 NA 或零方差环境因子
env_mat <- env_mat[, colSums(is.na(env_mat)) == 0, drop = FALSE]

env_sd <- apply(env_mat, 2, sd, na.rm = TRUE)
env_mat <- env_mat[, env_sd > 0, drop = FALSE]

if (scale_env) {
  env_mat <- scale(env_mat) %>%
    as.data.frame()
}


# 6. Prepare community wide matrix


comm_wide <- comm_tbl %>%
  mutate(
    feature_id = paste(domain, taxon, sep = "__")
  ) %>%
  group_by(sample_id, feature_id) %>%
  summarise(
    abundance = sum(abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = feature_id,
    values_from = abundance,
    values_fill = 0
  ) %>%
  arrange(sample_id)

# 与 env_mat 再次对齐
common_samples2 <- intersect(rownames(env_mat), comm_wide$sample_id)

env_mat <- env_mat[common_samples2, , drop = FALSE]

comm_mat <- comm_wide %>%
  filter(sample_id %in% common_samples2) %>%
  arrange(match(sample_id, common_samples2)) %>%
  column_to_rownames("sample_id") %>%
  as.data.frame()

# 删除全 0 taxa
comm_mat <- comm_mat[, colSums(comm_mat, na.rm = TRUE) > 0, drop = FALSE]

# Hellinger transformation
if (use_hellinger) {
  comm_mat <- vegan::decostand(comm_mat, method = "hellinger") %>%
    as.data.frame()
}


# 7. Build domain group index for Mantel test


spec_select <- list()

for (dm in target_domains) {
  idx <- which(str_starts(colnames(comm_mat), paste0(dm, "__")))
  if (length(idx) > 0) {
    spec_select[[dm]] <- idx
  }
}

if (length(spec_select) == 0) {
  stop("No domain-specific features found in community matrix.")
}

message("[INFO] Samples used: ", nrow(comm_mat))
message("[INFO] Environmental variables: ", ncol(env_mat))
message("[INFO] Community features: ", ncol(comm_mat))
message("[INFO] Domain groups: ", paste(names(spec_select), collapse = ", "))


# 8. Spearman correlation among environmental factors


env_cor <- linkET::correlate(
  env_mat,
  method = "spearman"
)

# 输出环境因子 Spearman 相关矩阵
env_cor_tbl <- cor(
  env_mat,
  method = "spearman",
  use = "pairwise.complete.obs"
) %>%
  as.data.frame() %>%
  rownames_to_column("env_var_1") %>%
  pivot_longer(
    cols = -env_var_1,
    names_to = "env_var_2",
    values_to = "spearman_r"
  )

readr::write_tsv(
  env_cor_tbl,
  file.path(out_dir, "environment_spearman_correlation.tsv")
)


# 9. Mantel test


mantel_res <- linkET::mantel_test(
  spec = comm_mat,
  env = env_mat,
  spec_select = spec_select,
  method = "spearman",
  permutations = n_perm,
  dist = comm_dist_method,
  env_dist = env_dist_method
)

# 整理 Mantel 结果
mantel_res2 <- mantel_res %>%
  mutate(
    mantel_r_abs = abs(r),
    mantel_sign = if_else(r >= 0, "Positive", "Negative"),
    mantel_r_class = cut(
      mantel_r_abs,
      breaks = c(-Inf, 0.2, 0.4, Inf),
      labels = c("< 0.2", "0.2–0.4", "≥ 0.4")
    ),
    mantel_p_class = cut(
      p,
      breaks = c(-Inf, 0.01, 0.05, Inf),
      labels = c("< 0.01", "0.01–0.05", "≥ 0.05")
    )
  )

readr::write_tsv(
  mantel_res2,
  file.path(out_dir, "mantel_test_results.tsv")
)


# 10. Nature-style color theme


# Nature-like diverging palette for Spearman correlation
cor_cols <- c(
  "#3B4CC0",  # blue
  "#6F91D2",
  "#BFD3E6",
  "#F7F7F7",
  "#F4A582",
  "#D6604D",
  "#8B1A1A"   # dark red
)

# Mantel p-value colors
mantel_p_cols <- c(
  "< 0.01"    = "#C56E33",  # orange-brown
  "0.01–0.05" = "#4C9A7A",  # green
  "≥ 0.05"    = "#BDBDBD"   # grey
)

# Mantel line width
mantel_size_vals <- c(
  "< 0.2"   = 0.35,
  "0.2–0.4" = 0.75,
  "≥ 0.4"   = 1.35
)

# Mantel line type
mantel_lty_vals <- c(
  "Positive" = "solid",
  "Negative" = "dashed"
)

theme_nature_mantel <- function(base_size = 11, base_family = "Arial") {
  theme(
    text = element_text(family = base_family, color = "grey15"),
    plot.title = element_text(
      size = base_size + 2,
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      size = base_size,
      hjust = 0.5,
      color = "grey35"
    ),
    axis.text = element_text(size = base_size - 1, color = "grey20"),
    axis.title = element_blank(),
    legend.title = element_text(size = base_size, face = "bold"),
    legend.text = element_text(size = base_size - 1),
    legend.key = element_blank(),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(8, 8, 8, 8)
  )
}


# 11. Plot


p_mantel <- qcorrplot(
  env_cor,
  type = "lower",
  diag = FALSE
) +
  geom_square() +
  geom_mark(
    size = 3,
    only_mark = TRUE,
    sig_level = c(0.05, 0.01, 0.001),
    sig_thres = 0.05
  ) +
  geom_couple(
    data = mantel_res2,
    aes(
      colour = mantel_p_class,
      size = mantel_r_class,
      linetype = mantel_sign
    ),
    curvature = nice_curvature()
  ) +
  scale_fill_gradientn(
    colours = cor_cols,
    limits = c(-1, 1),
    name = "Spearman's r"
  ) +
  scale_colour_manual(
    values = mantel_p_cols,
    name = "Mantel's p"
  ) +
  scale_size_manual(
    values = mantel_size_vals,
    name = "Mantel's |r|"
  ) +
  scale_linetype_manual(
    values = mantel_lty_vals,
    name = "Mantel's sign"
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(3.0, "cm"),
      barheight = unit(0.35, "cm")
    ),
    colour = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      override.aes = list(size = 1.2)
    ),
    size = guide_legend(
      title.position = "top",
      title.hjust = 0.5
    ),
    linetype = guide_legend(
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  labs(
    title = "Mantel test and environmental correlation",
    subtitle = paste0(
      "Community distance: ", comm_dist_method,
      "; environmental distance: ", env_dist_method,
      "; permutations: ", n_perm
    )
  ) +
  theme_nature_mantel(base_size = 11)

print(p_mantel)


# 12. Save figure


ggsave(
  filename = file.path(out_dir, "mantel_correlation_heatmap.pdf"),
  plot = p_mantel,
  width = 9.5,
  height = 6.8,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = file.path(out_dir, "mantel_correlation_heatmap.png"),
  plot = p_mantel,
  width = 9.5,
  height = 6.8,
  units = "in",
  dpi = 600
)

message("[DONE] ", file.path(out_dir, "mantel_correlation_heatmap.pdf"))
message("[DONE] ", file.path(out_dir, "mantel_correlation_heatmap.png"))
message("[DONE] ", file.path(out_dir, "mantel_test_results.tsv"))
message("[DONE] ", file.path(out_dir, "environment_spearman_correlation.tsv"))