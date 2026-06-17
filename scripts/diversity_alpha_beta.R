#!/usr/bin/env Rscript

.script_path <- tryCatch(
  if (requireNamespace("rstudioapi", quietly = TRUE)) rstudioapi::getActiveDocumentContext()$path else "",
  error = function(e) ""
)
if (!nzchar(.script_path)) .script_path <- "scripts"
.script_dir <- if (dir.exists(.script_path)) .script_path else dirname(normalizePath(.script_path, mustWork = FALSE))
source(file.path(.script_dir, "analysis_common.R"))

# RStudio mode: edit this block, then click Source/Run.
MC_CONFIG <- list(
  table = "",
  metadata = "",
  `table-format` = "sample_rows",
  `metadata-cols` = "",
  group = "",
  covariates = "",
  strata = "",
  outdir = "analysis_out/diversity_alpha_beta"
)

opt <- mc_parse_args(list(
  table = "", metadata = "", `table-format` = "sample_rows", `sample-col` = "sample_id",
  `feature-col` = "feature_id", `metadata-cols` = "", group = "", covariates = "",
  strata = "", outdir = "analysis_out/diversity_alpha_beta", normalize = "relative",
  distance = "bray", `min-prevalence` = "0.10", `min-total` = "0",
  permutations = "9999", seed = "1", `fig-format` = "both"
))

mc_need_packages(c("vegan", "ggplot2"))
set.seed(as.integer(opt$seed))
dirs <- mc_make_outdirs(opt$outdir)
log_file <- file.path(dirs$logs, "diversity_alpha_beta.log")

dat <- mc_read_feature_data(opt$table, opt$metadata, opt$`table-format`, opt$`sample-col`,
                            opt$`feature-col`, opt$`metadata-cols`)
mat0 <- mc_filter_features(dat$abund, as.numeric(opt$`min-prevalence`), as.numeric(opt$`min-total`))
if (ncol(mat0) < 2) stop("Fewer than two features remain after filtering.")
meta <- dat$metadata
sample_col <- dat$sample_col
if (opt$group != "" && !opt$group %in% names(meta)) stop("Group column not found: ", opt$group)

mat_alpha <- mc_normalize(mat0, opt$normalize)
alpha <- data.frame(
  sample_id = rownames(mat_alpha),
  observed = rowSums(mat_alpha > 0),
  shannon = vegan::diversity(mat_alpha, index = "shannon"),
  simpson = vegan::diversity(mat_alpha, index = "simpson"),
  invsimpson = vegan::diversity(mat_alpha, index = "invsimpson"),
  stringsAsFactors = FALSE
)
alpha$pielou <- ifelse(alpha$observed > 1, alpha$shannon / log(alpha$observed), NA_real_)
names(alpha)[1] <- sample_col
alpha <- merge(meta, alpha, by = sample_col, sort = FALSE)
alpha <- alpha[match(rownames(mat_alpha), alpha[[sample_col]]), , drop = FALSE]
mc_write_tsv(alpha, file.path(dirs$tables, "alpha_diversity.tsv"))

metrics <- c("observed", "shannon", "simpson", "invsimpson", "pielou")
if (opt$group != "") {
  stats <- do.call(rbind, lapply(metrics, function(m) {
    keep <- complete.cases(alpha[, c(m, opt$group)])
    g <- factor(alpha[[opt$group]][keep])
    p <- if (nlevels(g) > 1) stats::kruskal.test(alpha[[m]][keep] ~ g)$p.value else NA_real_
    data.frame(metric = m, test = "Kruskal-Wallis", p_value = p, stringsAsFactors = FALSE)
  }))
  stats$p_adjust <- p.adjust(stats$p_value, method = "BH")
  mc_write_tsv(stats, file.path(dirs$tables, "alpha_group_tests.tsv"))
}

alpha_cols <- c(sample_col, metrics[metrics %in% names(alpha)])
if (opt$group != "") alpha_cols <- c(sample_col, opt$group, metrics[metrics %in% names(alpha)])
alpha_long <- reshape(alpha[, alpha_cols, drop = FALSE],
                      varying = metrics[metrics %in% names(alpha)], v.names = "value", timevar = "metric",
                      times = metrics, direction = "long")
alpha_long$metric <- factor(alpha_long$metric, levels = metrics)
if (opt$group != "") {
  p_alpha <- ggplot2::ggplot(alpha_long, ggplot2::aes(x = .data[[opt$group]], y = .data$value, fill = .data[[opt$group]])) +
    ggplot2::geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.75) +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.12, height = 0), size = 1.8, alpha = 0.8) +
    ggplot2::facet_wrap(~ metric, scales = "free_y", nrow = 1) +
    ggplot2::scale_fill_manual(values = mc_palette(length(unique(alpha_long[[opt$group]])))) +
    ggplot2::labs(x = opt$group, y = "Alpha diversity", title = "Alpha diversity") +
    mc_theme() + ggplot2::theme(legend.position = "none")
} else {
  p_alpha <- ggplot2::ggplot(alpha_long, ggplot2::aes(x = metric, y = value)) +
    ggplot2::geom_boxplot(width = 0.65, outlier.shape = NA, fill = "#3B6FB6", alpha = 0.75) +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.12, height = 0), size = 1.8, alpha = 0.8) +
    ggplot2::labs(x = NULL, y = "Alpha diversity", title = "Alpha diversity") + mc_theme()
}
mc_save_plot(p_alpha, file.path(dirs$figures, "alpha_diversity"), 9, 4.8, opt$`fig-format`)

mat_beta <- mc_normalize(mat0, if (opt$distance == "bray") "relative" else opt$normalize)
dist_obj <- vegan::vegdist(mat_beta, method = opt$distance)
pcoa <- stats::cmdscale(dist_obj, k = 2, eig = TRUE)
eig <- pcoa$eig[pcoa$eig > 0]
var_exp <- if (length(eig) >= 2) round(100 * eig[1:2] / sum(eig), 1) else c(NA, NA)
ord <- data.frame(sample_id = rownames(pcoa$points), Axis1 = pcoa$points[, 1], Axis2 = pcoa$points[, 2],
                  stringsAsFactors = FALSE)
names(ord)[1] <- sample_col
ord <- merge(meta, ord, by = sample_col, sort = FALSE)
ord <- ord[match(rownames(mat_beta), ord[[sample_col]]), , drop = FALSE]
mc_write_tsv(ord, file.path(dirs$tables, paste0("pcoa_", opt$distance, ".tsv")))

if (opt$group != "") {
  p_ord <- ggplot2::ggplot(ord, ggplot2::aes(x = Axis1, y = Axis2, color = .data[[opt$group]])) +
    ggplot2::geom_point(size = 3, alpha = 0.9) +
    ggplot2::stat_ellipse(linewidth = 0.5, alpha = 0.8, show.legend = FALSE) +
    ggplot2::scale_color_manual(values = mc_palette(length(unique(ord[[opt$group]])))) +
    ggplot2::labs(x = paste0("PCoA1 (", var_exp[1], "%)"), y = paste0("PCoA2 (", var_exp[2], "%)"),
                  color = opt$group, title = paste("PCoA -", opt$distance))
} else {
  p_ord <- ggplot2::ggplot(ord, ggplot2::aes(x = Axis1, y = Axis2)) +
    ggplot2::geom_point(size = 3, alpha = 0.9, color = "#3B6FB6") +
    ggplot2::labs(x = paste0("PCoA1 (", var_exp[1], "%)"), y = paste0("PCoA2 (", var_exp[2], "%)"),
                  title = paste("PCoA -", opt$distance))
}
p_ord <- p_ord + mc_theme()
mc_save_plot(p_ord, file.path(dirs$figures, paste0("pcoa_", opt$distance)), 6.4, 5.2, opt$`fig-format`)

model_vars <- c(opt$group, mc_split(opt$covariates))
model_vars <- model_vars[model_vars != ""]
if (length(model_vars) > 0) {
  keep <- mc_complete_cases(meta, c(model_vars, opt$strata))
  dist_sub <- vegan::vegdist(mat_beta[keep, , drop = FALSE], method = opt$distance)
  meta_sub <- meta[keep, , drop = FALSE]
  form <- stats::as.formula(paste("dist_sub ~", paste(model_vars, collapse = " + ")))
  strata_vec <- if (opt$strata != "") meta_sub[[opt$strata]] else NULL
  perm <- as.integer(opt$permutations)
  ad <- vegan::adonis2(form, data = meta_sub, permutations = perm, by = "margin", strata = strata_vec)
  ad_df <- data.frame(term = rownames(as.data.frame(ad)), as.data.frame(ad), check.names = FALSE)
  mc_write_tsv(ad_df, file.path(dirs$tables, paste0("permanova_", opt$distance, ".tsv")))
}

if (opt$group != "") {
  g <- factor(meta[[opt$group]])
  if (nlevels(g) > 1) {
    bd <- vegan::betadisper(dist_obj, group = g)
    bd_perm <- vegan::permutest(bd, permutations = as.integer(opt$permutations))
    bd_df <- data.frame(term = rownames(as.data.frame(bd_perm$tab)), as.data.frame(bd_perm$tab), check.names = FALSE)
    mc_write_tsv(bd_df, file.path(dirs$tables, paste0("beta_dispersion_", opt$distance, ".tsv")))
  }
}

mc_log("Done. Output: ", opt$outdir, log_file = log_file)
