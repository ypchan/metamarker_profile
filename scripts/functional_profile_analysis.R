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
  `table-format` = "feature_rows",
  `feature-col` = "function_id",
  group = "",
  covariates = "",
  outdir = "analysis_out/functions"
)

opt <- mc_parse_args(list(
  table = "", metadata = "", `table-format` = "feature_rows", `sample-col` = "sample_id",
  `feature-col` = "function_id", `metadata-cols` = "", group = "", covariates = "",
  outdir = "analysis_out/functions", `min-prevalence` = "0.10", `min-total` = "0",
  `top-n` = "30", permutations = "9999", pseudocount = "1e-06", `fig-format` = "both"
))

mc_need_packages(c("vegan", "ggplot2"))
dirs <- mc_make_outdirs(opt$outdir)
log_file <- file.path(dirs$logs, "functional_profile_analysis.log")

dat <- mc_read_feature_data(opt$table, opt$metadata, opt$`table-format`, opt$`sample-col`,
                            opt$`feature-col`, opt$`metadata-cols`)
mat <- mc_filter_features(dat$abund, as.numeric(opt$`min-prevalence`), as.numeric(opt$`min-total`))
if (ncol(mat) < 2) stop("Fewer than two functions remain after filtering.")
meta <- dat$metadata
if (opt$group != "" && !opt$group %in% names(meta)) stop("Group column not found: ", opt$group)

rel <- mc_normalize(mat, "relative")
top_features <- mc_rank_columns(rel, as.integer(opt$`top-n`))
means <- data.frame(function_id = colnames(rel), mean_relative_abundance = colMeans(rel),
                    prevalence = colMeans(rel > 0), stringsAsFactors = FALSE)
means <- means[order(-means$mean_relative_abundance), ]
mc_write_tsv(means, file.path(dirs$tables, "function_summary.tsv"))

top <- means[means$function_id %in% top_features, ]
top$function_id <- factor(top$function_id, levels = rev(top$function_id))
p_bar <- ggplot2::ggplot(top, ggplot2::aes(x = function_id, y = mean_relative_abundance)) +
  ggplot2::geom_col(fill = "#3B6FB6", width = 0.75) +
  ggplot2::coord_flip() +
  ggplot2::labs(x = NULL, y = "Mean relative abundance", title = "Top predicted functions") +
  mc_theme()
mc_save_plot(p_bar, file.path(dirs$figures, "top_functions"), 8, 6, opt$`fig-format`)

dist_obj <- vegan::vegdist(rel, method = "bray")
pcoa <- stats::cmdscale(dist_obj, k = 2, eig = TRUE)
eig <- pcoa$eig[pcoa$eig > 0]
var_exp <- if (length(eig) >= 2) round(100 * eig[1:2] / sum(eig), 1) else c(NA, NA)
ord <- data.frame(sample_id = rownames(pcoa$points), Axis1 = pcoa$points[, 1], Axis2 = pcoa$points[, 2],
                  stringsAsFactors = FALSE)
names(ord)[1] <- dat$sample_col
ord <- merge(meta, ord, by = dat$sample_col, sort = FALSE)
ord <- ord[match(rownames(rel), ord[[dat$sample_col]]), , drop = FALSE]
mc_write_tsv(ord, file.path(dirs$tables, "function_pcoa_bray.tsv"))

if (opt$group != "") {
  p_ord <- ggplot2::ggplot(ord, ggplot2::aes(x = Axis1, y = Axis2, color = .data[[opt$group]])) +
    ggplot2::geom_point(size = 3, alpha = 0.9) +
    ggplot2::stat_ellipse(linewidth = 0.5, alpha = 0.8, show.legend = FALSE) +
    ggplot2::scale_color_manual(values = mc_palette(length(unique(ord[[opt$group]])))) +
    ggplot2::labs(x = paste0("PCoA1 (", var_exp[1], "%)"), y = paste0("PCoA2 (", var_exp[2], "%)"),
                  color = opt$group, title = "Predicted function PCoA")
} else {
  p_ord <- ggplot2::ggplot(ord, ggplot2::aes(x = Axis1, y = Axis2)) +
    ggplot2::geom_point(size = 3, alpha = 0.9, color = "#3B6FB6") +
    ggplot2::labs(x = paste0("PCoA1 (", var_exp[1], "%)"), y = paste0("PCoA2 (", var_exp[2], "%)"),
                  title = "Predicted function PCoA")
}
mc_save_plot(p_ord + mc_theme(), file.path(dirs$figures, "function_pcoa_bray"), 6.4, 5.2, opt$`fig-format`)

if (opt$group != "") {
  covars <- mc_split(opt$covariates)
  keep <- mc_complete_cases(meta, c(opt$group, covars))
  dist_sub <- vegan::vegdist(rel[keep, , drop = FALSE], method = "bray")
  meta_sub <- meta[keep, , drop = FALSE]
  form <- stats::as.formula(paste("dist_sub ~", paste(c(opt$group, covars), collapse = " + ")))
  ad <- vegan::adonis2(form, data = meta_sub, permutations = as.integer(opt$permutations), by = "margin")
  ad_df <- data.frame(term = rownames(as.data.frame(ad)), as.data.frame(ad), check.names = FALSE)
  mc_write_tsv(ad_df, file.path(dirs$tables, "function_permanova_bray.tsv"))

  clr <- mc_normalize(mat[keep, , drop = FALSE], "clr", as.numeric(opt$pseudocount))
  meta_lm <- meta[keep, , drop = FALSE]
  meta_lm[[opt$group]] <- factor(meta_lm[[opt$group]])
  rows <- list()
  for (j in seq_len(ncol(clr))) {
    df <- data.frame(y = clr[, j], meta_lm, check.names = FALSE)
    fit <- stats::lm(stats::as.formula(paste("y ~", paste(c(opt$group, covars), collapse = " + "))), data = df)
    cf <- as.data.frame(summary(fit)$coefficients)
    cf$term <- rownames(cf)
    cf <- cf[grepl(paste0("^", opt$group), cf$term), , drop = FALSE]
    if (nrow(cf) > 0) {
      rows[[length(rows) + 1]] <- data.frame(function_id = colnames(clr)[j],
                                             contrast = sub(paste0("^", opt$group), "", cf$term),
                                             estimate_clr = cf$Estimate,
                                             p_value = cf$`Pr(>|t|)`,
                                             stringsAsFactors = FALSE)
    }
  }
  if (length(rows) > 0) {
    diff <- do.call(rbind, rows)
    diff$q_value <- p.adjust(diff$p_value, method = "BH")
    diff <- diff[order(diff$q_value, diff$p_value), ]
    mc_write_tsv(diff, file.path(dirs$tables, "function_differential_clr_lm.tsv"))
  }
}

mc_log("Done. Output: ", opt$outdir, log_file = log_file)
