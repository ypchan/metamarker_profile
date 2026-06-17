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
  reference = "",
  covariates = "",
  outdir = "analysis_out/differential_abundance"
)

opt <- mc_parse_args(list(
  table = "", metadata = "", `table-format` = "sample_rows", `sample-col` = "sample_id",
  `feature-col` = "feature_id", `metadata-cols` = "", group = "", reference = "",
  covariates = "", outdir = "analysis_out/differential_abundance",
  `min-prevalence` = "0.10", `min-total` = "0", pseudocount = "1e-06",
  alpha = "0.05", `top-n` = "30", `fig-format` = "both"
))

mc_need_packages(c("ggplot2"))
if (opt$group == "") stop("Set MC_CONFIG$group before sourcing this template.")
dirs <- mc_make_outdirs(opt$outdir)
log_file <- file.path(dirs$logs, "differential_abundance_compositional.log")

dat <- mc_read_feature_data(opt$table, opt$metadata, opt$`table-format`, opt$`sample-col`,
                            opt$`feature-col`, opt$`metadata-cols`)
mat <- mc_filter_features(dat$abund, as.numeric(opt$`min-prevalence`), as.numeric(opt$`min-total`))
if (ncol(mat) < 2) stop("Fewer than two features remain after filtering.")
meta <- dat$metadata
sample_col <- dat$sample_col
if (!opt$group %in% names(meta)) stop("Group column not found: ", opt$group)
covars <- mc_split(opt$covariates)
missing_cov <- setdiff(covars, names(meta))
if (length(missing_cov) > 0) stop("Covariate column(s) not found: ", paste(missing_cov, collapse = ", "))

keep <- mc_complete_cases(meta, c(opt$group, covars))
meta <- meta[keep, , drop = FALSE]
mat <- mat[keep, , drop = FALSE]
group <- factor(meta[[opt$group]])
if (opt$reference != "") group <- stats::relevel(group, ref = opt$reference)
if (nlevels(group) < 2) stop("The group column must contain at least two levels after filtering.")
meta[[opt$group]] <- group

rel <- mc_normalize(mat, "relative")
clr <- mc_normalize(mat, "clr", as.numeric(opt$pseudocount))
feature_anno <- mc_feature_annotation(colnames(mat))

overall <- vector("list", ncol(clr))
coef_rows <- list()
for (j in seq_len(ncol(clr))) {
  df <- data.frame(y = clr[, j], meta, check.names = FALSE)
  form <- stats::as.formula(paste("y ~", paste(c(opt$group, covars), collapse = " + ")))
  fit <- stats::lm(form, data = df)
  a <- stats::anova(fit)
  p_group <- if (opt$group %in% rownames(a)) a[opt$group, "Pr(>F)"] else NA_real_
  overall[[j]] <- data.frame(feature_id = colnames(clr)[j], p_value = p_group, stringsAsFactors = FALSE)
  cf <- as.data.frame(summary(fit)$coefficients)
  cf$term <- rownames(cf)
  cf <- cf[grepl(paste0("^", opt$group), cf$term), , drop = FALSE]
  if (nrow(cf) > 0) {
    coef_rows[[length(coef_rows) + 1]] <- data.frame(
      feature_id = colnames(clr)[j],
      contrast = sub(paste0("^", opt$group), "", cf$term),
      estimate_clr = cf$Estimate,
      std_error = cf$`Std. Error`,
      statistic = cf$`t value`,
      p_value = cf$`Pr(>|t|)`,
      stringsAsFactors = FALSE
    )
  }
}

overall <- do.call(rbind, overall)
overall$q_value <- p.adjust(overall$p_value, method = "BH")
overall <- merge(overall, feature_anno, by = "feature_id", all.x = TRUE, sort = FALSE)
overall <- overall[order(overall$q_value, overall$p_value), ]
mc_write_tsv(overall, file.path(dirs$tables, "overall_group_tests.tsv"))

if (length(coef_rows) > 0) {
  coef_tab <- do.call(rbind, coef_rows)
  coef_tab$q_value <- p.adjust(coef_tab$p_value, method = "BH")
  means <- aggregate(rel, by = list(group = group), FUN = mean)
  mean_long <- reshape(means, varying = colnames(mat), v.names = "mean_relative_abundance",
                       timevar = "feature_id", times = colnames(mat), direction = "long")
  mean_long$group <- as.character(mean_long$group)
  mc_write_tsv(mean_long[, c("group", "feature_id", "mean_relative_abundance")],
               file.path(dirs$tables, "group_mean_relative_abundance.tsv"))
  coef_tab <- merge(coef_tab, feature_anno, by = "feature_id", all.x = TRUE, sort = FALSE)
  coef_tab <- coef_tab[order(coef_tab$q_value, coef_tab$p_value), ]
  mc_write_tsv(coef_tab, file.path(dirs$tables, "group_contrasts.tsv"))

  alpha <- as.numeric(opt$alpha)
  coef_tab$significant <- coef_tab$q_value <= alpha
  p_vol <- ggplot2::ggplot(coef_tab, ggplot2::aes(x = estimate_clr, y = -log10(q_value), color = significant)) +
    ggplot2::geom_point(size = 2, alpha = 0.8) +
    ggplot2::facet_wrap(~ contrast, scales = "free_x") +
    ggplot2::scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "#D55E00")) +
    ggplot2::labs(x = "CLR effect estimate", y = "-log10(FDR)", color = paste0("FDR <= ", alpha),
                  title = "Differential abundance contrasts") +
    mc_theme()
  mc_save_plot(p_vol, file.path(dirs$figures, "contrast_volcano"), 8, 5, opt$`fig-format`)

  top <- head(coef_tab[order(coef_tab$q_value, -abs(coef_tab$estimate_clr)), ], as.integer(opt$`top-n`))
  if (nrow(top) > 0) {
    top$label <- ifelse(is.na(top$taxon), top$feature_id, top$taxon)
    top$label <- factor(top$label, levels = rev(unique(top$label)))
    p_bar <- ggplot2::ggplot(top, ggplot2::aes(x = label, y = estimate_clr, fill = contrast)) +
      ggplot2::geom_col(width = 0.75) +
      ggplot2::coord_flip() +
      ggplot2::labs(x = NULL, y = "CLR effect estimate", fill = "Contrast",
                    title = "Top differential features") +
      ggplot2::scale_fill_manual(values = mc_palette(length(unique(top$contrast)))) +
      mc_theme()
    mc_save_plot(p_bar, file.path(dirs$figures, "top_differential_features"), 8, 6, opt$`fig-format`)
  }
}

mc_log("Done. Output: ", opt$outdir, log_file = log_file)
