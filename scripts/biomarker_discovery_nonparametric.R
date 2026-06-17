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
  outdir = "analysis_out/biomarkers"
)

opt <- mc_parse_args(list(
  table = "", metadata = "", `table-format` = "sample_rows", `sample-col` = "sample_id",
  `feature-col` = "feature_id", `metadata-cols` = "", group = "",
  outdir = "analysis_out/biomarkers", `min-prevalence` = "0.10", `min-total` = "0",
  pseudocount = "1e-06", alpha = "0.05", `min-abs-log2fc` = "1",
  `top-n` = "30", `fig-format` = "both"
))

mc_need_packages(c("ggplot2"))
if (opt$group == "") stop("Set MC_CONFIG$group before sourcing this template.")
dirs <- mc_make_outdirs(opt$outdir)
log_file <- file.path(dirs$logs, "biomarker_discovery_nonparametric.log")

dat <- mc_read_feature_data(opt$table, opt$metadata, opt$`table-format`, opt$`sample-col`,
                            opt$`feature-col`, opt$`metadata-cols`)
mat <- mc_filter_features(dat$abund, as.numeric(opt$`min-prevalence`), as.numeric(opt$`min-total`))
if (ncol(mat) < 2) stop("Fewer than two features remain after filtering.")
meta <- dat$metadata
if (!opt$group %in% names(meta)) stop("Group column not found: ", opt$group)
keep <- mc_complete_cases(meta, opt$group)
meta <- meta[keep, , drop = FALSE]
mat <- mat[keep, , drop = FALSE]
group <- factor(meta[[opt$group]])
if (nlevels(group) < 2) stop("The group column must contain at least two levels.")

rel <- mc_normalize(mat, "relative")
clr <- mc_normalize(mat, "clr", as.numeric(opt$pseudocount))
anno <- mc_feature_annotation(colnames(mat))
levels_g <- levels(group)

rows <- vector("list", ncol(rel))
for (j in seq_len(ncol(rel))) {
  x <- rel[, j]
  p <- stats::kruskal.test(x ~ group)$p.value
  means <- tapply(x, group, mean)
  top_group <- names(which.max(means))
  other_mean <- mean(x[group != top_group])
  log2fc <- log2((means[[top_group]] + as.numeric(opt$pseudocount)) /
                   (other_mean + as.numeric(opt$pseudocount)))
  rows[[j]] <- data.frame(
    feature_id = colnames(rel)[j],
    enriched_group = top_group,
    mean_enriched = means[[top_group]],
    mean_other = other_mean,
    log2fc_enriched_vs_other = log2fc,
    p_value = p,
    stringsAsFactors = FALSE
  )
}
tab <- do.call(rbind, rows)
tab$q_value <- p.adjust(tab$p_value, method = "BH")
tab <- merge(tab, anno, by = "feature_id", all.x = TRUE, sort = FALSE)
tab$is_biomarker <- tab$q_value <= as.numeric(opt$alpha) &
  abs(tab$log2fc_enriched_vs_other) >= as.numeric(opt$`min-abs-log2fc`)
tab <- tab[order(tab$q_value, -abs(tab$log2fc_enriched_vs_other)), ]
mc_write_tsv(tab, file.path(dirs$tables, "biomarker_candidates.tsv"))

if (nlevels(group) == 2) {
  pw <- do.call(rbind, lapply(seq_len(ncol(rel)), function(j) {
    wt <- stats::wilcox.test(rel[, j] ~ group, exact = FALSE)
    data.frame(feature_id = colnames(rel)[j], group_a = levels_g[[1]], group_b = levels_g[[2]],
               p_value = wt$p.value, stringsAsFactors = FALSE)
  }))
  pw$q_value <- p.adjust(pw$p_value, method = "BH")
  mc_write_tsv(pw, file.path(dirs$tables, "pairwise_wilcoxon.tsv"))
} else {
  pw_rows <- list()
  for (j in seq_len(ncol(rel))) {
    pw <- stats::pairwise.wilcox.test(rel[, j], group, p.adjust.method = "BH", exact = FALSE)
    pmat <- as.data.frame(as.table(pw$p.value))
    names(pmat) <- c("group_a", "group_b", "p_adjust")
    pmat <- pmat[!is.na(pmat$p_adjust), ]
    if (nrow(pmat) > 0) {
      pmat$feature_id <- colnames(rel)[j]
      pw_rows[[length(pw_rows) + 1]] <- pmat[, c("feature_id", "group_a", "group_b", "p_adjust")]
    }
  }
  if (length(pw_rows) > 0) mc_write_tsv(do.call(rbind, pw_rows), file.path(dirs$tables, "pairwise_wilcoxon.tsv"))
}

top <- head(tab, as.integer(opt$`top-n`))
if (nrow(top) > 0 && mc_make_figures(opt$`fig-format`)) {
  top$label <- ifelse(is.na(top$taxon), top$feature_id, top$taxon)
  top$label <- factor(top$label, levels = rev(unique(top$label)))
  p_bar <- ggplot2::ggplot(top, ggplot2::aes(x = label, y = log2fc_enriched_vs_other, fill = enriched_group)) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = mc_palette(length(unique(top$enriched_group)))) +
    ggplot2::labs(x = NULL, y = "log2 fold change vs other groups", fill = opt$group,
                  title = "Top biomarker candidates") +
    mc_theme()
  mc_save_plot(p_bar, file.path(dirs$figures, "top_biomarker_effects"), 8, 6, opt$`fig-format`)

  features <- top$feature_id
  heat <- as.data.frame(clr[, features, drop = FALSE], check.names = FALSE)
  heat[[dat$sample_col]] <- rownames(clr)
  heat[[opt$group]] <- group
  heat_long <- reshape(heat, varying = features, v.names = "clr_abundance",
                       timevar = "feature_id", times = features, direction = "long")
  heat_long$feature_label <- ifelse(heat_long$feature_id %in% tab$feature_id,
                                    tab$taxon[match(heat_long$feature_id, tab$feature_id)],
                                    heat_long$feature_id)
  heat_long$feature_label[is.na(heat_long$feature_label)] <- heat_long$feature_id[is.na(heat_long$feature_label)]
  heat_long$feature_label <- factor(heat_long$feature_label, levels = rev(unique(top$label)))
  heat_long[[dat$sample_col]] <- factor(heat_long[[dat$sample_col]],
                                        levels = heat[[dat$sample_col]][order(group)])
  p_heat <- ggplot2::ggplot(heat_long, ggplot2::aes(x = .data[[dat$sample_col]], y = feature_label, fill = clr_abundance)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#D55E00", midpoint = 0) +
    ggplot2::labs(x = "Sample", y = NULL, fill = "CLR", title = "Top biomarker CLR heatmap") +
    mc_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7))
  mc_save_plot(p_heat, file.path(dirs$figures, "top_biomarker_heatmap"), 9, 7, opt$`fig-format`)
}

mc_log("Done. Output: ", opt$outdir, log_file = log_file)
