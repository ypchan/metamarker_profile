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
  `long-table` = "",
  metadata = "",
  taxa = "",
  `time-col` = "",
  `space-col` = "",
  `group-col` = "",
  `target-rank` = "genus",
  outdir = "analysis_out/taxa_time_rank_spacetime"
)

opt <- mc_parse_args(list(
  `long-table` = "", metadata = "", `sample-col` = "sample_id", `taxon-col` = "taxon",
  `rank-col` = "rank", `value-col` = "marker_rpm", `marker-col` = "marker",
  `domain-col` = "domain", `time-col` = "", `space-col` = "", `group-col` = "",
  taxa = "", `target-rank` = "genus", outdir = "analysis_out/taxa_time_rank_spacetime",
  `top-n` = "20", `fig-format` = "both"
))

mc_need_packages(c("vegan", "ggplot2"))
dirs <- mc_make_outdirs(opt$outdir)
log_file <- file.path(dirs$logs, "taxa_time_rank_spacetime.log")

x <- mc_read_table(opt$`long-table`)
required <- c(opt$`sample-col`, opt$`taxon-col`, opt$`rank-col`, opt$`value-col`)
missing <- setdiff(required, names(x))
if (length(missing) > 0) stop("Missing required long-table column(s): ", paste(missing, collapse = ", "))
x[[opt$`value-col`]] <- as.numeric(x[[opt$`value-col`]])
x[[opt$`value-col`]][is.na(x[[opt$`value-col`]])] <- 0

if (opt$metadata != "") {
  meta <- mc_read_table(opt$metadata)
  if (!opt$`sample-col` %in% names(meta)) stop("Sample column not found in metadata: ", opt$`sample-col`)
  x <- merge(x, meta, by = opt$`sample-col`, all.x = TRUE, sort = FALSE)
}

for (v in c(opt$`time-col`, opt$`space-col`, opt$`group-col`)) {
  if (v != "" && !v %in% names(x)) stop("Column not found: ", v)
}

mc_write_tsv(x, file.path(dirs$tables, "long_table_with_metadata.tsv"))

patterns <- mc_split(opt$taxa)
if (length(patterns) > 0) {
  hit <- Reduce(`|`, lapply(patterns, function(p) grepl(p, x[[opt$`taxon-col`]], ignore.case = TRUE)))
  target <- x[hit, , drop = FALSE]
  mc_write_tsv(target, file.path(dirs$tables, "target_taxa_long.tsv"))
  if (nrow(target) > 0) {
    group_vars <- c(opt$`taxon-col`, opt$`time-col`, opt$`space-col`, opt$`group-col`)
    group_vars <- group_vars[group_vars != ""]
    agg <- aggregate(target[[opt$`value-col`]], target[group_vars], mean)
    names(agg)[ncol(agg)] <- "mean_abundance"
    mc_write_tsv(agg, file.path(dirs$tables, "target_taxa_summary.tsv"))
    xvar <- if (opt$`time-col` != "") opt$`time-col` else opt$`taxon-col`
    p_taxa <- ggplot2::ggplot(agg, ggplot2::aes(x = .data[[xvar]], y = mean_abundance,
                                                color = .data[[opt$`taxon-col`]],
                                                group = .data[[opt$`taxon-col`]])) +
      ggplot2::geom_line(linewidth = 0.6, alpha = 0.85) +
      ggplot2::geom_point(size = 2, alpha = 0.9) +
      ggplot2::labs(x = xvar, y = paste("Mean", opt$`value-col`), color = "Taxon",
                    title = "Target taxa dynamics") +
      ggplot2::scale_color_manual(values = mc_palette(length(unique(agg[[opt$`taxon-col`]])))) +
      mc_theme()
    if (opt$`space-col` != "") p_taxa <- p_taxa + ggplot2::facet_wrap(stats::as.formula(paste("~", opt$`space-col`)))
    mc_save_plot(p_taxa, file.path(dirs$figures, "target_taxa_dynamics"), 8, 5, opt$`fig-format`)
  }
}

rank_rows <- list()
for (rk in unique(x[[opt$`rank-col`]])) {
  xr <- x[x[[opt$`rank-col`]] == rk, , drop = FALSE]
  wide <- stats::xtabs(stats::as.formula(paste(opt$`value-col`, "~", opt$`sample-col`, "+", opt$`taxon-col`)), data = xr)
  mat <- as.matrix(wide)
  rank_rows[[length(rank_rows) + 1]] <- data.frame(
    sample_id = rownames(mat),
    rank = rk,
    observed = rowSums(mat > 0),
    shannon = vegan::diversity(mat, index = "shannon"),
    simpson = vegan::diversity(mat, index = "simpson"),
    stringsAsFactors = FALSE
  )
}
rank_div <- do.call(rbind, rank_rows)
names(rank_div)[1] <- opt$`sample-col`
sample_meta_cols <- unique(c(opt$`sample-col`, opt$`time-col`, opt$`space-col`, opt$`group-col`))
sample_meta_cols <- sample_meta_cols[sample_meta_cols != "" & sample_meta_cols %in% names(x)]
sample_meta <- unique(x[sample_meta_cols])
rank_div <- merge(rank_div, sample_meta, by = opt$`sample-col`, all.x = TRUE, sort = FALSE)
mc_write_tsv(rank_div, file.path(dirs$tables, "rank_diversity.tsv"))

y_group <- if (opt$`group-col` != "") opt$`group-col` else "rank"
p_rank <- ggplot2::ggplot(rank_div, ggplot2::aes(x = rank, y = shannon, fill = .data[[y_group]])) +
  ggplot2::geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.75) +
  ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.12, height = 0), size = 1.6, alpha = 0.75) +
  ggplot2::labs(x = "Rank", y = "Shannon diversity", fill = y_group, title = "Diversity across taxonomic ranks") +
  ggplot2::scale_fill_manual(values = mc_palette(length(unique(rank_div[[y_group]])))) +
  mc_theme()
mc_save_plot(p_rank, file.path(dirs$figures, "rank_shannon_diversity"), 7.5, 5, opt$`fig-format`)

comp <- x[x[[opt$`rank-col`]] == opt$`target-rank`, , drop = FALSE]
if (nrow(comp) > 0) {
  top_taxa <- names(sort(tapply(comp[[opt$`value-col`]], comp[[opt$`taxon-col`]], sum), decreasing = TRUE))
  top_taxa <- head(top_taxa, as.integer(opt$`top-n`))
  comp$taxon_plot <- ifelse(comp[[opt$`taxon-col`]] %in% top_taxa, comp[[opt$`taxon-col`]], "Other")
  group_vars <- c("taxon_plot", opt$`time-col`, opt$`space-col`, opt$`group-col`)
  group_vars <- group_vars[group_vars != ""]
  comp_sum <- aggregate(comp[[opt$`value-col`]], comp[group_vars], mean)
  names(comp_sum)[ncol(comp_sum)] <- "mean_abundance"
  mc_write_tsv(comp_sum, file.path(dirs$tables, paste0("composition_", opt$`target-rank`, ".tsv")))
  xvar <- if (opt$`time-col` != "") opt$`time-col` else if (opt$`group-col` != "") opt$`group-col` else "taxon_plot"
  p_comp <- ggplot2::ggplot(comp_sum, ggplot2::aes(x = .data[[xvar]], y = mean_abundance, fill = taxon_plot)) +
    ggplot2::geom_col(position = "fill", width = 0.75) +
    ggplot2::labs(x = xvar, y = "Relative contribution", fill = "Taxon",
                  title = paste("Top", opt$`target-rank`, "composition")) +
    ggplot2::scale_fill_manual(values = mc_palette(length(unique(comp_sum$taxon_plot)))) +
    mc_theme()
  if (opt$`space-col` != "") p_comp <- p_comp + ggplot2::facet_wrap(stats::as.formula(paste("~", opt$`space-col`)))
  mc_save_plot(p_comp, file.path(dirs$figures, paste0("composition_", opt$`target-rank`)), 9, 5.5, opt$`fig-format`)
}

mc_log("Done. Output: ", opt$outdir, log_file = log_file)
