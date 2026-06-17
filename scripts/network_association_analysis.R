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
  outdir = "analysis_out/network",
  transform = "clr",
  `cor-method` = "spearman"
)

opt <- mc_parse_args(list(
  table = "", metadata = "", `table-format` = "sample_rows", `sample-col` = "sample_id",
  `feature-col` = "feature_id", `metadata-cols` = "", outdir = "analysis_out/network",
  transform = "clr", `cor-method` = "spearman", `min-prevalence` = "0.20",
  `min-total` = "0", `top-n` = "100", `min-abs-cor` = "0.60", alpha = "0.05",
  pseudocount = "1e-06", `fig-format` = "both"
))

mc_need_packages(c("ggplot2"))
dirs <- mc_make_outdirs(opt$outdir)
log_file <- file.path(dirs$logs, "network_association_analysis.log")

dat <- mc_read_feature_data(opt$table, opt$metadata, opt$`table-format`, opt$`sample-col`,
                            opt$`feature-col`, opt$`metadata-cols`)
mat <- mc_filter_features(dat$abund, as.numeric(opt$`min-prevalence`), as.numeric(opt$`min-total`))
if (ncol(mat) > as.integer(opt$`top-n`)) {
  mat <- mat[, mc_rank_columns(mat, as.integer(opt$`top-n`)), drop = FALSE]
}
if (ncol(mat) < 3) stop("At least three features are required for network analysis.")
if (nrow(mat) < 5) stop("At least five samples are recommended for correlation analysis.")

x <- mc_normalize(mat, opt$transform, as.numeric(opt$pseudocount))
n <- ncol(x)
cor_mat <- stats::cor(x, method = opt$`cor-method`, use = "pairwise.complete.obs")

pairs <- utils::combn(seq_len(n), 2)
edge_rows <- vector("list", ncol(pairs))
for (k in seq_len(ncol(pairs))) {
  i <- pairs[1, k]
  j <- pairs[2, k]
  ct <- suppressWarnings(stats::cor.test(x[, i], x[, j], method = opt$`cor-method`, exact = FALSE))
  edge_rows[[k]] <- data.frame(
    source = colnames(x)[i],
    target = colnames(x)[j],
    correlation = unname(ct$estimate),
    p_value = ct$p.value,
    stringsAsFactors = FALSE
  )
}
edges_all <- do.call(rbind, edge_rows)
edges_all$q_value <- p.adjust(edges_all$p_value, method = "BH")
edges_all$abs_correlation <- abs(edges_all$correlation)
edges_all$sign <- ifelse(edges_all$correlation >= 0, "positive", "negative")
edges_all <- edges_all[order(edges_all$q_value, -edges_all$abs_correlation), ]
mc_write_tsv(edges_all, file.path(dirs$tables, "all_pairwise_associations.tsv"))

edges <- edges_all[edges_all$q_value <= as.numeric(opt$alpha) &
                     edges_all$abs_correlation >= as.numeric(opt$`min-abs-cor`), ]
mc_write_tsv(edges, file.path(dirs$tables, "network_edges.tsv"))

anno <- mc_feature_annotation(colnames(x))
nodes <- anno
nodes$degree <- 0L
nodes$strength <- 0
if (nrow(edges) > 0) {
  for (v in nodes$feature_id) {
    e <- edges[edges$source == v | edges$target == v, , drop = FALSE]
    nodes$degree[nodes$feature_id == v] <- nrow(e)
    nodes$strength[nodes$feature_id == v] <- sum(abs(e$correlation))
  }
}
nodes <- nodes[order(-nodes$degree, -nodes$strength), ]
mc_write_tsv(nodes, file.path(dirs$tables, "network_nodes.tsv"))

if (mc_make_figures(opt$`fig-format`)) {
  plot_features <- head(nodes$feature_id, min(50, nrow(nodes)))
  heat <- as.data.frame(as.table(cor_mat[plot_features, plot_features, drop = FALSE]))
  names(heat) <- c("feature_a", "feature_b", "correlation")
  p_heat <- ggplot2::ggplot(heat, ggplot2::aes(x = feature_a, y = feature_b, fill = correlation)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#D55E00", midpoint = 0) +
    ggplot2::labs(x = NULL, y = NULL, fill = "r", title = "Feature association heatmap") +
    mc_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
                   axis.text.y = ggplot2::element_text(size = 6))
  mc_save_plot(p_heat, file.path(dirs$figures, "association_heatmap"), 9, 8, opt$`fig-format`)

  top_nodes <- head(nodes[nodes$degree > 0, ], 30)
  if (nrow(top_nodes) > 0) {
    top_nodes$label <- ifelse(is.na(top_nodes$taxon), top_nodes$feature_id, top_nodes$taxon)
    top_nodes$label <- factor(top_nodes$label, levels = rev(top_nodes$label))
    p_deg <- ggplot2::ggplot(top_nodes, ggplot2::aes(x = label, y = degree, fill = marker)) +
      ggplot2::geom_col(width = 0.75) +
      ggplot2::coord_flip() +
      ggplot2::labs(x = NULL, y = "Degree", fill = "Marker", title = "Top network hubs") +
      ggplot2::scale_fill_manual(values = mc_palette(length(unique(top_nodes$marker)))) +
      mc_theme()
    mc_save_plot(p_deg, file.path(dirs$figures, "top_network_hubs"), 8, 6, opt$`fig-format`)
  }

  if (nrow(edges) > 0 && requireNamespace("igraph", quietly = TRUE)) {
    g <- igraph::graph_from_data_frame(edges[, c("source", "target", "correlation", "q_value", "sign")],
                                       directed = FALSE, vertices = nodes)
    igraph::V(g)$size <- 4 + 2 * log1p(igraph::degree(g))
    igraph::E(g)$color <- ifelse(igraph::E(g)$sign == "positive", "#D55E00", "#3B6FB6")
    igraph::E(g)$width <- 1 + 3 * abs(igraph::E(g)$correlation)
    layout <- igraph::layout_with_fr(g)
    fmts <- if (opt$`fig-format` == "both") c("pdf", "png") else mc_split(opt$`fig-format`)
    for (fmt in fmts) {
      file <- file.path(dirs$figures, paste0("association_network.", fmt))
      if (fmt == "pdf") grDevices::pdf(file, width = 8, height = 7) else grDevices::png(file, width = 2400, height = 2100, res = 300)
      plot(g, layout = layout, vertex.label = NA, vertex.color = "#F2B701",
           main = "Association network")
      grDevices::dev.off()
    }
  } else {
    mc_log("igraph is not installed or no edges passed filters; network plot skipped.", log_file = log_file)
  }
}

mc_log("Done. Output: ", opt$outdir, log_file = log_file)
