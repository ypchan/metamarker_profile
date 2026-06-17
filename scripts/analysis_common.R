#!/usr/bin/env Rscript

mc_parse_args <- function(defaults = list()) {
  if (!exists("MC_CONFIG", envir = parent.frame(), inherits = FALSE)) {
    stop("MC_CONFIG was not found. Edit the MC_CONFIG block near the top of this RStudio template.")
  }
  cfg <- get("MC_CONFIG", envir = parent.frame(), inherits = FALSE)
  if (!is.list(cfg)) stop("MC_CONFIG must be a named list.")
  utils::modifyList(defaults, cfg)
}

mc_bool <- function(x) {
  tolower(as.character(x)) %in% c("1", "true", "t", "yes", "y")
}

mc_split <- function(x) {
  if (is.null(x) || is.na(x) || x == "") character(0) else trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

mc_need_packages <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) > 0) {
    stop(
      "Missing R package(s): ", paste(miss, collapse = ", "),
      "\nInstall them before running this script. Example:\n",
      "install.packages(c(", paste(sprintf('\"%s\"', miss), collapse = ", "), "))",
      call. = FALSE
    )
  }
}

mc_read_table <- function(path) {
  if (is.null(path) || path == "") stop("Missing input path.")
  if (!file.exists(path)) stop("Input file not found: ", path)
  read.delim(path, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE, quote = "")
}

mc_write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "NA")
}

mc_make_outdirs <- function(outdir) {
  dirs <- file.path(outdir, c("tables", "figures", "logs"))
  invisible(vapply(dirs, dir.create, logical(1), recursive = TRUE, showWarnings = FALSE))
  list(tables = dirs[[1]], figures = dirs[[2]], logs = dirs[[3]])
}

mc_log <- function(..., log_file = NULL) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
  message(msg)
  if (!is.null(log_file)) cat(msg, "\n", file = log_file, append = TRUE)
}

mc_numeric_columns <- function(df, exclude = character(0)) {
  nms <- setdiff(names(df), exclude)
  keep <- vapply(df[nms], function(x) {
    y <- suppressWarnings(as.numeric(x))
    mean(!is.na(y)) >= 0.95
  }, logical(1))
  nms[keep]
}

mc_read_feature_data <- function(table, metadata = "", table_format = "sample_rows",
                                 sample_col = "sample_id", feature_col = "feature_id",
                                 metadata_cols = "") {
  x <- mc_read_table(table)
  meta_cols <- unique(c(sample_col, mc_split(metadata_cols)))

  if (table_format == "sample_rows") {
    if (!sample_col %in% names(x)) stop("Sample column not found in abundance table: ", sample_col)
    feature_cols <- setdiff(mc_numeric_columns(x, exclude = meta_cols), sample_col)
    if (length(feature_cols) == 0) stop("No numeric feature columns were detected.")
    abund <- as.matrix(data.frame(lapply(x[feature_cols], function(v) as.numeric(v)), check.names = FALSE))
    rownames(abund) <- x[[sample_col]]
    meta <- x[setdiff(names(x), feature_cols)]
  } else if (table_format == "feature_rows") {
    if (!feature_col %in% names(x)) stop("Feature column not found: ", feature_col)
    feature_ids <- x[[feature_col]]
    sample_cols <- setdiff(names(x), feature_col)
    abund <- t(as.matrix(data.frame(lapply(x[sample_cols], function(v) as.numeric(v)), check.names = FALSE)))
    colnames(abund) <- feature_ids
    rownames(abund) <- sample_cols
    meta <- data.frame(sample_id = rownames(abund), stringsAsFactors = FALSE)
    names(meta)[1] <- sample_col
  } else {
    stop("MC_CONFIG[['table-format']] must be sample_rows or feature_rows.")
  }

  if (!is.null(metadata) && metadata != "") {
    m <- mc_read_table(metadata)
    if (!sample_col %in% names(m)) stop("Sample column not found in metadata: ", sample_col)
    meta <- merge(meta, m, by = sample_col, all.x = TRUE, sort = FALSE)
    meta <- meta[match(rownames(abund), meta[[sample_col]]), , drop = FALSE]
  }

  storage.mode(abund) <- "numeric"
  abund[is.na(abund)] <- 0
  if (any(abund < 0)) stop("Abundance table contains negative values.")
  list(abund = abund, metadata = meta, sample_col = sample_col)
}

mc_filter_features <- function(mat, min_prevalence = 0.10, min_total = 0, min_mean = 0) {
  prevalence <- colMeans(mat > 0)
  total <- colSums(mat)
  mean_abund <- colMeans(mat)
  keep <- prevalence >= min_prevalence & total >= min_total & mean_abund >= min_mean
  mat[, keep, drop = FALSE]
}

mc_normalize <- function(mat, method = "relative", pseudocount = 1e-06) {
  method <- tolower(method)
  if (method == "none") return(mat)
  rs <- rowSums(mat)
  rs[rs == 0] <- 1
  rel <- sweep(mat, 1, rs, "/")
  if (method == "relative") return(rel)
  if (method == "cpm") return(rel * 1e6)
  if (method == "hellinger") return(sqrt(rel))
  if (method == "clr") {
    z <- rel + pseudocount
    lx <- log(z)
    return(sweep(lx, 1, rowMeans(lx), "-"))
  }
  stop("Unknown normalization method: ", method)
}

mc_feature_annotation <- function(features) {
  out <- data.frame(feature_id = features, domain = NA_character_, marker = NA_character_,
                    rank = NA_character_, taxon = features, stringsAsFactors = FALSE)
  has_pipe <- grepl("\\|", features)
  parts <- strsplit(features[has_pipe], "\\|")
  out$domain[has_pipe] <- vapply(parts, function(x) x[[1]], character(1))
  out$marker[has_pipe] <- vapply(parts, function(x) if (length(x) >= 2) x[[2]] else NA_character_, character(1))
  rank_taxon <- vapply(parts, function(x) if (length(x) >= 3) x[[3]] else NA_character_, character(1))
  out$rank[has_pipe] <- sub("__.*$", "", rank_taxon)
  out$taxon[has_pipe] <- sub("^[^_]+__", "", rank_taxon)
  out
}

mc_theme <- function(base_size = 9) {
  ggplot2::theme_classic(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      axis.line = ggplot2::element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.30, colour = "black"),
      axis.ticks.length = grid::unit(2, "pt"),
      axis.text = ggplot2::element_text(colour = "black", size = base_size * 0.9),
      axis.title = ggplot2::element_text(colour = "black", size = base_size),
      plot.title = ggplot2::element_text(face = "plain", hjust = 0, size = base_size * 1.05),
      legend.title = ggplot2::element_text(face = "plain", size = base_size * 0.9),
      legend.text = ggplot2::element_text(size = base_size * 0.85),
      legend.key = ggplot2::element_rect(fill = "white", colour = NA),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "plain", colour = "black", size = base_size * 0.9)
    )
}

mc_palette <- function(n) {
  base <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9",
            "#000000", "#999999", "#332288", "#88CCEE", "#44AA99", "#AA4499",
            "#882255", "#DDCC77", "#117733")
  if (n <= length(base)) base[seq_len(n)] else grDevices::hcl.colors(n, palette = "Dark 3", rev = FALSE)
}

mc_save_plot <- function(plot, file_base, width = 7, height = 5, formats = "both", dpi = 320) {
  if (tolower(formats) == "none") return(invisible(NULL))
  fmts <- if (formats == "both") c("pdf", "png") else mc_split(formats)
  for (fmt in fmts) {
    path <- paste0(file_base, ".", fmt)
    ggplot2::ggsave(path, plot, width = width, height = height, dpi = dpi, limitsize = FALSE)
  }
}

mc_make_figures <- function(formats) {
  tolower(formats) != "none"
}

mc_complete_cases <- function(meta, vars) {
  vars <- vars[vars != "" & !is.na(vars)]
  if (length(vars) == 0) return(rep(TRUE, nrow(meta)))
  stats::complete.cases(meta[, vars, drop = FALSE])
}

mc_rank_columns <- function(mat, top_n = 30) {
  ord <- order(colMeans(mat), decreasing = TRUE)
  colnames(mat)[ord[seq_len(min(top_n, length(ord)))]]
}
