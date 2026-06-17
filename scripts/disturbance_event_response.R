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
  `site-col` = "site",
  `time-col` = "year",
  `event-col` = "",
  `event-time-col` = "",
  `env-cols` = "",
  `focal-taxa` = "",
  outdir = "analysis_out/disturbance_event"
)

opt <- mc_parse_args(list(
  table = "", metadata = "", `table-format` = "sample_rows", `sample-col` = "sample_id",
  `feature-col` = "feature_id", `metadata-cols` = "", `site-col` = "site",
  `time-col` = "year", `event-col` = "", `event-time-col` = "", `env-cols` = "",
  `focal-taxa` = "", outdir = "analysis_out/disturbance_event", `min-prevalence` = "0.10",
  `top-n` = "100", permutations = "9999", pseudocount = "1e-06", `fig-format` = "both"
))

mc_need_packages(c("vegan", "ggplot2"))
dirs <- mc_make_outdirs(opt$outdir)
log_file <- file.path(dirs$logs, "disturbance_event_response.log")

dat <- mc_read_feature_data(opt$table, opt$metadata, opt$`table-format`, opt$`sample-col`,
                            opt$`feature-col`, opt$`metadata-cols`)
mat <- mc_filter_features(dat$abund, as.numeric(opt$`min-prevalence`), 0)
if (ncol(mat) < 3) stop("At least three features are required after filtering.")
meta <- dat$metadata
for (v in c(opt$`site-col`, opt$`time-col`)) {
  if (!v %in% names(meta)) stop("Required metadata column not found: ", v)
}
env_cols <- mc_split(opt$`env-cols`)
missing_env <- setdiff(env_cols, names(meta))
if (length(missing_env) > 0) stop("Environmental column(s) not found: ", paste(missing_env, collapse = ", "))

meta$.__time <- suppressWarnings(as.numeric(meta[[opt$`time-col`]]))
if (any(is.na(meta$.__time))) stop("--time-col must be numeric or coercible to numeric for this template.")
if (opt$`event-col` != "") {
  if (!opt$`event-col` %in% names(meta)) stop("Event column not found: ", opt$`event-col`)
  meta$.__event <- meta[[opt$`event-col`]]
} else if (opt$`event-time-col` != "") {
  if (!opt$`event-time-col` %in% names(meta)) stop("Event-time column not found: ", opt$`event-time-col`)
  event_time <- suppressWarnings(as.numeric(meta[[opt$`event-time-col`]]))
  meta$.__event <- ifelse(meta$.__time >= event_time, "post", "pre")
} else {
  stop("Set either MC_CONFIG[['event-col']] or MC_CONFIG[['event-time-col']] before sourcing this template.")
}
meta$.__event <- factor(meta$.__event)
if (nlevels(meta$.__event) < 2) stop("Event column must contain at least two levels.")
meta$.__site <- factor(meta[[opt$`site-col`]])

rel <- mc_normalize(mat, "relative")
clr <- mc_normalize(mat, "clr", as.numeric(opt$pseudocount))
anno <- mc_feature_annotation(colnames(mat))

dist_obj <- vegan::vegdist(rel, method = "bray")
pcoa <- stats::cmdscale(dist_obj, k = 2, eig = TRUE)
eig <- pcoa$eig[pcoa$eig > 0]
var_exp <- if (length(eig) >= 2) round(100 * eig[1:2] / sum(eig), 1) else c(NA, NA)
ord <- data.frame(sample_id = rownames(pcoa$points), Axis1 = pcoa$points[, 1], Axis2 = pcoa$points[, 2],
                  stringsAsFactors = FALSE)
names(ord)[1] <- dat$sample_col
ord <- merge(meta, ord, by = dat$sample_col, sort = FALSE)
ord <- ord[match(rownames(rel), ord[[dat$sample_col]]), , drop = FALSE]
mc_write_tsv(ord, file.path(dirs$tables, "community_pcoa_bray.tsv"))

model_vars <- c(".__event", ".__time", env_cols)
form <- stats::as.formula(paste("dist_obj ~", paste(model_vars, collapse = " + ")))
strata_vec <- if (nlevels(meta$.__site) > 1) meta$.__site else NULL
ad <- vegan::adonis2(form, data = meta, permutations = as.integer(opt$permutations),
                     by = "margin", strata = strata_vec)
ad_df <- data.frame(term = rownames(as.data.frame(ad)), as.data.frame(ad), check.names = FALSE)
mc_write_tsv(ad_df, file.path(dirs$tables, "community_event_permanova.tsv"))

extract_terms <- function(fit, pattern, response) {
  cf <- as.data.frame(summary(fit)$coefficients)
  cf$term <- rownames(cf)
  cf <- cf[grepl(pattern, cf$term), , drop = FALSE]
  if (nrow(cf) == 0) return(NULL)
  data.frame(response = response, term = cf$term, estimate = cf$Estimate,
             std_error = cf$`Std. Error`, statistic = cf$`t value`,
             p_value = cf$`Pr(>|t|)`, stringsAsFactors = FALSE)
}

env_rows <- list()
for (env in env_cols) {
  y <- suppressWarnings(as.numeric(meta[[env]]))
  if (mean(!is.na(y)) < 0.8) next
  df <- data.frame(y = y, event = meta$.__event, time = meta$.__time, site = meta$.__site)
  fit <- stats::lm(y ~ event + time + site, data = df)
  env_rows[[length(env_rows) + 1]] <- extract_terms(fit, "^event", env)
}
if (length(env_rows) > 0) {
  env_tab <- do.call(rbind, env_rows)
  env_tab$q_value <- p.adjust(env_tab$p_value, method = "BH")
  mc_write_tsv(env_tab, file.path(dirs$tables, "environment_event_tests.tsv"))
}

patterns <- mc_split(opt$`focal-taxa`)
focal_df <- NULL
if (length(patterns) > 0) {
  hit <- Reduce(`|`, lapply(patterns, function(p) grepl(p, colnames(rel), ignore.case = TRUE)))
  if (any(hit)) {
    focal <- rowSums(rel[, hit, drop = FALSE])
    focal_df <- data.frame(sample_id = rownames(rel), focal_taxa_abundance = focal,
                           log_focal_taxa = log(focal + as.numeric(opt$pseudocount)),
                           stringsAsFactors = FALSE)
    names(focal_df)[1] <- dat$sample_col
    focal_df <- merge(meta, focal_df, by = dat$sample_col, sort = FALSE)
    focal_df <- focal_df[match(rownames(rel), focal_df[[dat$sample_col]]), , drop = FALSE]
    mc_write_tsv(focal_df, file.path(dirs$tables, "focal_taxa_abundance.tsv"))

    fit_focal <- stats::lm(stats::as.formula(paste("log_focal_taxa ~ .__event + .__time + .__site",
                                                   if (length(env_cols) > 0) paste("+", paste(env_cols, collapse = " + ")) else "")),
                           data = focal_df)
    focal_event <- extract_terms(fit_focal, "^.__event", "focal_taxa")
    focal_env <- if (length(env_cols) > 0) do.call(rbind, lapply(env_cols, function(e) extract_terms(fit_focal, paste0("^", e, "$"), "focal_taxa"))) else NULL
    focal_tests <- rbind(focal_event, focal_env)
    if (!is.null(focal_tests)) {
      focal_tests$q_value <- p.adjust(focal_tests$p_value, method = "BH")
      mc_write_tsv(focal_tests, file.path(dirs$tables, "focal_taxa_models.tsv"))
    }
  } else {
    mc_log("No focal taxa matched --focal-taxa patterns.", log_file = log_file)
  }
}

axis_df <- data.frame(Axis1 = ord$Axis1, Axis2 = ord$Axis2, meta)
if (!is.null(focal_df)) axis_df$log_focal_taxa <- focal_df$log_focal_taxa
axis_rows <- list()
axis_predictors <- c(".__event", ".__time", ".__site", env_cols, if (!is.null(focal_df)) "log_focal_taxa")
for (axis in c("Axis1", "Axis2")) {
  fit <- stats::lm(stats::as.formula(paste(axis, "~", paste(axis_predictors, collapse = " + "))), data = axis_df)
  axis_rows[[length(axis_rows) + 1]] <- data.frame(response = axis, term = rownames(summary(fit)$coefficients),
                                                   as.data.frame(summary(fit)$coefficients),
                                                   check.names = FALSE)
}
axis_tab <- do.call(rbind, axis_rows)
axis_tab$q_value <- p.adjust(axis_tab$`Pr(>|t|)`, method = "BH")
mc_write_tsv(axis_tab, file.path(dirs$tables, "community_axis_driver_screen.tsv"))

if (!is.null(focal_df)) {
  top_features <- mc_rank_columns(rel, as.integer(opt$`top-n`))
  lag_base <- data.frame(sample_id = rownames(rel), site = meta$.__site, time = meta$.__time,
                         event = meta$.__event, log_focal = focal_df$log_focal_taxa,
                         stringsAsFactors = FALSE)
  names(lag_base)[1] <- dat$sample_col
  lag_rows <- list()
  for (feat in setdiff(top_features, colnames(rel)[hit])) {
    d <- lag_base
    d$y <- log(rel[, feat] + as.numeric(opt$pseudocount))
    d <- d[order(d$site, d$time), ]
    d$lag_y <- ave(d$y, d$site, FUN = function(z) c(NA, head(z, -1)))
    d$lag_focal <- ave(d$log_focal, d$site, FUN = function(z) c(NA, head(z, -1)))
    keep <- stats::complete.cases(d[, c("y", "lag_y", "lag_focal", "event", "time", "site")])
    if (sum(keep) < 6) next
    fit <- stats::lm(y ~ lag_y + lag_focal + event + time + site, data = d[keep, ])
    cf <- as.data.frame(summary(fit)$coefficients)
    if ("lag_focal" %in% rownames(cf)) {
      lag_rows[[length(lag_rows) + 1]] <- data.frame(
        responder_feature = feat,
        lag_focal_estimate = cf["lag_focal", "Estimate"],
        lag_focal_p_value = cf["lag_focal", "Pr(>|t|)"],
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(lag_rows) > 0) {
    lag_tab <- do.call(rbind, lag_rows)
    lag_tab$q_value <- p.adjust(lag_tab$lag_focal_p_value, method = "BH")
    lag_tab <- merge(lag_tab, anno, by.x = "responder_feature", by.y = "feature_id", all.x = TRUE, sort = FALSE)
    lag_tab <- lag_tab[order(lag_tab$q_value, lag_tab$lag_focal_p_value), ]
    mc_write_tsv(lag_tab, file.path(dirs$tables, "microbe_lag_response_screen.tsv"))
  }
}

if (mc_make_figures(opt$`fig-format`)) {
  p_ord <- ggplot2::ggplot(ord, ggplot2::aes(x = Axis1, y = Axis2, color = .__event, shape = .__site)) +
    ggplot2::geom_point(size = 3, alpha = 0.9) +
    ggplot2::scale_color_manual(values = mc_palette(nlevels(meta$.__event))) +
    ggplot2::labs(x = paste0("PCoA1 (", var_exp[1], "%)"), y = paste0("PCoA2 (", var_exp[2], "%)"),
                  color = "Event", shape = opt$`site-col`, title = "Community response to disturbance") +
    mc_theme()
  mc_save_plot(p_ord, file.path(dirs$figures, "community_event_pcoa"), 7, 5.5, opt$`fig-format`)

  if (!is.null(focal_df)) {
    p_focal <- ggplot2::ggplot(focal_df, ggplot2::aes(x = .__time, y = focal_taxa_abundance,
                                                      color = .__event, group = .__site)) +
      ggplot2::geom_line(alpha = 0.7) +
      ggplot2::geom_point(size = 2.4) +
      ggplot2::facet_wrap(~ .__site) +
      ggplot2::scale_color_manual(values = mc_palette(nlevels(meta$.__event))) +
      ggplot2::labs(x = opt$`time-col`, y = "Focal taxa relative abundance",
                    color = "Event", title = "Focal taxa response") +
      mc_theme()
    mc_save_plot(p_focal, file.path(dirs$figures, "focal_taxa_event_trajectory"), 8, 5.5, opt$`fig-format`)
  }
}

mc_log("Done. Outputs are screening evidence, not causal proof. Output: ", opt$outdir, log_file = log_file)
