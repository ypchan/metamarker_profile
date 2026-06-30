# ============================================================
# Random forest importance + ridgeline abundance plot
# Nature-style template
#
# Purpose:
#   This script visualizes key taxa selected by Random Forest.
#   Left panel: variable importance bar plot.
#   Right panel: ridgeline/joyplot showing abundance patterns
#                of the selected taxa across groups or ordered samples.
#
# Required input files:
#
#   1. rf_importance.tsv
#      Random Forest variable importance table.
#
#      Required columns:
#        genus                    Character. Taxon name used as RF feature.
#        phylum                   Character. Phylum annotation for color.
#        mean_decrease_accuracy   Numeric. RF importance score.
#
#      Example:
#        genus                                 phylum             mean_decrease_accuracy
#        unclassified_Peptostreptococcales     Firmicutes         18.0
#        Permianibacter                        Proteobacteria     14.2
#        Xanthobacter                          Proteobacteria     9.8
#        Giesbergeria                          Proteobacteria     8.6
#        GCA_003244245                         Actinobacteriota   7.2
#
#      Notes:
#        - genus names must match the genus column in abundance_long.tsv.
#        - mean_decrease_accuracy must be numeric.
#        - duplicated genus names are not recommended.
#        - the top N taxa are selected according to mean_decrease_accuracy.
#
#
#   2. abundance_long.tsv
#      Long-format abundance table for plotting selected taxa.
#
#      Required columns:
#        sample_id    Character. Sample name.
#        group        Character. Group or stage name, e.g. D1, D3, D6.
#        x_value      Numeric. Ordered x-axis position.
#        genus        Character. Taxon name. Must match rf_importance.tsv.
#        phylum       Character. Phylum annotation.
#        abundance    Numeric. Taxon abundance.
#
#      Example:
#        sample_id    group   x_value   genus             phylum           abundance
#        S01          D1      1         Permianibacter    Proteobacteria   12.5
#        S02          D1      1         Permianibacter    Proteobacteria   15.1
#        S03          D3      2         Permianibacter    Proteobacteria   4.2
#        S04          D6      3         Permianibacter    Proteobacteria   8.8
#        S01          D1      1         Xanthobacter      Proteobacteria   3.1
#        S02          D3      2         Xanthobacter      Proteobacteria   20.2
#
#      Notes:
#        - sample_id does not need to be unique because this is long format.
#        - each row should represent one sample × genus observation.
#        - abundance can be RPM, TPM, relative abundance, read count,
#          CLR-transformed abundance, or other quantitative values.
#        - For sequencing count data, RPM/TPM or relative abundance is preferred.
#        - x_value controls the horizontal position of each sample/group.
#          For simple groups:
#              D1 = 1, D3 = 2, D6 = 3
#          For time-series data:
#              x_value can be day, month, year, sampling order, etc.
#
#
# Optional input format:
#
#   If your abundance table is wide format:
#
#        sample_id   group   x_value   Permianibacter   Xanthobacter   Giesbergeria
#        S01         D1      1         12.5             3.1            0.2
#        S02         D1      1         15.1             2.9            0.0
#        S03         D3      2         4.2              20.2           1.3
#
#   Convert it to long format before running this script:
#
#        abundance_long <- abundance_wide %>%
#          pivot_longer(
#            cols = -c(sample_id, group, x_value),
#            names_to = "genus",
#            values_to = "abundance"
#          ) %>%
#          left_join(taxon_annotation, by = "genus")
#
#   taxon_annotation should contain:
#
#        genus             phylum
#        Permianibacter    Proteobacteria
#        Xanthobacter      Proteobacteria
#        Giesbergeria      Proteobacteria
#
#
# Output files:
#
#   rf_ridgeline_plot/
#     rf_importance_ridgeline.pdf
#     rf_importance_ridgeline.png
#     top_genus_rf_importance.tsv
#     top_genus_abundance_plot_data.tsv
#
#
# Main parameters to modify:
#
#   abund_file  <- "abundance_long.tsv"
#   imp_file    <- "rf_importance.tsv"
#   top_n       <- 16
#
#   group levels:
#     group = factor(group, levels = c("D1", "D3", "D6"))
#
#   x-axis labels:
#     scale_x_continuous(
#       breaks = c(1, 2, 3),
#       labels = c("D1", "D3", "D6")
#     )
#
# ============================================================
# ============================================================
# Random forest importance + ridgeline abundance plot
# Nature-style template
# ============================================================

library(tidyverse)
library(ggridges)
library(patchwork)
library(scales)
library(grid)

# -----------------------------
# 1. Input files and parameters
# -----------------------------

abund_file <- "abundance_long.tsv"
imp_file   <- "rf_importance.tsv"

out_dir <- "rf_ridgeline_plot"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

top_n <- 16

# Group order shown on x-axis
group_levels <- c("D1", "D3", "D6")

# If TRUE, use log1p(abundance) for ridgeline height
use_log1p <- TRUE

# -----------------------------
# 2. Read input data
# -----------------------------

abund_raw <- readr::read_tsv(abund_file, show_col_types = FALSE)
imp_raw   <- readr::read_tsv(imp_file, show_col_types = FALSE)

# -----------------------------
# 3. Required columns
# -----------------------------

required_abund_cols <- c(
  "sample_id",
  "group",
  "x_value",
  "genus",
  "phylum",
  "abundance"
)

required_imp_cols <- c(
  "genus",
  "phylum",
  "mean_decrease_accuracy"
)

missing_abund_cols <- setdiff(required_abund_cols, colnames(abund_raw))
missing_imp_cols   <- setdiff(required_imp_cols, colnames(imp_raw))

if (length(missing_abund_cols) > 0) {
  stop(
    "abundance_long.tsv is missing required column(s): ",
    paste(missing_abund_cols, collapse = ", ")
  )
}

if (length(missing_imp_cols) > 0) {
  stop(
    "rf_importance.tsv is missing required column(s): ",
    paste(missing_imp_cols, collapse = ", ")
  )
}

# -----------------------------
# 4. Format checks
# -----------------------------

# Convert key columns
abund_raw <- abund_raw %>%
  mutate(
    sample_id = as.character(sample_id),
    group = as.character(group),
    x_value = as.numeric(x_value),
    genus = as.character(genus),
    phylum = as.character(phylum),
    abundance = as.numeric(abundance)
  )

imp_raw <- imp_raw %>%
  mutate(
    genus = as.character(genus),
    phylum = as.character(phylum),
    mean_decrease_accuracy = as.numeric(mean_decrease_accuracy)
  )

# Check numeric columns
if (any(is.na(abund_raw$x_value))) {
  stop("Column x_value in abundance_long.tsv contains NA or non-numeric values.")
}

if (any(is.na(abund_raw$abundance))) {
  stop("Column abundance in abundance_long.tsv contains NA or non-numeric values.")
}

if (any(is.na(imp_raw$mean_decrease_accuracy))) {
  stop("Column mean_decrease_accuracy in rf_importance.tsv contains NA or non-numeric values.")
}

# Check duplicated RF features
dup_genus <- imp_raw %>%
  count(genus) %>%
  filter(n > 1)

if (nrow(dup_genus) > 0) {
  warning(
    "Duplicated genus names found in rf_importance.tsv: ",
    paste(dup_genus$genus, collapse = ", "),
    "\nOnly the highest mean_decrease_accuracy will be kept for each genus."
  )

  imp_raw <- imp_raw %>%
    arrange(genus, desc(mean_decrease_accuracy)) %>%
    distinct(genus, .keep_all = TRUE)
}

# Check shared genus names
shared_genera <- intersect(imp_raw$genus, abund_raw$genus)

if (length(shared_genera) == 0) {
  stop("No shared genus names between rf_importance.tsv and abundance_long.tsv.")
}

message("[INFO] Number of RF genera: ", n_distinct(imp_raw$genus))
message("[INFO] Number of abundance genera: ", n_distinct(abund_raw$genus))
message("[INFO] Shared genera: ", length(shared_genera))

# -----------------------------
# 5. Select top important genera
# -----------------------------

top_genus_tbl <- imp_raw %>%
  filter(genus %in% shared_genera) %>%
  arrange(desc(mean_decrease_accuracy)) %>%
  slice_head(n = top_n) %>%
  mutate(
    genus = factor(genus, levels = rev(genus))
  )

top_genera <- as.character(top_genus_tbl$genus)

# -----------------------------
# 6. Prepare abundance table
# -----------------------------

plot_tbl <- abund_raw %>%
  filter(genus %in% top_genera) %>%
  left_join(
    top_genus_tbl %>%
      mutate(genus = as.character(genus)) %>%
      select(genus, mean_decrease_accuracy, phylum_rf = phylum),
    by = "genus"
  ) %>%
  mutate(
    phylum = coalesce(phylum, phylum_rf),
    group = factor(group, levels = group_levels),
    genus = factor(genus, levels = levels(top_genus_tbl$genus))
  )

if (use_log1p) {
  plot_tbl <- plot_tbl %>%
    mutate(abundance_plot = log1p(abundance))
} else {
  plot_tbl <- plot_tbl %>%
    mutate(abundance_plot = abundance)
}

# -----------------------------
# 7. Colors
# -----------------------------

pal_phylum <- c(
  "Proteobacteria"   = "#D7C84A",
  "Firmicutes"       = "#A8D08D",
  "Actinobacteriota" = "#A8AED3",
  "Bacteroidota"     = "#8ECFC9",
  "Chloroflexi"      = "#FFBE7A",
  "Acidobacteriota"  = "#FA7F6F",
  "Myxococcota"      = "#BEB8DC",
  "Others"           = "#BDBDBD"
)

plot_tbl <- plot_tbl %>%
  mutate(
    phylum_plot = if_else(phylum %in% names(pal_phylum), phylum, "Others")
  )

top_genus_tbl <- top_genus_tbl %>%
  mutate(
    phylum_plot = if_else(phylum %in% names(pal_phylum), phylum, "Others")
  )

group_bg <- tibble(
  group = factor(group_levels, levels = group_levels),
  xmin = seq_along(group_levels) - 0.5,
  xmax = seq_along(group_levels) + 0.5,
  fill = c("#F6D4D6", "#DDEBF2", "#E2F0E4")[seq_along(group_levels)]
)

# -----------------------------
# 8. Theme
# -----------------------------

theme_nature <- function(base_size = 10, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(color = "grey15"),
      axis.line = element_line(color = "grey20", linewidth = 0.35),
      axis.ticks = element_line(color = "grey20", linewidth = 0.3),
      axis.text = element_text(color = "grey20"),
      axis.title = element_text(color = "grey15"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(face = "bold"),
      legend.key = element_blank(),
      plot.margin = margin(4, 4, 4, 4)
    )
}

# -----------------------------
# 9. Left panel: RF importance bar
# -----------------------------

p_bar <- top_genus_tbl %>%
  ggplot(aes(
    x = mean_decrease_accuracy,
    y = genus,
    fill = phylum_plot
  )) +
  geom_col(width = 0.72, color = NA) +
  scale_x_reverse(
    name = "Mean Decrease Accuracy",
    position = "top",
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_fill_manual(
    values = pal_phylum,
    name = NULL,
    drop = FALSE
  ) +
  ylab(NULL) +
  theme_nature(base_size = 10) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    legend.position = c(0.06, 0.05),
    legend.justification = c(0, 0),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.32, "cm"),
    plot.margin = margin(4, 4, 4, 4)
  )

# -----------------------------
# 10. Middle panel: ridgeline abundance pattern
# -----------------------------

p_ridge <- ggplot() +
  geom_rect(
    data = group_bg,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = group),
    inherit.aes = FALSE,
    alpha = 0.65,
    color = NA
  ) +
  geom_ridgeline(
    data = plot_tbl,
    aes(
      x = x_value,
      y = genus,
      height = abundance_plot,
      group = genus
    ),
    scale = 0.85,
    fill = "#D9C1B2",
    color = "grey30",
    linewidth = 0.35,
    alpha = 0.85,
    min_height = 0
  ) +
  scale_fill_manual(
    values = setNames(group_bg$fill, group_bg$group),
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = seq_along(group_levels),
    labels = group_levels,
    position = "top",
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_discrete(position = "right") +
  labs(x = NULL, y = NULL) +
  theme_nature(base_size = 10) +
  theme(
    axis.text.x = element_text(face = "bold", size = 12),
    axis.text.y.right = element_text(
      size = 9,
      face = "italic",
      color = "grey15"
    ),
    axis.text.y.left = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(4, 4, 4, 4)
  )

# -----------------------------
# 11. Combine panels
# -----------------------------

p_final <- p_bar + p_ridge +
  plot_layout(widths = c(1.0, 1.75)) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(face = "bold", size = 16)
  )

print(p_final)

# -----------------------------
# 12. Save output
# -----------------------------

ggsave(
  file.path(out_dir, "rf_importance_ridgeline.pdf"),
  p_final,
  width = 8.5,
  height = 4.5,
  units = "in",
  device = cairo_pdf
)

ggsave(
  file.path(out_dir, "rf_importance_ridgeline.png"),
  p_final,
  width = 8.5,
  height = 4.5,
  units = "in",
  dpi = 600
)

readr::write_tsv(
  top_genus_tbl,
  file.path(out_dir, "top_genus_rf_importance.tsv")
)

readr::write_tsv(
  plot_tbl,
  file.path(out_dir, "top_genus_abundance_plot_data.tsv")
)

message("[DONE] ", file.path(out_dir, "rf_importance_ridgeline.pdf"))
message("[DONE] ", file.path(out_dir, "rf_importance_ridgeline.png"))
message("[DONE] ", file.path(out_dir, "top_genus_rf_importance.tsv"))
message("[DONE] ", file.path(out_dir, "top_genus_abundance_plot_data.tsv"))