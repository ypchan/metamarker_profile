# meta_marker_count

`meta_marker_count` counts marker-derived reads from paired-end clean reads. It uses BBDuk to extract candidate 16S/18S/ITS reads, minimap2 to align them to marker references, and a taxonomy table to produce marker RPM tables.

The default marker set is `16S,ITS`. Add 18S explicitly with `--markers 16S,18S,ITS`.

## Installation

Clone the repository and install the command wrappers:

```bash
git clone <repo-url>
cd meta_marker_count

mamba env create -f environment.yml
mamba activate meta-marker-count

make install PREFIX="$HOME/.local"
export PATH="$HOME/.local/bin:$PATH"
```

Installed commands:

- `meta_marker_count`: main pipeline.
- `meta_marker_build_refs`: converts SILVA/UNITE FASTA files into the default reference layout and minimap2 indexes.

Update an existing installation:

```bash
cd meta_marker_count
git pull --ff-only
make install PREFIX="$HOME/.local"
```

## Default Reference Directory

The software has built-in reference paths. By default, reference files are stored in:

```text
~/.local/share/meta_marker_count/ref
```

Override this directory with either:

```bash
export META_MARKER_COUNT_REF_DIR=/path/to/meta_marker_count_ref
```

or at runtime:

```bash
meta_marker_count --ref-dir /path/to/meta_marker_count_ref ...
```

After references are built in the default directory, normal pipeline runs do not need `--ref-16s`, `--index-16s`, `--ref-its`, `--taxonomy`, or related options.

Default files expected by `meta_marker_count`:

```text
~/.local/share/meta_marker_count/ref/SILVA_138.2_SSURef_tax_silva.dna.arc_bac.shortid.fasta
~/.local/share/meta_marker_count/ref/SILVA_138.2_SSURef_tax_silva.dna.arc_bac.shortid.fasta.mmi
~/.local/share/meta_marker_count/ref/SILVA_138.2_SSURef_tax_silva.dna.euk.shortid.fasta
~/.local/share/meta_marker_count/ref/SILVA_138.2_SSURef_tax_silva.dna.euk.shortid.fasta.mmi
~/.local/share/meta_marker_count/ref/UNITE_public_19.02.2025.shortid.fasta
~/.local/share/meta_marker_count/ref/UNITE_public_19.02.2025.shortid.fasta.mmi
~/.local/share/meta_marker_count/ref/ref_taxonomy.tsv
```

`meta_marker_count --help` also prints these default paths.

## Download Raw Reference Files

Download SILVA SSU directly:

```bash
curl -L \
  -o SILVA_138.2_SSURef_tax_silva.fasta.gz \
  "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.2_SSURef_tax_silva.fasta.gz"
```

Download UNITE ITS from the repository page:

```text
https://unite.ut.ee/repository.php
```

Select the public FASTA gzip release and save it as:

```text
UNITE_public_19.02.2025.fasta.gz
```

## Build And Configure References

Run this once after downloading the raw reference files:

```bash
meta_marker_build_refs \
  --silva SILVA_138.2_SSURef_tax_silva.fasta.gz \
  --unite UNITE_public_19.02.2025.fasta.gz \
  --force
```

By default this writes to `~/.local/share/meta_marker_count/ref` and builds `.mmi` indexes next to the generated FASTA files. The main pipeline will automatically use these paths.

If minimap2 indexes should be created later, use:

```bash
meta_marker_build_refs \
  --silva SILVA_138.2_SSURef_tax_silva.fasta.gz \
  --unite UNITE_public_19.02.2025.fasta.gz \
  --skip-index \
  --force
```

Then build indexes manually:

```bash
REF_DIR="${META_MARKER_COUNT_REF_DIR:-$HOME/.local/share/meta_marker_count/ref}"

minimap2 -d "$REF_DIR/SILVA_138.2_SSURef_tax_silva.dna.arc_bac.shortid.fasta.mmi" \
  "$REF_DIR/SILVA_138.2_SSURef_tax_silva.dna.arc_bac.shortid.fasta"

minimap2 -d "$REF_DIR/SILVA_138.2_SSURef_tax_silva.dna.euk.shortid.fasta.mmi" \
  "$REF_DIR/SILVA_138.2_SSURef_tax_silva.dna.euk.shortid.fasta"

minimap2 -d "$REF_DIR/UNITE_public_19.02.2025.shortid.fasta.mmi" \
  "$REF_DIR/UNITE_public_19.02.2025.shortid.fasta"
```

Generated reference files:

```text
UNITE_public_19.02.2025.shortid.fasta
UNITE_public_19.02.2025.taxonomy.tsv
SILVA_138.2_SSURef_tax_silva.dna.arc_bac.shortid.fasta
SILVA_138.2_SSURef_tax_silva.dna.arc_bac.taxonomy.tsv
SILVA_138.2_SSURef_tax_silva.dna.euk.shortid.fasta
SILVA_138.2_SSURef_tax_silva.dna.euk.taxonomy.tsv
ref_taxonomy.tsv
```

Combined taxonomy format:

```tsv
ref_id	marker	domain	phylum	class	order	family	genus	species
AB000393.1.1510	16S	Bacteria	Pseudomonadota	Gammaproteobacteria	Enterobacterales	Vibrionaceae	Vibrio	Vibrio_halioticoli
```

Reference parsing rules:

- UNITE ITS headers such as `UDB016649|k__Fungi;...|SH1281904.10FU` use `UDB016649` as `ref_id` and marker `ITS`.
- SILVA `Bacteria` and `Archaea` records are written to the arc/bac FASTA and marked as `16S`.
- SILVA `Eukaryota` records are written to the euk FASTA and marked as `18S`.
- SILVA eukaryotic lineages can have more than seven levels; the last two levels are retained as `genus` and `species`.
- SILVA records with an exact lineage item `Mitochondria` or `Chloroplast` are skipped by default. Use `--keep-organelles` to keep them.

## Sample Input

Multi-sample input is a TSV file with a header. Columns are matched by name, not by order.

Required columns:

- `sample_id`
- `r1_path`
- `r2_path`

Optional columns:

- `year`
- `month`
- `depth`
- any other sample metadata column

Example:

```tsv
sample_id	year	month	depth	site_type	r1_path	r2_path
202311_MF1_00-10	2023	11	00-10	MF	/data/202311_MF1_00-10_R1.fq.gz	/data/202311_MF1_00-10_R2.fq.gz
202311_MF1_10-20	2023	11	10-20	MF	/data/202311_MF1_10-20_R1.fq.gz	/data/202311_MF1_10-20_R2.fq.gz
```

## Run The Pipeline

Default 16S + ITS run:

```bash
meta_marker_count \
  --input data_path.tsv \
  --outdir marker_count_out \
  --markers 16S,ITS \
  --rank genus \
  --jobs 8 \
  --threads-per-sample 4
```

16S + 18S + ITS run:

```bash
meta_marker_count \
  --input data_path.tsv \
  --outdir marker_count_out \
  --markers 16S,18S,ITS \
  --rank genus \
  --jobs 8 \
  --threads-per-sample 4
```

Single-sample run:

```bash
meta_marker_count \
  --sample-id 202311_MF1_00-10 \
  --r1 /data/202311_MF1_00-10_R1.fq.gz \
  --r2 /data/202311_MF1_00-10_R2.fq.gz \
  --outdir marker_count_out
```

Dry run:

```bash
meta_marker_count \
  --input data_path.tsv \
  --outdir marker_count_out \
  --markers 16S,ITS \
  --dry-run
```

## Key Options

- `--steps all`: run all steps. Subsets are allowed: `reads_count,extract,align,abundance`.
- `--markers 16S,ITS`: marker list. Allowed markers are `16S`, `18S`, `ITS`.
- `--rank genus`: output rank. Allowed ranks are `domain,phylum,class,order,family,genus,species,all`.
- `--ref-dir DIR`: use a non-default reference directory.
- `--jobs INT`: samples processed in parallel.
- `--threads-per-sample INT`: BBDuk/minimap2 threads per sample.
- `--bbmap-mem 64G`: Java heap for BBDuk.
- `--force`: overwrite and rerun existing outputs.
- `--no-resume`: ignore checkpoints.
- `--clean-tmp`: remove temporary directories after a successful run.

## Output

Default output directory:

```text
marker_count_out/
  00_manifest/samples.tsv
  00_logs/marker_count.YYYYmmdd_HHMMSS.log
  00_logs/run_config.tsv
  01_reads_count/clean_reads.tsv
  02_marker_reads/
  03_align/
  04_abundance/
    marker_rpm.<rank>.long.tsv
    marker_rpm.<rank>.matrix.tsv
    marker_rpm.<rank>.domain_total.tsv
    marker_rpm.<rank>.assignment_stats.tsv
  99_checkpoints/
```

Important result tables:

- `marker_rpm.<rank>.long.tsv`: one row per sample/marker/domain/taxon.
- `marker_rpm.<rank>.matrix.tsv`: sample by feature matrix.
- `marker_rpm.<rank>.domain_total.tsv`: total marker reads/RPM per domain and marker.
- `marker_rpm.<rank>.assignment_stats.tsv`: alignment filtering and ambiguity statistics.

## R Analysis Scripts

The scripts in `scripts/` are RStudio templates, not command-line tools. Open a script in RStudio, edit the `MC_CONFIG <- list(...)` block near the top, then click **Source** or run the whole file.

Use the repository root as the RStudio working directory so relative paths such as `marker_count_out/...` resolve predictably.

Each template writes the same output layout:

```text
OUTDIR/
  tables/
  figures/
  logs/
```

Figures use the shared theme and palette from `scripts/analysis_common.R` and are saved as PDF and PNG by default. Set `fig-format = "none"` inside `MC_CONFIG` when only tables are needed.

### R Dependencies

The conda environment includes the required R packages:

```bash
mamba env update -f environment.yml
```

For a Windows R/RStudio installation, install the required CRAN packages in RStudio:

```r
install.packages(c("ggplot2", "vegan", "igraph"))
```

`igraph` is optional for the network layout figure; edge/node tables are still produced without it.

### Input Formats For R Templates

Most templates accept `sample_rows` input, matching `meta_marker_count` matrix output:

```tsv
sample_id	year	month	depth	Bacteria|16S|genus__Vibrio	Fungi|ITS|genus__Fusarium
S1	2023	11	00-10	12.5	3.2
S2	2023	11	10-20	0.0	8.1
```

Set non-feature columns in `metadata-cols`, for example `"year,month,depth,site_type"`.

Functional prediction tables often use `feature_rows`:

```tsv
function_id	S1	S2	S3
PWY-1234	10	4	8
PWY-5678	0	9	3
```

Long abundance tables are used by the time/rank/spatiotemporal template:

```tsv
sample_id	marker	domain	rank	taxon	taxon_marker_reads	clean_reads_total	marker_rpm
S1	16S	Bacteria	genus	Vibrio	20	1000000	20
```

### Template Groups

| Category | RStudio template | Main outputs | Default statistical approach |
| --- | --- | --- | --- |
| Diversity | `scripts/diversity_alpha_beta.R` | alpha diversity, PCoA, PERMANOVA, beta dispersion | Bray-Curtis on relative abundance, 9999 permutations |
| Differential abundance | `scripts/differential_abundance_compositional.R` | overall tests, group contrasts, volcano plot, top feature bars | CLR transform, linear model with optional covariates, BH-FDR |
| Biomarkers | `scripts/biomarker_discovery_nonparametric.R` | biomarker candidates, pairwise Wilcoxon tests, effect-size plot, CLR heatmap | Kruskal-Wallis, BH-FDR, minimum log2 fold-change |
| Correlation and network | `scripts/network_association_analysis.R` | all correlations, filtered edge table, node metrics, heatmap, optional network plot | Spearman correlation on CLR abundance, global BH-FDR |
| Functional prediction | `scripts/functional_profile_analysis.R` | top predicted functions, PCoA, PERMANOVA, CLR differential tests | Bray-Curtis PERMANOVA and CLR linear models |
| Other analyses | `scripts/taxa_time_rank_spacetime.R` | target taxa dynamics, rank diversity, time/space composition | long-table summaries, Shannon/Simpson by rank |
| Other analyses | `scripts/disturbance_event_response.R` | event PERMANOVA, environment response tests, focal-taxon models, lagged microbe screens | longitudinal event screening with site/time adjustment |

### Figure Style

All R templates use the same plotting style from `scripts/analysis_common.R`. The defaults are intentionally close to common Nature-family figure conventions:

- white background with no decorative grid
- thin black axes and ticks
- compact sans-serif text
- restrained color-blind-aware palette based on Okabe-Ito-style colors
- PDF and 320-dpi PNG output by default
- consistent output folders and file names for multi-panel figure assembly

The scripts produce analysis-ready figures, but final journal figures may still need panel lettering, exact physical sizing, and small manual adjustments in Illustrator, Inkscape, PowerPoint, or R.

### Statistical Defaults And Limitations

These templates prioritize defensible defaults over speed or decorative plots. They are suitable starting points for publication-grade analysis, but they do not replace study-design review, biological interpretation, or sensitivity analysis.

| Template | Default methods | Key parameters | Main limitations |
| --- | --- | --- | --- |
| `diversity_alpha_beta.R` | Observed richness, Shannon, Simpson, inverse Simpson, Pielou; Bray-Curtis PCoA; PERMANOVA with `vegan::adonis2`; beta dispersion with `betadisper` | `distance = "bray"`, `normalize = "relative"`, `permutations = 9999`, `min-prevalence = 0.10` | PERMANOVA can be significant because of dispersion differences, so inspect `beta_dispersion_*`. Bray-Curtis is abundance-sensitive and does not model compositional uncertainty. Rarefaction is not performed by default. |
| `differential_abundance_compositional.R` | Relative abundance followed by CLR transform; per-feature linear models; optional covariates; BH-FDR | `pseudocount = 1e-06`, `min-prevalence = 0.10`, `alpha = 0.05` | CLR models require enough samples per group and sensible pseudocount choice. Strong zero inflation, batch effects, and unbalanced designs need sensitivity checks. This is not a replacement for ANCOM-BC2, ALDEx2, MaAsLin2, or DESeq2 when those are required by a study. |
| `biomarker_discovery_nonparametric.R` | Kruskal-Wallis screening; pairwise Wilcoxon tests; BH-FDR; enriched-group log2 fold change | `alpha = 0.05`, `min-abs-log2fc = 1`, `min-prevalence = 0.10` | Nonparametric tests do not adjust for covariates in this template. Biomarker calls are association markers, not mechanistic proof. |
| `network_association_analysis.R` | CLR transform; pairwise Spearman correlation; global BH-FDR; edge filtering by FDR and absolute correlation; optional igraph layout | `transform = "clr"`, `cor-method = "spearman"`, `min-prevalence = 0.20`, `top-n = 100`, `min-abs-cor = 0.60`, `alpha = 0.05` | Correlation networks are not causal networks. Compositional data can induce indirect correlations. Results depend strongly on prevalence filtering, transformation, and sample size. Use networks as hypotheses, not proof. |
| `functional_profile_analysis.R` | Top function summaries; Bray-Curtis PCoA; PERMANOVA; CLR linear models for group contrasts | `permutations = 9999`, `min-prevalence = 0.10`, `pseudocount = 1e-06` | Functional prediction accuracy depends on the external tool and database. This script analyzes predicted tables; it does not infer functions from reads. Treat unvalidated predictions cautiously. |
| `taxa_time_rank_spacetime.R` | Target-taxon summaries; Shannon/Simpson by rank; top-taxon composition across time/space/group | `target-rank = "genus"`, `top-n = 20` | Primarily descriptive. It does not model temporal autocorrelation or repeated measures. Use mixed models or time-series models for formal inference when needed. |
| `disturbance_event_response.R` | Bray-Curtis PERMANOVA with site strata; event effects on environmental variables; focal-taxon linear models; PCoA-axis driver screening; lagged focal-microbe association screen | `permutations = 9999`, `min-prevalence = 0.10`, `top-n = 100`, `pseudocount = 1e-06` | This is evidence screening, not causal proof. Event, environment, focal fungi, bacteria, and archaea may be confounded by time, site, season, and unmeasured drivers. Lagged associations require enough repeated time points per site; small time series are unstable. |

Recommended reporting:

- report sample size per group/site/time point
- report feature filtering thresholds
- report distance metric, transformation, pseudocount, and permutation count
- report FDR method as Benjamini-Hochberg
- for PERMANOVA, report whether dispersion tests were checked
- for disturbance analysis, avoid causal wording unless supported by experimental design or external evidence

### Diversity Template

Open `scripts/diversity_alpha_beta.R` and edit:

```r
MC_CONFIG <- list(
  table = "marker_count_out/04_abundance/marker_rpm.genus.matrix.tsv",
  metadata = "",
  `table-format` = "sample_rows",
  `metadata-cols` = "year,month,depth,site_type",
  group = "site_type",
  covariates = "year,month",
  strata = "",
  outdir = "analysis_out/diversity_genus"
)
```

### Differential Abundance Template

Open `scripts/differential_abundance_compositional.R` and edit:

```r
MC_CONFIG <- list(
  table = "marker_count_out/04_abundance/marker_rpm.genus.matrix.tsv",
  metadata = "",
  `table-format` = "sample_rows",
  `metadata-cols` = "year,month,depth,site_type",
  group = "site_type",
  reference = "",
  covariates = "year,month",
  outdir = "analysis_out/differential_genus"
)
```

### Biomarker Template

Open `scripts/biomarker_discovery_nonparametric.R` and edit:

```r
MC_CONFIG <- list(
  table = "marker_count_out/04_abundance/marker_rpm.genus.matrix.tsv",
  metadata = "",
  `table-format` = "sample_rows",
  `metadata-cols` = "year,month,depth,site_type",
  group = "site_type",
  outdir = "analysis_out/biomarkers_genus"
)
```

### Correlation And Network Template

Open `scripts/network_association_analysis.R` and edit:

```r
MC_CONFIG <- list(
  table = "marker_count_out/04_abundance/marker_rpm.genus.matrix.tsv",
  metadata = "",
  `table-format` = "sample_rows",
  `metadata-cols` = "year,month,depth",
  outdir = "analysis_out/network_genus",
  transform = "clr",
  `cor-method` = "spearman"
)
```

### Functional Prediction Template

Use this after external tools such as PICRUSt2, Tax4Fun, or FUNGuild produce a function abundance table. Open `scripts/functional_profile_analysis.R` and edit:

```r
MC_CONFIG <- list(
  table = "picrust2_pathways.tsv",
  metadata = "metadata.tsv",
  `table-format` = "feature_rows",
  `feature-col` = "function_id",
  group = "site_type",
  covariates = "year,month",
  outdir = "analysis_out/functions_pathways"
)
```

### Target Taxa, Rank Diversity, And Spatiotemporal Template

Open `scripts/taxa_time_rank_spacetime.R` and edit:

```r
MC_CONFIG <- list(
  `long-table` = "marker_count_out/04_abundance/marker_rpm.genus.long.tsv",
  metadata = "metadata.tsv",
  taxa = "Vibrio,Fusarium",
  `time-col` = "month",
  `space-col` = "site",
  `group-col` = "site_type",
  `target-rank` = "genus",
  outdir = "analysis_out/taxa_time_space"
)
```

### Disturbance Or Restoration Event Response Template

Use `scripts/disturbance_event_response.R` for repeated sampling at one or more sites where a disturbance, restoration activity, or management action occurs at known time points. The template screens evidence in a defensible order:

1. Did the event change community composition?
2. Which environmental variables changed after the event?
3. Do those environmental variables explain the focal microbial response?
4. After accounting for time/site/event, do focal microbes show lagged associations with other microbial groups?

This is an evidence-screening workflow, not proof of causality. Strong causal claims still require experimental design, replication, and domain justification.

Metadata example:

```tsv
sample_id	site	year	event	pH	moisture	organic_matter
S1	A	2019	pre	6.1	20.1	3.4
S2	A	2020	pre	6.2	19.7	3.3
S3	A	2021	post	6.8	24.3	4.1
S4	A	2022	post	6.9	25.0	4.5
```

RStudio configuration:

```r
MC_CONFIG <- list(
  table = "marker_count_out/04_abundance/marker_rpm.genus.matrix.tsv",
  metadata = "metadata.tsv",
  `table-format` = "sample_rows",
  `metadata-cols` = "year,site,event,pH,moisture,organic_matter",
  `site-col` = "site",
  `time-col` = "year",
  `event-col` = "event",
  `env-cols` = "pH,moisture,organic_matter",
  `focal-taxa` = "Fusarium,Ascomycota",
  outdir = "analysis_out/disturbance_event"
)
```

Core output tables:

- `community_event_permanova.tsv`: community-level event/time/environment effects.
- `environment_event_tests.tsv`: environmental variables that changed after the event.
- `focal_taxa_models.tsv`: event and environment terms explaining focal taxa.
- `community_axis_driver_screen.tsv`: event, environment, and focal taxa associations with main community gradients.
- `microbe_lag_response_screen.tsv`: lagged focal-microbe associations with other microbial responders.

## Checks

Run syntax checks:

```bash
make check
```

Show command defaults:

```bash
meta_marker_count --print-defaults
meta_marker_count --help
```
