# meta_marker_count

<p align="center">
  <b>Fast marker-read quantification for metagenomic paired-end clean reads</b><br>
  16S / 18S / ITS candidate-read extraction, reference alignment, taxonomy-aware assignment, and marker RPM/RPKM tables.
</p>

<p align="center">
  <img src="workflow.png" alt="meta_marker_count workflow" width="920">
</p>

<p align="center">
  <a href="#installation">Installation</a> ·
  <a href="#reference-database">Reference database</a> ·
  <a href="#input-format">Input format</a> ·
  <a href="#run-the-pipeline">Run</a> ·
  <a href="#output-layout">Output</a> ·
  <a href="#r-analysis-templates">R analysis</a> ·
  <a href="#license">License</a>
</p>

---

## What is `meta_marker_count`?

`meta_marker_count` is a Bash-based bioinformatics workflow for estimating marker-derived reads from paired-end metagenomic clean reads. It is designed for ecological comparison of bacterial, archaeal, and fungal marker signals across many samples.

The pipeline:

1. counts clean reads from paired FASTQ files;
2. extracts candidate marker reads with **BBDuk**;
3. aligns extracted reads to marker references with **minimap2**;
4. assigns reads to taxonomy using a combined reference taxonomy table;
5. reports marker abundance as **RPM** and **RPKM-like length-normalized marker abundance**.

Default markers:

```text
16S,ITS
```

Add eukaryotic SSU explicitly when needed:

```bash
--markers 16S,18S,ITS
```

### Marker meaning

| Marker | Primary target | Default reference source | Typical use |
|---|---|---|---|
| `16S` | Bacteria + Archaea SSU rRNA | SILVA SSURef NR99 | prokaryotic community signal |
| `18S` | Eukaryotic SSU rRNA | SILVA SSURef NR99 | broad eukaryotic marker signal |
| `ITS` | Fungal ITS | UNITE public release | fungal community signal |

---

## Repository layout

Recommended repository structure:

```text
meta_marker_count/
├── meta_marker_count.sh              # main pipeline
├── scripts/
│   ├── build_ref_files.py            # SILVA/UNITE reference builder
│   ├── analysis_common.R             # shared R plotting/stat helpers
│   ├── diversity_alpha_beta.R
│   ├── differential_abundance_compositional.R
│   ├── biomarker_discovery_nonparametric.R
│   ├── network_association_analysis.R
│   ├── functional_profile_analysis.R
│   ├── taxa_time_rank_spacetime.R
│   └── disturbance_event_response.R
├── refs/                             # generated reference FASTA/index/taxonomy files
├── docs/
│   └── meta_marker_count_workflow.svg
├── Makefile
├── README.md
└── LICENSE
```

---

## Installation

### 1. Clone the repository

```bash
gh repo clone ypchan/meta_marker_count
cd meta_marker_count
```

or:

```bash
git clone https://github.com/ypchan/meta_marker_count.git
cd meta_marker_count
```

### 2. Install command links into the user bin directory

The recommended installation mode creates symbolic links in `~/bin` and keeps the real scripts inside the repository. This makes updates simple: after `git pull`, the linked commands immediately use the updated scripts.

```bash
make install
```

By default, references are expected under the cloned repository:

```text
<repo>/refs
```

`make install` also writes `<repo>/.meta_marker_count_ref_dir` so both `meta_marker_count` and `meta_marker_build_refs` resolve the same reference directory without editing source code.

Expected links:

```text
~/bin/meta_marker_count       -> <repo>/meta_marker_count.sh
~/bin/meta_marker_build_refs  -> <repo>/scripts/build_ref_files.py
```

Check them:

```bash
ls -l ~/bin/meta_marker_count ~/bin/meta_marker_build_refs
which meta_marker_count
meta_marker_count --help
```

If `~/bin` is not in your `PATH`, add it once:

```bash
echo 'export PATH="$HOME/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

### 3. Syntax check

```bash
make check
```

This checks the Bash pipeline, the Python reference builder, and optionally R templates if `Rscript` is available.

---

## Update

Because installation uses symbolic links, updating is usually just:

```bash
cd meta_marker_count
git pull --ff-only
make check
```

If the repository was moved to a new directory, rerun:

```bash
make install
```

This refreshes links and the local reference-directory config.

---

## Dependencies

### Required runtime tools

| Tool | Used in | Purpose |
|---|---:|---|
| `bash` | all steps | workflow driver |
| `awk` | abundance | PAF parsing and table generation |
| `sort`, `find`, `xargs` | all steps | file handling and parallel task dispatch |
| `seqkit` | Step 01 | clean-read counting |
| `bbduk.sh` | Step 02 | marker candidate-read extraction |
| `minimap2` | Step 03 | alignment to marker reference indexes |
| `python3` | reference setup | SILVA/UNITE reference conversion |

### Optional tools

| Tool | Purpose |
|---|---|
| `Rscript` | syntax checks for R analysis templates |
| `mamba` / `conda` | environment management |
| `aria2c` | faster reference download |

### Example conda environment

```bash
mamba create -n meta_marker_count -c conda-forge -c bioconda \
  seqkit minimap2 bbmap python r-base r-ggplot2 r-vegan r-igraph

mamba activate meta_marker_count
```

---

## Reference database

`meta_marker_count` expects a compact reference directory containing processed FASTA files, minimap2 indexes, and a combined taxonomy table.

Default reference directory after installation:

```text
<repo>/refs
```

Reference directory resolution order:

1. `META_MARKER_COUNT_REF_DIR`
2. `<repo>/.meta_marker_count_ref_dir`
3. `<repo>/refs`

You can always override it at runtime:

```bash
meta_marker_count --ref-dir /path/to/refs ...
```

### Expected reference files

```text
refs/
├── SILVA_138.2_SSURef_NR99_tax_silva.dna.arc_bac.shortid.fasta
├── SILVA_138.2_SSURef_NR99_tax_silva.dna.arc_bac.shortid.fasta.mmi
├── SILVA_138.2_SSURef_NR99_tax_silva.dna.euk.shortid.fasta
├── SILVA_138.2_SSURef_NR99_tax_silva.dna.euk.shortid.fasta.mmi
├── UNITE_public_19.02.2025.shortid.fasta
├── UNITE_public_19.02.2025.shortid.fasta.mmi
└── ref_taxonomy.tsv
```

### Download raw SILVA SSU

```bash
mkdir -p raw_refs
cd raw_refs

aria2c \
  -c \
  -x 16 \
  -s 16 \
  -k 1M \
  --file-allocation=none \
  -o SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz \
  "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz"

gzip -t SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz
```

If the SILVA FASTA contains RNA bases, convert `U/u` to `T/t` safely without overwriting the input in-place:

```bash
zcat SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz \
  | awk '
      /^>/ { print; next }
      { gsub(/U/, "T"); gsub(/u/, "t"); print }
    ' \
  | gzip -c > SILVA_138.2_SSURef_NR99_tax_silva.dna.fasta.gz
```

### Download raw UNITE ITS

Download the public FASTA gzip release from the UNITE repository page and save it as:

```text
UNITE_public_19.02.2025.fasta.gz
```

### Build processed references

From the repository root:

```bash
meta_marker_build_refs \
  --silva raw_refs/SILVA_138.2_SSURef_NR99_tax_silva.dna.fasta.gz \
  --unite raw_refs/UNITE_public_19.02.2025.fasta.gz \
  --force
```

The builder writes processed FASTA files, taxonomy files, the combined `ref_taxonomy.tsv`, and minimap2 `.mmi` indexes into the configured reference directory.

### Build FASTA first, indexes later

Useful on clusters where index building should be submitted as a separate job:

```bash
meta_marker_build_refs \
  --silva raw_refs/SILVA_138.2_SSURef_NR99_tax_silva.dna.fasta.gz \
  --unite raw_refs/UNITE_public_19.02.2025.fasta.gz \
  --skip-index \
  --force
```

Then:

```bash
REF_DIR="refs"

minimap2 -d "$REF_DIR/SILVA_138.2_SSURef_NR99_tax_silva.dna.arc_bac.shortid.fasta.mmi" \
  "$REF_DIR/SILVA_138.2_SSURef_NR99_tax_silva.dna.arc_bac.shortid.fasta"

minimap2 -d "$REF_DIR/SILVA_138.2_SSURef_NR99_tax_silva.dna.euk.shortid.fasta.mmi" \
  "$REF_DIR/SILVA_138.2_SSURef_NR99_tax_silva.dna.euk.shortid.fasta"

minimap2 -d "$REF_DIR/UNITE_public_19.02.2025.shortid.fasta.mmi" \
  "$REF_DIR/UNITE_public_19.02.2025.shortid.fasta"
```

### Combined taxonomy format

`ref_taxonomy.tsv` must contain:

```tsv
ref_id	marker	domain	phylum	class	order	family	genus	species
AB000393.1.1510	16S	Bacteria	Pseudomonadota	Gammaproteobacteria	Enterobacterales	Vibrionaceae	Vibrio	Vibrio_halioticoli
UDB016649	ITS	Fungi	Ascomycota	Sordariomycetes	Hypocreales	Nectriaceae	Fusarium	Fusarium_sp
```

Reference parsing rules:

- SILVA `Bacteria` and `Archaea` records are treated as `16S`.
- SILVA `Eukaryota` records are treated as `18S`.
- UNITE records are treated as `ITS`.
- SILVA records with exact lineage item `Mitochondria` or `Chloroplast` are removed by default.
- Use the reference builder option `--keep-organelles` only when organellar SSU reads are part of the intended analysis.

---

## Input format

Multi-sample mode uses a tab-delimited file with a header. Columns are matched by name, not by position.

### Required columns

| Column | Meaning |
|---|---|
| `sample_id` | unique sample ID; avoid spaces and `/` |
| `r1_path` | path to clean R1 FASTQ or FASTQ.gz |
| `r2_path` | path to clean R2 FASTQ or FASTQ.gz |

### Optional metadata columns

`year`, `month`, `depth`, `site`, `site_type`, `group`, `habitat`, `replicate`, or any other metadata columns are allowed. Extra metadata columns are carried into long abundance outputs.

Example:

```tsv
sample_id	year	month	depth	site_type	r1_path	r2_path
202311_MF1_00-10	2023	11	00-10	MF	/data/202311_MF1_00-10_R1.fq.gz	/data/202311_MF1_00-10_R2.fq.gz
202311_MF1_10-20	2023	11	10-20	MF	/data/202311_MF1_10-20_R1.fq.gz	/data/202311_MF1_10-20_R2.fq.gz
```

Check the parsed input format:

```bash
meta_marker_count --print-format
```

---

## Run the pipeline

### Default 16S + ITS run

```bash
meta_marker_count \
  --input data_path.tsv \
  --outdir marker_count_out \
  --markers 16S,ITS \
  --rank genus \
  --jobs 4 \
  --threads-per-sample 10
```

### 16S + 18S + ITS run

```bash
meta_marker_count \
  --input data_path.tsv \
  --outdir marker_count_out \
  --markers 16S,18S,ITS \
  --rank genus \
  --jobs 4 \
  --threads-per-sample 10
```

### Single-sample run

```bash
meta_marker_count \
  --sample-id 202311_MF1_00-10 \
  --r1 /data/202311_MF1_00-10_R1.fq.gz \
  --r2 /data/202311_MF1_00-10_R2.fq.gz \
  --outdir marker_count_out
```

### Dry run

Use this before long cluster runs:

```bash
meta_marker_count \
  --input data_path.tsv \
  --outdir marker_count_out \
  --markers 16S,ITS \
  --dry-run
```

### Run only selected steps

```bash
# only count reads
meta_marker_count --input data_path.tsv --outdir marker_count_out --steps reads_count

# extract + align only
meta_marker_count --input data_path.tsv --outdir marker_count_out --steps extract,align

# recompute abundance from existing PAF files
meta_marker_count --input data_path.tsv --outdir marker_count_out --steps abundance --rank genus
```

---

## Step-by-step parameters

### Step 01 — counting reads

Counts reads from clean FASTQ files with `seqkit` and writes `reads_stat.tsv`.

| Option | Default | Meaning |
|---|---:|---|
| `--seqkit-threads` | `auto` | threads used by `seqkit stats` per sample |
| `--seqkit-max-threads` | `4` | cap for seqkit threads because gzip/I/O usually bottlenecks first |
| `--jobs` | `4` | number of samples processed in parallel |

Output:

```text
reads_stat.tsv
```

Core columns:

```text
sample_id  year  month  depth  r1_path  r2_path  read_pairs  clean_reads_total
```

Normalization denominator:

```text
clean_reads_total = R1 read count × 2
```

### Step 02 — extract marker candidate reads

Uses `BBDuk` to extract paired reads matching marker references.

| Option | 16S default | 18S default | ITS default | Meaning |
|---|---:|---:|---:|---|
| `--k-*` | `31` | `31` | `25` | k-mer length |
| `--hdist-*` | `0` | `0` | `0` | Hamming distance for k-mer matching |
| `--mkh-*` | `1` | `1` | `1` | minimum k-mer hits |
| `--mink-*` | `0` | `0` | `0` | shorter terminal k-mers; `0` disables |
| `--bbmap-mem` | colspan | colspan | `64G` | Java heap passed to BBDuk |
| `--bbduk-threads` | colspan | colspan | `auto` | BBDuk threads per sample |
| `--bbduk-max-threads` | colspan | colspan | `12` | cap for BBDuk threads |

Output:

```text
02_marker_reads/
├── 16S/<sample>.16S.R1.fq.gz
├── 16S/<sample>.16S.R2.fq.gz
├── ITS/<sample>.ITS.R1.fq.gz
├── ITS/<sample>.ITS.R2.fq.gz
└── stats/<sample>.<marker>.bbduk.stats.txt
```

### Step 03 — align marker reads

Uses `minimap2` to align extracted marker reads to marker indexes and writes PAF files.

| Option | Default | Meaning |
|---|---:|---|
| `--minimap2-preset` | `sr` | minimap2 preset for short reads |
| `--minimap2-N` | `10` | retain up to N secondary alignments |
| `--minimap2-threads` | `auto` | minimap2 threads per sample |
| `--minimap2-max-threads` | `12` | cap for minimap2 threads |

Output:

```text
03_align/
├── 16S/<sample>.16S.R1.paf
├── 16S/<sample>.16S.R2.paf
├── ITS/<sample>.ITS.R1.paf
└── ITS/<sample>.ITS.R2.paf
```

### Step 04 — assign taxonomy and calculate abundance

Parses PAF files, filters alignments, assigns each read pair to the single best PAF hit, and writes long/stat tables.

| Option | 16S default | 18S default | ITS default | Meaning |
|---|---:|---:|---:|---|
| `--min-identity-*` | `0.97` | `0.97` | `0.95` | minimum alignment identity |
| `--min-aln-len-*` | `80` | `80` | `80` | minimum aligned length |
| `--min-qcov-*` | `0.60` | `0.60` | `0.60` | minimum query coverage |
| `--min-mapq` | colspan | colspan | `0` | minimum MAPQ |
| `--top-identity-diff` | colspan | colspan | `0.010` | legacy compatibility option |
| `--top-aln-len-diff` | colspan | colspan | `10` | legacy compatibility option |
| `--rank` | colspan | colspan | `genus` | `domain,phylum,class,order,family,genus,species,all` |

Best-hit assignment order:

1. identity
2. aligned length
3. matching bases
4. MAPQ
5. query coverage
6. target length
7. mate (`R1` before `R2` only as a deterministic final biological tie-break)
8. target ID

PE handling is read-level counting with pair-level arbitration:

- if both mates support the same lineage, both reads are counted for that lineage;
- if only one mate has a valid hit, that one read is counted;
- if mates support different lineages, both reads are assigned to the better mate classification;
- discordant pairs are counted in `discordant_pairs`.
- reads whose own mate-level classification was overridden are counted in `reassigned_discordant_reads`.

Abundance formulas:

```text
marker_rpm  = taxon_marker_reads / clean_reads_total × 1,000,000
marker_rpkm = taxon_marker_reads / (clean_reads_total / 1,000,000) / (mean_target_length_bp / 1,000)
```

Output:

```text
all.marker_rpm.<rank>.long.tsv
all.marker_rpm.<rank>.assignment_stats.tsv
```

---

## Main command options

### Input and workflow

| Option | Default | Description |
|---|---:|---|
| `--input FILE` | none | multi-sample TSV input |
| `--sample-id ID` | none | single-sample ID |
| `--r1 FILE` | none | single-sample R1 FASTQ |
| `--r2 FILE` | none | single-sample R2 FASTQ |
| `--outdir DIR` | `marker_count_out` | output directory |
| `--steps LIST` | `all` | `reads_count,extract,align,abundance` or `all` |
| `--markers LIST` | `16S,ITS` | marker set: `16S`, `18S`, `ITS` |
| `--rank RANK` | `genus` | output taxonomy rank or `all` |

### Parallelism and HPC resources

| Option | Default | Description |
|---|---:|---|
| `--jobs INT` | `4` | samples processed in parallel |
| `--threads-per-sample INT` | `4` | per-sample CPU budget before tool caps |
| `--seqkit-threads INT\|auto` | `auto` | seqkit threads per sample |
| `--seqkit-max-threads INT` | `4` | cap for seqkit threads |
| `--bbduk-threads INT\|auto` | `auto` | BBDuk threads per sample |
| `--bbduk-max-threads INT` | `12` | cap for BBDuk threads |
| `--minimap2-threads INT\|auto` | `auto` | minimap2 threads per sample |
| `--minimap2-max-threads INT` | `12` | cap for minimap2 threads |
| `--bbmap-mem STR` | `64G` | BBDuk Java heap |

Thread scheduling:

```text
CPU budget = jobs × threads-per-sample
effective tool threads = min(tool request, tool cap)
step sample parallelism = floor(CPU budget / effective tool threads)
```

Example for a 40-core node:

```bash
meta_marker_count \
  --input data_path.tsv \
  --outdir marker_count_out \
  --jobs 4 \
  --threads-per-sample 10 \
  --bbmap-mem 80G
```

With this example, seqkit uses up to 4 threads, while BBDuk/minimap2 use up to 10 threads. If `--threads-per-sample 24` is used, BBDuk/minimap2 are capped at 12 and the remaining CPU budget is used for more sample-level parallelism.

### Resume, overwrite, and logs

| Option | Meaning |
|---|---|
| `--force` | rerun and overwrite existing outputs |
| `--no-resume` | ignore checkpoints |
| `--clean-tmp` | remove temporary files after success |
| `--dry-run` | validate inputs and print run plan only |
| `--check-deps` | check external tools and exit |
| `--keep-task-logs` | keep per-sample logs even when tasks succeed |
| `--progress-interval INT` | seconds between progress refreshes; default `5` |

---

## Output layout

The current output layout keeps important result files directly under `OUTDIR` and stores large intermediate files in numbered directories.

```text
marker_count_out/
├── sample_manifest.tsv
├── run_config.tsv
├── meta_marker_count.YYYYmmdd_HHMMSS.log
├── commands/
│   ├── reads_count/
│   ├── extract/
│   ├── align/
│   └── abundance/
├── reads_stat.tsv
├── all.marker_rpm.genus.long.tsv
├── all.marker_rpm.genus.assignment_stats.tsv
├── 02_marker_reads/
│   ├── 16S/
│   ├── 18S/
│   ├── ITS/
│   ├── stats/
│   └── tmp/
├── 03_align/
│   ├── 16S/
│   ├── 18S/
│   ├── ITS/
│   └── tmp/
├── .checkpoints/
└── .tmp/
```

### Important output tables

| File | Description |
|---|---|
| `sample_manifest.tsv` | parsed and standardized sample manifest |
| `run_config.tsv` | run-time parameters and resolved reference paths |
| `commands/` | re-runnable command lines grouped by stage and sample |
| `reads_stat.tsv` | read-pair counts and clean-read totals |
| `all.marker_rpm.<rank>.long.tsv` | tidy long abundance table |
| `all.marker_rpm.<rank>.assignment_stats.tsv` | alignment filtering, paired-end discordance, reassigned discordant reads, tie-breaking, and assignment statistics |

### Long-table columns

```text
sample_id
marker
domain
rank
lineage
taxon_marker_reads
clean_reads_total
marker_rpm
marker_rpkm
mean_target_length_bp
<additional sample metadata columns>
```

---

## Logging and progress

The screen output is designed for HPC runs:

- dependency checks are printed once at the beginning;
- long-running steps use progress bars;
- re-runnable command lines are always preserved under `commands/`;
- successful task logs are removed by default;
- failed task logs are preserved under `.tmp/task_logs/`;
- the main log is kept as `meta_marker_count.YYYYmmdd_HHMMSS.log`.

Example:

```text
[2026-06-18 17:24:12] [RUN]  counting reads  [======------------------]  25%  16/65  elapsed=3m10s
[2026-06-18 17:35:43] [DONE] counting reads  [========================] 100%  65/65  elapsed=11m31s
[2026-06-18 17:35:43] [RUN]  extract         [====--------------------]  18%  24/130 elapsed=8m04s
```

Keep all per-sample logs for debugging:

```bash
meta_marker_count \
  --input data_path.tsv \
  --outdir marker_count_out \
  --keep-task-logs
```

---

## R analysis templates

The R scripts in `scripts/` are analysis templates. They are intended for RStudio or interactive R use, not as strict command-line tools.

General workflow:

1. open one template in RStudio;
2. edit the `MC_CONFIG <- list(...)` block near the top;
3. run the script;
4. collect outputs under the configured `outdir`.

### Shared conventions

Most templates write:

```text
OUTDIR/
├── tables/
├── figures/
└── logs/
```

Figures use a shared clean theme and restrained palette from `scripts/analysis_common.R`, with PDF and 320-dpi PNG output by default.

### Template summary

| Category | Template | Main outputs | Default approach |
|---|---|---|---|
| Diversity | `scripts/diversity_alpha_beta.R` | alpha diversity, PCoA, PERMANOVA, beta dispersion | Bray-Curtis, 9999 permutations |
| Differential abundance | `scripts/differential_abundance_compositional.R` | contrasts, volcano, top feature bars | CLR transform, linear models, BH-FDR |
| Biomarkers | `scripts/biomarker_discovery_nonparametric.R` | biomarker candidates, effect sizes, heatmap | Kruskal-Wallis, Wilcoxon, BH-FDR |
| Correlation / network | `scripts/network_association_analysis.R` | edge table, node metrics, network plot | CLR + Spearman + BH-FDR |
| Functional profile | `scripts/functional_profile_analysis.R` | top functions, PCoA, differential tests | Bray-Curtis + CLR models |
| Taxa over time/space | `scripts/taxa_time_rank_spacetime.R` | target taxa dynamics, rank diversity | long-table summaries |
| Disturbance response | `scripts/disturbance_event_response.R` | event PERMANOVA, env tests, focal taxa models | longitudinal screening |

### Matrix input example

```tsv
sample_id	year	month	depth	Bacteria|16S|genus__Vibrio	Fungi|ITS|genus__Fusarium
S1	2023	11	00-10	12.5	3.2
S2	2023	11	10-20	0.0	8.1
```

### Long-table input example

```tsv
sample_id	marker	domain	rank	lineage	taxon_marker_reads	clean_reads_total	marker_rpm	marker_rpkm
S1	16S	Bacteria	genus	Bacteria;p__Pseudomonadota;g__Vibrio	20	1000000	20	0.13
```

### Example: diversity analysis

This example assumes a sample-by-feature matrix prepared from `all.marker_rpm.genus.long.tsv`.

```r
MC_CONFIG <- list(
  table = "marker_count_out/all.marker_rpm.genus.matrix.tsv",
  metadata = "",
  `table-format` = "sample_rows",
  `metadata-cols` = "year,month,depth,site_type",
  group = "site_type",
  covariates = "year,month",
  strata = "",
  outdir = "analysis_out/diversity_genus"
)
```

### Example: taxa time-space analysis

```r
MC_CONFIG <- list(
  `long-table` = "marker_count_out/all.marker_rpm.genus.long.tsv",
  metadata = "metadata.tsv",
  taxa = "Vibrio,Fusarium",
  `time-col` = "month",
  `space-col` = "site",
  `group-col` = "site_type",
  `target-rank` = "genus",
  outdir = "analysis_out/taxa_time_space"
)
```

---

## Biological interpretation notes

`meta_marker_count` reports marker-derived read signals, not absolute cell counts.

Important caveats:

- RPM normalizes by sequencing depth but not by marker copy number.
- `marker_rpkm` adds reference-length normalization, but it is still not a genome-copy correction.
- 16S, 18S, and ITS have different biological and technical properties; compare markers carefully.
- ITS length and database coverage are uneven, so genus-level summaries are usually safer than species-level interpretation.
- Paired reads are arbitrated before counting; R1/R2 conflicts are counted in `discordant_pairs`, and both mates follow the better mate classification.
- For time series or spatial designs, use metadata-aware statistical models rather than simple group means only.

Recommended reporting:

- marker set and reference versions;
- identity, alignment length, and qcov thresholds;
- read-depth normalization formula;
- number of samples and reads per group;
- rank used for ecological comparison;
- whether organellar references were removed;
- how paired-end discordant reads and deterministic tie-breaks were handled.

---

## Troubleshooting

### `seqkit` seems slow or CPU usage is low

`seqkit stats` may be I/O or gzip-decompression limited. Prefer sample-level parallelism first:

```bash
--jobs 4 --threads-per-sample 8
```

By default, seqkit is capped at 4 threads and extra CPU budget is used for more sample-level parallelism. Avoid forcing too many threads onto a single compressed FASTQ file.

### BBDuk reports Java heap errors

Increase BBDuk memory:

```bash
--bbmap-mem 120G
```

On shared clusters, request matching memory from the scheduler.

### The run appears stuck during extract or align

Use progress output and inspect failed logs:

```bash
ls marker_count_out/.tmp/task_logs
```

For debugging:

```bash
--keep-task-logs --progress-interval 2
```

### Reference files are missing

Check resolved paths:

```bash
meta_marker_count --print-defaults
cat marker_count_out/run_config.tsv
```

Then rebuild references or run with:

```bash
--ref-dir /path/to/refs
```

### R scripts cannot find output tables

Use the current flat output paths. Long-table based scripts can point directly to:

```r
`long-table` = "marker_count_out/all.marker_rpm.genus.long.tsv"
```

Matrix-based templates need a sample-by-feature matrix prepared from the long table. Do not use the old nested path `marker_count_out/04_abundance/...`.

---

## Checks and diagnostics

```bash
make check
meta_marker_count --check-deps
meta_marker_count --print-defaults
meta_marker_count --help
```

---

## License

This project is released under the MIT License. See [`LICENSE`](LICENSE) for details.
