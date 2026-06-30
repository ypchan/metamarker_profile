# metamarker_profile

<p align="center">
  <b>Fast marker-read quantification for metagenomic paired-end clean reads</b><br>
  16S / 18S / ITS candidate-read extraction, reference alignment, taxonomy-aware assignment, and marker RPM/RPKM tables.
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

## What is `metamarker_profile`?

`metamarker_profile` is a Python3-based bioinformatics workflow for estimating marker-derived reads from paired-end metagenomic clean reads. It is designed for ecological comparison of bacterial, archaeal, and fungal marker signals across many samples.

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
metamarker_profile/
├── pyproject.toml                    # Python package metadata
├── metamarker_profile/
│   ├── cli.py                        # main pipeline
│   └── __init__.py
├── scripts/
│   ├── build_ref_files.py            # optional reference preparation helper
│   ├── analysis_common.R             # shared R plotting/stat helpers
│   ├── extract_lineage_reads_assemble.py
│   ├── diversity_alpha_beta.R
│   ├── differential_abundance_compositional.R
│   ├── biomarker_discovery_nonparametric.R
│   ├── network_association_analysis.R
│   ├── functional_profile_analysis.R
│   ├── taxa_time_rank_spacetime.R
│   └── disturbance_event_response.R
├── refs/                             # optional local reference FASTA/index/taxonomy files
├── README.md
└── LICENSE
```

---

## Installation

### 1. Optional: create an environment

Use a dedicated Python environment so Python dependencies do not mix with system packages:

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -U pip
```

### 2. Install directly from GitHub

For normal use, install the latest repository version directly with pip:

```bash
python3 -m pip install "git+https://github.com/ypchan/metamarker_profile.git"
```

### 3. Or clone and install locally

```bash
gh repo clone ypchan/metamarker_profile
cd metamarker_profile
python3 -m pip install .
```

or:

```bash
git clone https://github.com/ypchan/metamarker_profile.git
cd metamarker_profile
python3 -m pip install .
```

For development, use editable mode from the cloned repository:

```bash
python3 -m pip install -e .
```

Installed commands:

```bash
which metamarker_profile
metamarker_profile --help
```

The Python package install covers Python dependencies. External tools such as `bbduk.sh`, `minimap2`, and optionally `seqkit` still need to be available in `PATH`.

### 4. Syntax check

```bash
python3 -m compileall metamarker_profile
```

This checks the installed Python package source. If `Rscript` is available, R analysis templates can be parsed separately:

```bash
Rscript -e 'files <- list.files("scripts", pattern="[.]R$", full.names=TRUE); invisible(lapply(files, parse))'
```

---

## Update

For direct GitHub installs, upgrade with:

```bash
python3 -m pip install -U "git+https://github.com/ypchan/metamarker_profile.git"
```

For a cloned regular pip install, update after pulling new code:

```bash
cd metamarker_profile
git pull --ff-only
python3 -m pip install -U .
```

For an editable install, pull the repository and rerun checks:

```bash
git pull --ff-only
python3 -m compileall metamarker_profile
```

---

## Dependencies

### Required runtime tools

| Tool | Used in | Purpose |
|---|---:|---|
| `python3` | all steps | workflow driver, scheduling, PAF parsing, table generation |
| `polars` | Step 04 | high-throughput PAF parsing and abundance aggregation |
| `rich`, `rich-argparse` | CLI | help text, progress display, structured logging |
| `bbduk.sh` | Step 02 | marker candidate-read extraction |
| `minimap2` | Step 03 | alignment to marker reference indexes |

`seqkit` is optional. With `--read-count-method auto`, the pipeline uses `seqkit` when available and falls back to native Python FASTQ counting otherwise.

### Optional tools

| Tool | Purpose |
|---|---|
| `miniasm` | assemble lineage-extracted `merged_segments.fq.gz` from minimap2 all-vs-all overlaps |
| `Rscript` | syntax checks for R analysis templates |
| `mamba` / `conda` | environment management |
| `aria2c` | faster reference download |

### Example conda environment

```bash
mamba create -n metamarker_profile -c conda-forge -c bioconda \
  python polars rich rich-argparse seqkit minimap2 miniasm bbmap \
  r-base r-ggplot2 r-vegan r-igraph

mamba activate metamarker_profile
```

---

## Reference database

`metamarker_profile` expects a compact reference directory containing processed FASTA files, minimap2 indexes, and a combined taxonomy table.

Default reference directory:

```text
./refs
```

Reference directory resolution order:

1. `METAMARKER_PROFILE_REF_DIR`
2. `METAMARKER_PROFILE_REF_CONFIG`, defaulting to `~/.config/metamarker_profile/ref_dir`
3. `./refs` in the current working directory

You can always override it at runtime:

```bash
metamarker_profile --ref-dir /path/to/refs ...
```

To set a persistent default reference directory:

```bash
mkdir -p ~/.config/metamarker_profile
printf '%s\n' /path/to/refs > ~/.config/metamarker_profile/ref_dir
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

### Download raw SILVA SSURef NR99

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

### Build minimap2 indexes

If the `.mmi` index files are not already present, create them from the processed FASTA files:

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
ref_id	marker	domain	euk_group	phylum	class	order	family	genus	species
AB000393.1.1510	16S	Bacteria	NA	Pseudomonadota	Gammaproteobacteria	Enterobacterales	Vibrionaceae	Vibrio	Vibrio_halioticoli
UDB016649	ITS	Fungi	Fungi	Ascomycota	Sordariomycetes	Hypocreales	Nectriaceae	Fusarium	Fusarium_sp
```

Reference parsing rules:

- SILVA `Bacteria` and `Archaea` records are treated as `16S`.
- SILVA `Eukaryota` records are treated as `18S`.
- UNITE records are treated as `ITS`.
- `euk_group` is metadata, not an eighth taxonomy rank. Use it to separate fungal eukaryotes from other eukaryotic groups. Recommended values include `Fungi`, `Metazoa`, `Viridiplantae`, `Alveolata`, `Stramenopiles`, `Rhizaria`, and `Amoebozoa`; use `NA` for non-eukaryotic records.
- Taxonomy files without `euk_group` are rejected. Rebuild refs with `scripts/build_ref_files.py` so every record has an explicit `euk_group`.
- SILVA records with exact lineage item `Mitochondria` or `Chloroplast` are removed by default.

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

---

## Run the pipeline

### Default 16S + ITS run

```bash
metamarker_profile \
  --input data_path.tsv \
  --outdir metamarker_profile_out \
  --markers 16S,ITS \
  --rank genus \
  --jobs 4 \
  --threads-per-sample 10
```

### 16S + 18S + ITS run

```bash
metamarker_profile \
  --input data_path.tsv \
  --outdir metamarker_profile_out \
  --markers 16S,18S,ITS \
  --rank genus \
  --jobs 4 \
  --threads-per-sample 10
```

### Single-sample run

```bash
metamarker_profile \
  --sample-id 202311_MF1_00-10 \
  --r1 /data/202311_MF1_00-10_R1.fq.gz \
  --r2 /data/202311_MF1_00-10_R2.fq.gz \
  --outdir metamarker_profile_out
```

### Dry run

Use this before long cluster runs:

```bash
metamarker_profile \
  --input data_path.tsv \
  --outdir metamarker_profile_out \
  --markers 16S,ITS \
  --dry-run
```

### Run only selected steps

```bash
# only count reads
metamarker_profile --input data_path.tsv --outdir metamarker_profile_out --steps reads_count

# extract + align only
metamarker_profile --input data_path.tsv --outdir metamarker_profile_out --steps extract,align

# recompute abundance from existing PAF files
metamarker_profile --input data_path.tsv --outdir metamarker_profile_out --steps abundance --rank genus
```

---

## Step-by-step parameters

### Step 01 — counting reads

Counts reads from clean FASTQ files with `seqkit` and writes `reads_stat.tsv`. R1/R2 may contain different numbers of reads; the workflow records both counts and continues.

| Option | Default | Meaning |
|---|---:|---|
| `--seqkit-threads` | `auto` | threads used by `seqkit stats` per sample |
| `--seqkit-max-threads` | `4` | cap for seqkit threads because gzip/I/O usually bottlenecks first |
| `--jobs` | `4` | maximum concurrent workflow tasks |

Output:

```text
reads_stat.tsv
```

Core columns:

```text
sample_id  year  month  depth  r1_path  r2_path  r1_reads  r2_reads  read_count_delta  read_pairs  clean_reads_total
```

Normalization denominator:

```text
clean_reads_total = R1 read count + R2 read count
```

If R1/R2 read counts differ, `read_count_delta` records the difference. Downstream pairing is based on read IDs, not record order.
`read_pairs` is retained as a compatibility estimate (`min(r1_reads, r2_reads)`); abundance normalization uses `clean_reads_total`.

### Step 02 — extract marker candidate reads

Uses `BBDuk` to extract marker-matching R1 and R2 reads independently. This allows input files with different read counts; downstream alignment restores mate relationships by read ID.

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

Uses one `minimap2` process per sample-marker task to align extracted R1/R2 reads together and writes a mate-labeled PAF file. Query IDs are normalized internally to `/1` and `/2` mate suffixes so the Polars parser can pair mates by read ID and arbitrate mates from one PAF.

PAF is used because this workflow quantifies marker-read support and taxonomic assignment, not insert-size or genome-coordinate paired-end evidence.

| Option | Default | Meaning |
|---|---:|---|
| `--minimap2-preset` | `sr` | minimap2 preset for short reads |
| `--minimap2-N` | `10` | retain up to N secondary alignments |
| `--minimap2-threads` | `auto` | minimap2 threads per sample |
| `--minimap2-max-threads` | `12` | cap for minimap2 threads |

Output:

```text
03_align/
├── 16S/<sample>.16S.pe.paf
└── ITS/<sample>.ITS.pe.paf
```

### Step 04 — assign taxonomy and calculate abundance

Parses paired PAF files, filters alignments, assigns each read pair to the single best PAF hit, and writes long/stat tables.

| Option | 16S default | 18S default | ITS default | Meaning |
|---|---:|---:|---:|---|
| `--min-identity-*` | `0.97` | `0.97` | `0.95` | minimum alignment identity |
| `--min-aln-len-*` | `80` | `80` | `80` | minimum aligned length |
| `--min-qcov-*` | `0.60` | `0.60` | `0.60` | minimum query coverage |
| `--min-mapq` | colspan | colspan | `0` | minimum MAPQ |
| `--rank` | colspan | colspan | `species` | `domain,phylum,class,order,family,genus,species,all` |

Best-hit assignment order:

1. identity
2. aligned length
3. matching bases
4. MAPQ
5. query coverage
6. target length
7. mate (`R1` before `R2` only as a deterministic final biological tie-break)
8. target ID

Mate handling is read-level counting with pair-level arbitration by normalized read ID:

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
all.marker_rpm.assignment_qc.tsv
```

`assignment_stats` preserves raw assignment counts per requested rank. `assignment_qc` is a rank-independent diagnostic table with rates derived from the same counts; it is not used to filter results.

---

## Extract Reads For Lineage Validation

Use `scripts/extract_lineage_reads_assemble.py` to pull out all marker reads assigned to one lineage and prepare sequence-level validation inputs.

Example:

```bash
python3 scripts/extract_lineage_reads_assemble.py \
  --outdir metamarker_profile_out \
  --rank genus \
  --lineage "Targetus" \
  --markers 16S \
  --jobs 8 \
  --assembler miniasm \
  --assembly-threads 8
```

The script reuses the same PAF filters and R1/R2 arbitration rules as the abundance step. In default `--mode assigned`, reads are selected from pairs assigned to the requested lineage after pair arbitration. This matches the abundance table: single-end reads are retained, and discordant R1/R2 pairs follow the better mate. For a stricter annotation check, use `--mode mate-hit` to keep only mates whose own best hit matches the lineage.

`--lineage` accepts either the short value at `--rank` or the complete lineage path from `all.marker_rpm.<rank>.long.tsv`. `selected_reads.tsv` reports both short rank labels and complete lineage paths for mate-level and assigned classifications.

Main outputs:

```text
lineage_extract/<rank>_<lineage>/
├── reads.R1.fq.gz
├── reads.R2.fq.gz
├── merged_segments.fq.gz
├── merged_segments.fa
├── miniasm.overlaps.paf
├── miniasm.gfa
├── miniasm_contigs.fa
├── assembly_summary.tsv
├── selected_reads.tsv
├── summary.tsv
├── run_status.tsv
└── extract_lineage_reads_assemble.log
```

`merged_segments.fq.gz` is the required assembly input. It contains the selected R1/R2 read pool, oriented by the best PAF strand before assembly. With `--assembler miniasm`, the script runs `minimap2` all-vs-all overlap search on `merged_segments.fq.gz`, sends the resulting PAF to `miniasm`, and exports GFA segments to `miniasm_contigs.fa`.

The script uses Rich help, colored run logs, progress bars, parallel PAF parsing with `--jobs`, and minimap2 threading with `--assembly-threads`. Resume is enabled by default: a stage is skipped only when the run log contains its completion marker and the expected output files exist. Use `--force` to recompute, `--no-resume` to disable checkpoint checks, `--help_input` for exact input format requirements, and `--help_default` for default parameters.

You can also assemble an existing merged segments file directly:

```bash
python3 scripts/extract_lineage_reads_assemble.py \
  --segments-fastq metamarker_profile_out/lineage_extract/genus_Targetus/merged_segments.fq.gz \
  --assembler miniasm \
  --overlap-preset ava-ont \
  --miniasm-opts "-m 40 -s 40"
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
| `--outdir DIR` | `metamarker_profile_out` | output directory |
| `--steps LIST` | `all` | `reads_count,extract,align,abundance` or `all` |
| `--markers LIST` | `16S,ITS` | marker set: `16S`, `18S`, `ITS` |
| `--rank RANK` | `species` | output taxonomy rank or `all` |

### Parallelism and HPC resources

| Option | Default | Description |
|---|---:|---|
| `--jobs INT` | `4` | maximum concurrent workflow tasks |
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
step parallelism = min(jobs, task count, floor(CPU budget / effective tool threads))
extract/align task count = samples × markers
abundance task count = samples
```

Example for a 40-core node:

```bash
metamarker_profile \
  --input data_path.tsv \
  --outdir metamarker_profile_out \
  --jobs 4 \
  --threads-per-sample 10 \
  --bbmap-mem 80G
```

With this example, seqkit uses up to 4 threads, while BBDuk/minimap2 use up to 10 threads. If `--threads-per-sample 24` is used, BBDuk/minimap2 are capped at 12; increase `--jobs` explicitly to run more workflow tasks at the same time.

### Resume, overwrite, and logs

| Option | Meaning |
|---|---|
| `--force` | rerun and overwrite existing outputs |
| `--no-resume` | ignore checkpoints |
| `--clean-tmp` | remove temporary files after success |
| `--dry-run` | validate inputs and print run plan only |
| `--check-deps` | check external tools and exit |
| `--keep-task-logs` | keep per-sample logs even when tasks succeed |

---

## Output layout

The current output layout keeps important result files directly under `OUTDIR` and stores large intermediate files in numbered directories.

```text
metamarker_profile_out/
├── sample_manifest.tsv
├── run_config.tsv
├── metamarker_profile.YYYYmmdd_HHMMSS.log
├── commands/
│   ├── reads_count/
│   ├── extract/
│   ├── align/
│   └── abundance/
├── reads_stat.tsv
├── all.marker_rpm.genus.long.tsv
├── all.marker_rpm.genus.assignment_stats.tsv
├── all.marker_rpm.assignment_qc.tsv
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
| `commands/` | executed command lines grouped by stage and sample |
| `reads_stat.tsv` | R1/R2 read counts, count deltas, legacy read-pair estimate, and clean-read totals |
| `all.marker_rpm.<rank>.long.tsv` | tidy long abundance table |
| `all.marker_rpm.<rank>.assignment_stats.tsv` | alignment filtering, paired-end discordance, reassigned discordant reads, tie-breaking, and assignment statistics |
| `all.marker_rpm.assignment_qc.tsv` | rank-independent diagnostic rates for alignment pass rate, assignment rate, and paired-end discordance |

### Long-table columns

```text
sample_id
marker
domain
euk_group
rank
lineage
rank_lineage
taxon_marker_reads
clean_reads_total
marker_rpm
marker_rpkm
mean_target_length_bp
<additional sample metadata columns>
```

`lineage` is the complete taxonomy path through the requested `rank`, formatted with rank prefixes.
The default `--rank species` writes seven ranks:
`d__<domain>;p__<phylum>;c__<class>;o__<order>;f__<family>;g__<genus>;s__<species>`.
For `--rank genus`, `lineage` writes six ranks and omits species:
`d__<domain>;p__<phylum>;c__<class>;o__<order>;f__<family>;g__<genus>`.
`euk_group` is reported separately and is not embedded in `lineage`.
`rank_lineage` is the short value at the requested `rank`. Taxonomy values with empty rank fields are reported as `unidentified`; missing domains are reported as `NA`.

---

## Logging and progress

The screen output is designed for HPC runs:

- dependency checks are printed once at the beginning;
- long-running steps use progress bars;
- executed command lines are preserved under `commands/`;
- successful task logs are removed by default;
- failed sample or sample-marker task logs are preserved under `.tmp/task_logs/`;
- the main log is kept as `metamarker_profile.YYYYmmdd_HHMMSS.log`.

Example:

```text
[2026-06-18 17:24:12] [RUN]  counting reads  [======------------------]  25%  16/65  elapsed=3m10s
[2026-06-18 17:35:43] [DONE] counting reads  [========================] 100%  65/65  elapsed=11m31s
[2026-06-18 17:35:43] [RUN]  extract         [====--------------------]  18%  24/130 elapsed=8m04s
```

Keep all per-sample logs for debugging:

```bash
metamarker_profile \
  --input data_path.tsv \
  --outdir metamarker_profile_out \
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

### Long-table input example

```tsv
sample_id	marker	domain	euk_group	rank	lineage	rank_lineage	taxon_marker_reads	clean_reads_total	marker_rpm	marker_rpkm
S1	16S	Bacteria	NA	genus	d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Vibrionales;f__Vibrionaceae;g__Vibrio	Vibrio	20	1000000	20	0.13
```

### Example: taxa time-space analysis

```r
MC_CONFIG <- list(
  `long-table` = "metamarker_profile_out/all.marker_rpm.genus.long.tsv",
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

`metamarker_profile` reports marker-derived read signals, not absolute cell counts.

Important caveats:

- RPM normalizes by sequencing depth but not by marker copy number.
- `marker_rpkm` is an RPKM-like reference-length normalization, not a genome-copy correction.
- 16S, 18S, and ITS have different biological and technical properties; interpret abundances within marker type rather than as direct cross-marker abundance comparisons.
- ITS length and database coverage are uneven, so genus-level summaries are usually safer than species-level interpretation.
- Paired reads are arbitrated before counting; R1/R2 conflicts are counted in `discordant_pairs`, and both mates follow the better mate classification.
- `all.marker_rpm.assignment_qc.tsv` reports diagnostic rates only. It does not apply or recommend ecological filtering.
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

### BBDuk reports different numbers of reads in paired input files

Older versions used BBDuk paired input mode and failed when R1/R2 read counts differed. Current versions extract R1 and R2 independently, so this error usually means the installed package is older than the source tree. Check both files with:

```bash
seqkit stats -T fq/ERR10446187_1.fastq.gz fq/ERR10446187_2.fastq.gz
```

If the read counts differ, rerun with the current `metamarker_profile` version and use `--force` or a fresh output directory so the extract step is regenerated. Pair-level arbitration later uses read IDs, so equal file lengths are not required.

### The run appears stuck during extract or align

Use progress output and inspect failed logs. Failed external commands also print the tail of their task log in the main error message:

```bash
ls metamarker_profile_out/.tmp/task_logs
```

For debugging, rerun with `--keep-task-logs`.

### Reference files are missing

Check resolved paths:

```bash
metamarker_profile --input data_path.tsv --dry-run
cat metamarker_profile_out/run_config.tsv
```

Then rebuild references or run with:

```bash
--ref-dir /path/to/refs
```

### R scripts cannot find output tables

Use the current flat output paths. The pipeline writes long tables, not sample-by-feature matrices. Long-table based scripts can point directly to:

```r
`long-table` = "metamarker_profile_out/all.marker_rpm.genus.long.tsv"
```

Some external analyses may require their own input conversion. This pipeline does not write matrix files.

---

## Checks and diagnostics

```bash
python3 -m compileall metamarker_profile
metamarker_profile --check-deps
metamarker_profile --help
```

---

## License

This project is released under the MIT License. See [`LICENSE`](LICENSE) for details.
