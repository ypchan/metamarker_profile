#!/usr/bin/env bash
set -Eeuo pipefail

PROGRAM="$(basename "$0")"

# ==============================================================================
# Defaults
# ==============================================================================
INPUT=""
SAMPLE_ID=""
R1_PATH=""
R2_PATH=""
# Single-sample metadata are intentionally not exposed as CLI options.
# They are recorded as NA when --sample-id/--r1/--r2 mode is used.
YEAR="NA"
MONTH="NA"
DEPTH="NA"

OUTDIR="marker_count_out"
STEPS="all"
MARKERS="16S,ITS"

REF_16S=""
REF_18S=""
REF_ITS=""
INDEX_16S=""
INDEX_18S=""
INDEX_ITS=""
TAXONOMY=""

RANK="genus"
JOBS=4
THREADS_PER_SAMPLE=4
SEQKIT_THREADS=1
BBMAP_MEM="64G"

# BBDuk candidate-read extraction defaults.
K_16S=31
HDIST_16S=0
MKH_16S=1
MINK_16S=0
K_18S=31
HDIST_18S=0
MKH_18S=1
MINK_18S=0
K_ITS=25
HDIST_ITS=0
MKH_ITS=1
MINK_ITS=0

# minimap2 short-read alignment defaults.
MINIMAP2_PRESET="sr"
MINIMAP2_N=10

# Abundance filtering defaults for genus-level ecological comparison.
MIN_IDENTITY_16S="0.97"
MIN_IDENTITY_18S="0.97"
MIN_IDENTITY_ITS="0.95"
MIN_ALN_LEN_16S=80
MIN_ALN_LEN_18S=80
MIN_ALN_LEN_ITS=80
MIN_QCOV_16S="0.60"
MIN_QCOV_18S="0.60"
MIN_QCOV_ITS="0.60"
MIN_MAPQ=0
TOP_IDENTITY_DIFF="0.010"
TOP_ALN_LEN_DIFF=10

FORCE=0
RESUME=1
CLEAN_TMP=0
DRY_RUN=0
CHECK_DEPS_ONLY=0

# Runtime paths; initialized in init_paths().
MANIFEST=""
LOG_DIR=""
STATUS_DIR=""
CLEAN_DIR=""
CLEAN_OUT=""
CLEAN_TMP_DIR=""
MARKER_DIR=""
ALIGN_DIR=""
ABUND_DIR=""
SCRIPT_DIR=""
MAIN_LOG=""

# ==============================================================================
# Logging and colors
# ==============================================================================
BOLD=""; DIM=""; RED=""; GREEN=""; YELLOW=""; BLUE=""; MAGENTA=""; CYAN=""; RESET=""

init_colors() {
    local use_color=0
    # NO_COLOR disables ANSI color. CLICOLOR_FORCE=1 or FORCE_COLOR>0 forces it.
    if [[ -n "${NO_COLOR:-}" ]]; then
        use_color=0
    elif [[ "${CLICOLOR_FORCE:-0}" == "1" || "${FORCE_COLOR:-0}" =~ ^[1-9] ]]; then
        use_color=1
    elif [[ -t 2 && "${TERM:-dumb}" != "dumb" ]]; then
        use_color=1
    fi

    if [[ "${use_color}" -eq 1 ]]; then
        BOLD=$'\033[1m'
        DIM=$'\033[2m'
        RED=$'\033[31m'
        GREEN=$'\033[32m'
        YELLOW=$'\033[33m'
        BLUE=$'\033[34m'
        MAGENTA=$'\033[35m'
        CYAN=$'\033[36m'
        RESET=$'\033[0m'
    fi
}

ts() { date '+%Y-%m-%d %H:%M:%S'; }

_log_plain() {
    local level="$1"; shift
    [[ -n "${MAIN_LOG:-}" ]] && printf '[%s] [%s] %s\n' "$(ts)" "${level}" "$*" >> "${MAIN_LOG}" || true
}

log() {
    printf '%b[%s] [%bINFO%b] %s%b\n' "${DIM}" "$(ts)" "${GREEN}" "${RESET}${DIM}" "$*" "${RESET}" >&2
    _log_plain INFO "$*"
}

log_step() {
    printf '%b[%s] [%bSTEP%b] %s%b\n' "${DIM}" "$(ts)" "${CYAN}" "${RESET}${DIM}" "$*" "${RESET}" >&2
    _log_plain STEP "$*"
}

log_warn() {
    printf '%b[%s] [%bWARN%b] %s%b\n' "${DIM}" "$(ts)" "${YELLOW}" "${RESET}${DIM}" "$*" "${RESET}" >&2
    _log_plain WARN "$*"
}

log_error() {
    printf '%b[%s] [%bERROR%b] %s%b\n' "${DIM}" "$(ts)" "${RED}" "${RESET}${DIM}" "$*" "${RESET}" >&2
    _log_plain ERROR "$*"
}

die() { log_error "$*"; exit 1; }

# ==============================================================================
# Help
# ==============================================================================
show_input_format() {
    cat <<'EOF_FORMAT'
Input format 1: multiple samples, TSV with header
--------------------------------------------------
Header is required. Columns are matched by column name, not by column number.
Column order can be arbitrary.

Required columns:
  sample_id
  r1_path
  r2_path

Optional metadata columns:
  year
  month
  depth
  any other columns, for example site_id, site_type, group, habitat, replicate

Notes:
  - year/month/depth are optional. Missing values are filled with NA.
  - Extra metadata columns are retained in the internal manifest.
  - Extra metadata columns are also carried into abundance output tables.
  - r1_path and r2_path must point to existing FASTQ(.gz) files.
  - Do not add inline comments after file paths; they will be treated as part of the path.
  - The examples below use real TAB characters and can be copied directly into a .tsv file.

Minimal valid TSV:
sample_id	r1_path	r2_path
S1	/path/S1_R1.fq.gz	/path/S1_R2.fq.gz
S2	/path/S2_R1.fq.gz	/path/S2_R2.fq.gz

Recommended TSV with metadata:
sample_id	year	month	depth	site_type	r1_path	r2_path
202311_MF1_00-10	2023	11	00-10	MF	/path/202311_MF1_00-10_R1.fq.gz	/path/202311_MF1_00-10_R2.fq.gz
202311_MF1_10-20	2023	11	10-20	MF	/path/202311_MF1_10-20_R1.fq.gz	/path/202311_MF1_10-20_R2.fq.gz

Also valid because columns are read by name:
r2_path	sample_id	depth	r1_path	month	year	site_type
/path/202311_MF1_00-10_R2.fq.gz	202311_MF1_00-10	00-10	/path/202311_MF1_00-10_R1.fq.gz	11	2023	MF

Input format 2: single sample from command line
-----------------------------------------------
marker_count.sh \
  --sample-id 202311_MF1_00-10 \
  --r1 /path/202311_MF1_00-10_R1.fq.gz \
  --r2 /path/202311_MF1_00-10_R2.fq.gz \
  --ref-16s ref/SILVA_16S.arc_bac.fasta \
  --ref-its ref/UNITE_ITS.fungi.fasta \
  --index-16s index/SILVA_16S.arc_bac.mmi \
  --index-its index/UNITE_ITS.fungi.mmi \
  --taxonomy ref_taxonomy.tsv

Taxonomy format for abundance step
----------------------------------
ref_id	marker	domain	phylum	class	order	family	genus	species
SILVA_16S_REF001	16S	Bacteria	p__Proteobacteria	c__...	o__...	f__...	g__...	s__...
SILVA_18S_REF001	18S	Eukaryota	p__Ascomycota	c__...	o__...	f__...	g__...	s__...
UNITE_ITS_REF001	ITS	Fungi	p__Ascomycota	c__...	o__...	f__...	g__...	s__...
EOF_FORMAT
}
show_defaults() {
    cat <<EOF_DEFAULTS
${PROGRAM} default parameters

Default markers
---------------
--markers ${MARKERS}
  Default runs 16S + ITS.
  18S is optional and can be added with --markers 16S,18S,ITS.

BBDuk candidate-read extraction
-------------------------------
16S: k=${K_16S}, hdist=${HDIST_16S}, mkh=${MKH_16S}, mink=${MINK_16S}
18S: k=${K_18S}, hdist=${HDIST_18S}, mkh=${MKH_18S}, mink=${MINK_18S}
ITS: k=${K_ITS}, hdist=${HDIST_ITS}, mkh=${MKH_ITS}, mink=${MINK_ITS}

Reasoning:
  16S and 18S are SSU markers, so k=31 gives stable specificity.
  ITS is more variable, so k=25 improves candidate-read recovery.
  hdist=0 avoids excessive memory use and reduces random k-mer matches.
  mkh=1 keeps the extraction step sensitive; later PAF filters control false positives.

minimap2 alignment
------------------
--minimap2-preset ${MINIMAP2_PRESET}
--minimap2-N ${MINIMAP2_N}

Reasoning:
  -x sr is appropriate for PE short reads.
  -N 10 preserves multiple good hits, allowing the abundance step to mark ambiguous reads.

Abundance filters
-----------------
16S: identity>=${MIN_IDENTITY_16S}, aln_len>=${MIN_ALN_LEN_16S}, qcov>=${MIN_QCOV_16S}
18S: identity>=${MIN_IDENTITY_18S}, aln_len>=${MIN_ALN_LEN_18S}, qcov>=${MIN_QCOV_18S}
ITS: identity>=${MIN_IDENTITY_ITS}, aln_len>=${MIN_ALN_LEN_ITS}, qcov>=${MIN_QCOV_ITS}
mapq>=${MIN_MAPQ}
near-top identity difference <= ${TOP_IDENTITY_DIFF}
near-top alignment-length difference <= ${TOP_ALN_LEN_DIFF}

Reasoning:
  Genus is the default rank because it is more stable than species for PE marker reads.
  16S/18S use 0.97 identity for conservative genus-level ecological comparison.
  ITS uses 0.95 because fungal ITS is more variable and database coverage is uneven.
  Near-top hits are compared before assignment; reads with conflicting near-best taxa are marked ambiguous.

Reads-count normalization
-------------------------
marker_rpm = taxon_marker_reads / clean_reads_total * 1,000,000
clean_reads_total = R1 read count * 2

EOF_DEFAULTS
}

show_help() {
    cat <<EOF_HELP
${BOLD}${PROGRAM}${RESET}

Integrated marker-read counting pipeline for paired-end clean reads.

Workflow:
  1. Count reads from clean FASTQ files with seqkit.
  2. Extract marker candidate reads with BBDuk.
  3. Align marker reads with minimap2 and write PAF.
  4. Parse PAF with shell/awk and calculate marker RPM.

${BOLD}Input TSV columns:${RESET}
  Multi-sample TSV must have a header. Columns are matched by name, not by position.
  Required column names: sample_id, r1_path, r2_path
  Optional metadata columns: year, month, depth, and any other sample metadata columns.
  Missing year/month/depth are filled with NA. Extra metadata columns are retained and carried into abundance outputs.

${BOLD}Required arguments - sample input:${RESET}
  ${CYAN}-i, --input${RESET} FILE                Multi-sample TSV input.
  ${DIM}or single-sample mode:${RESET}
  ${CYAN}--sample-id${RESET} ID                  Single-sample ID.
  ${CYAN}--r1${RESET} FILE                       Single-sample R1 FASTQ(.gz).
  ${CYAN}--r2${RESET} FILE                       Single-sample R2 FASTQ(.gz).

${BOLD}Required arguments - marker references:${RESET}
  ${CYAN}--ref-16s${RESET} FILE                  16S reference FASTA for BBDuk, required when 16S is selected.
  ${CYAN}--ref-18s${RESET} FILE                  18S reference FASTA for BBDuk, required when 18S is selected.
  ${CYAN}--ref-its${RESET} FILE                  ITS reference FASTA for BBDuk, required when ITS is selected.
  ${CYAN}--index-16s${RESET} FILE                minimap2 index for 16S, required when 16S is selected.
  ${CYAN}--index-18s${RESET} FILE                minimap2 index for 18S, required when 18S is selected.
  ${CYAN}--index-its${RESET} FILE                minimap2 index for ITS, required when ITS is selected.
  ${CYAN}--taxonomy${RESET} FILE                 ref_id/marker/domain/ranks taxonomy TSV for abundance.

${BOLD}Optional arguments - workflow control:${RESET}
  ${CYAN}-o, --outdir${RESET} DIR                Output directory. default: ${OUTDIR}
  ${CYAN}--steps${RESET} LIST                    all or comma list: reads_count,extract,align,abundance. default: ${STEPS}
                                      Compatibility alias: clean -> reads_count.
  ${CYAN}--markers${RESET} LIST                  Marker list: 16S,18S,ITS in any combination. default: ${MARKERS}
  ${CYAN}--rank${RESET} RANK                     domain,phylum,class,order,family,genus,species,all. default: ${RANK}
${BOLD}Optional arguments - parallelism and resources:${RESET}
  ${CYAN}-j, --jobs${RESET} INT                  Number of samples processed in parallel. default: ${JOBS}
  ${CYAN}-p, --threads-per-sample${RESET} INT    Threads used by BBDuk/minimap2 per sample. default: ${THREADS_PER_SAMPLE}
  ${CYAN}--seqkit-threads${RESET} INT            Threads used by seqkit per sample. default: ${SEQKIT_THREADS}
  ${CYAN}--bbmap-mem${RESET} STR                 BBDuk Java heap, e.g. 50G, 250G. default: ${BBMAP_MEM}

${BOLD}Optional arguments - BBDuk extraction parameters:${RESET}
  ${CYAN}--k-16s${RESET} INT                     default: ${K_16S}
  ${CYAN}--hdist-16s${RESET} INT                 default: ${HDIST_16S}
  ${CYAN}--mkh-16s${RESET} INT                   default: ${MKH_16S}
  ${CYAN}--mink-16s${RESET} INT                  0 disables. default: ${MINK_16S}
  ${CYAN}--k-18s${RESET} INT                     default: ${K_18S}
  ${CYAN}--hdist-18s${RESET} INT                 default: ${HDIST_18S}
  ${CYAN}--mkh-18s${RESET} INT                   default: ${MKH_18S}
  ${CYAN}--mink-18s${RESET} INT                  0 disables. default: ${MINK_18S}
  ${CYAN}--k-its${RESET} INT                     default: ${K_ITS}
  ${CYAN}--hdist-its${RESET} INT                 default: ${HDIST_ITS}
  ${CYAN}--mkh-its${RESET} INT                   default: ${MKH_ITS}
  ${CYAN}--mink-its${RESET} INT                  0 disables. default: ${MINK_ITS}

${BOLD}Optional arguments - minimap2 parameters:${RESET}
  ${CYAN}--minimap2-preset${RESET} STR           default: ${MINIMAP2_PRESET}
  ${CYAN}--minimap2-N${RESET} INT                default: ${MINIMAP2_N}

${BOLD}Optional arguments - abundance filters:${RESET}
  ${CYAN}--min-identity-16s${RESET} FLOAT        default: ${MIN_IDENTITY_16S}
  ${CYAN}--min-identity-18s${RESET} FLOAT        default: ${MIN_IDENTITY_18S}
  ${CYAN}--min-identity-its${RESET} FLOAT        default: ${MIN_IDENTITY_ITS}
  ${CYAN}--min-aln-len-16s${RESET} INT           default: ${MIN_ALN_LEN_16S}
  ${CYAN}--min-aln-len-18s${RESET} INT           default: ${MIN_ALN_LEN_18S}
  ${CYAN}--min-aln-len-its${RESET} INT           default: ${MIN_ALN_LEN_ITS}
  ${CYAN}--min-qcov-16s${RESET} FLOAT            default: ${MIN_QCOV_16S}
  ${CYAN}--min-qcov-18s${RESET} FLOAT            default: ${MIN_QCOV_18S}
  ${CYAN}--min-qcov-its${RESET} FLOAT            default: ${MIN_QCOV_ITS}
  ${CYAN}--min-mapq${RESET} INT                  default: ${MIN_MAPQ}
  ${CYAN}--top-identity-diff${RESET} FLOAT       default: ${TOP_IDENTITY_DIFF}
  ${CYAN}--top-aln-len-diff${RESET} INT          default: ${TOP_ALN_LEN_DIFF}

${BOLD}Optional arguments - resume, overwrite and logging:${RESET}
  ${CYAN}--force${RESET}                         Re-run and overwrite outputs.
  ${CYAN}--no-resume${RESET}                     Ignore done flags and output checks.
  ${CYAN}--clean-tmp${RESET}                     Remove temporary files after successful run.
  ${CYAN}--dry-run${RESET}                       Check inputs and print run plan only.
  ${CYAN}--check-deps${RESET}                    Check required external tools and exit.
  ${CYAN}--print-format${RESET}                  Print input format examples and exit.
  ${CYAN}--print-defaults${RESET}                Explain default parameters and exit.
  ${CYAN}-h, --help${RESET}                      Show this help.

${BOLD}Example:${RESET}
  ${CYAN}${PROGRAM}${RESET} \\
    --input data_path.tsv \\
    --outdir marker_count_out \\
    --markers 16S,ITS \\
    --ref-16s ref/SILVA_16S.arc_bac.fasta \\
    --ref-its ref/UNITE_ITS.fungi.fasta \\
    --index-16s index/SILVA_16S.arc_bac.mmi \\
    --index-its index/UNITE_ITS.fungi.mmi \\
    --taxonomy ref_taxonomy.tsv \\
    --jobs 8 \\
    --threads-per-sample 4

${BOLD}Structured outputs:${RESET}
  OUTDIR/00_manifest/samples.tsv
  OUTDIR/00_logs/marker_count.YYYYmmdd_HHMMSS.log
  OUTDIR/01_reads_count/clean_reads.tsv
  OUTDIR/02_marker_reads/{16S,18S,ITS,stats,logs,tmp}/
  OUTDIR/03_align/{16S,18S,ITS,logs,tmp}/
  OUTDIR/04_abundance/marker_rpm.<rank>.long.tsv
  OUTDIR/04_abundance/marker_rpm.<rank>.matrix.tsv
  OUTDIR/04_abundance/marker_rpm.<rank>.domain_total.tsv
  OUTDIR/04_abundance/marker_rpm.<rank>.assignment_stats.tsv
  OUTDIR/99_checkpoints/{reads_count,extract,align,abundance}/
    Internal checkpoint files for --resume. These are not result tables.

EOF_HELP
}

require_value() {
    local opt="$1" value="${2:-}"
    [[ -n "${value}" && "${value}" != -* ]] || die "Option ${opt} requires a value."
}

is_positive_int() { [[ "$1" =~ ^[1-9][0-9]*$ ]]; }
is_nonnegative_int() { [[ "$1" =~ ^[0-9]+$ ]]; }

normalize_marker_token() {
    echo "$1" | tr '[:lower:]' '[:upper:]'
}

normalize_markers() {
    local raw token out="" seen16=0 seen18=0 seenits=0
    raw="$(echo "$1" | tr -d ' ' | sed 's/;/,/g')"
    IFS=',' read -r -a toks <<< "${raw}"
    for token in "${toks[@]}"; do
        token="$(normalize_marker_token "${token}")"
        [[ -n "${token}" ]] || continue
        case "${token}" in
            16S) [[ "${seen16}" -eq 0 ]] && { out+="${out:+,}16S"; seen16=1; } ;;
            18S) [[ "${seen18}" -eq 0 ]] && { out+="${out:+,}18S"; seen18=1; } ;;
            ITS) [[ "${seenits}" -eq 0 ]] && { out+="${out:+,}ITS"; seenits=1; } ;;
            *) die "Unsupported marker: ${token}. Allowed: 16S,18S,ITS" ;;
        esac
    done
    [[ -n "${out}" ]] || die "--markers is empty. Allowed: 16S,18S,ITS"
    echo "${out}"
}

list_has_marker() {
    local marker="$1"
    [[ ",${MARKERS}," == *",${marker},"* ]]
}

for_each_marker() {
    local m
    IFS=',' read -r -a _marker_array <<< "${MARKERS}"
    for m in "${_marker_array[@]}"; do
        echo "${m}"
    done
}

normalize_steps() {
    local raw token step out="" seen_reads_count=0 seen_extract=0 seen_align=0 seen_abundance=0
    raw="$(echo "$1" | tr -d ' ' | tr '[:upper:]' '[:lower:]' | sed 's/;/,/g')"
    [[ "${raw}" == "all" ]] && { echo "all"; return 0; }

    IFS=',' read -r -a toks <<< "${raw}"
    for token in "${toks[@]}"; do
        [[ -n "${token}" ]] || continue
        case "${token}" in
            reads_count|read_count|reads-count|read-count|count|counts) step="reads_count" ;;
            clean) step="reads_count"; log_warn "Deprecated step name 'clean' detected; use 'reads_count' instead." ;;
            extract|align|abundance) step="${token}" ;;
            *) die "Unsupported step: ${token}. Allowed: all or reads_count,extract,align,abundance" ;;
        esac

        case "${step}" in
            reads_count) [[ "${seen_reads_count}" -eq 0 ]] && { out+="${out:+,}${step}"; seen_reads_count=1; } ;;
            extract) [[ "${seen_extract}" -eq 0 ]] && { out+="${out:+,}${step}"; seen_extract=1; } ;;
            align) [[ "${seen_align}" -eq 0 ]] && { out+="${out:+,}${step}"; seen_align=1; } ;;
            abundance) [[ "${seen_abundance}" -eq 0 ]] && { out+="${out:+,}${step}"; seen_abundance=1; } ;;
        esac
    done
    [[ -n "${out}" ]] || die "--steps is empty. Allowed: all or reads_count,extract,align,abundance"
    echo "${out}"
}

should_run_step() {
    local step="$1"
    [[ "${STEPS}" == "all" || ",${STEPS}," == *",${step},"* ]]
}

parse_args() {
    if [[ "$#" -eq 0 ]]; then
        show_help
        exit 0
    fi

    while [[ "$#" -gt 0 ]]; do
        case "$1" in
            -i|--input) require_value "$1" "${2:-}"; INPUT="$2"; shift 2 ;;
            --input=*) INPUT="${1#*=}"; shift ;;
            --sample-id) require_value "$1" "${2:-}"; SAMPLE_ID="$2"; shift 2 ;;
            --sample-id=*) SAMPLE_ID="${1#*=}"; shift ;;
            --r1) require_value "$1" "${2:-}"; R1_PATH="$2"; shift 2 ;;
            --r1=*) R1_PATH="${1#*=}"; shift ;;
            --r2) require_value "$1" "${2:-}"; R2_PATH="$2"; shift 2 ;;
            --r2=*) R2_PATH="${1#*=}"; shift ;;

            -o|--outdir) require_value "$1" "${2:-}"; OUTDIR="$2"; shift 2 ;;
            --outdir=*) OUTDIR="${1#*=}"; shift ;;
            --steps) require_value "$1" "${2:-}"; STEPS="$2"; shift 2 ;;
            --steps=*) STEPS="${1#*=}"; shift ;;
            --markers) require_value "$1" "${2:-}"; MARKERS="$2"; shift 2 ;;
            --markers=*) MARKERS="${1#*=}"; shift ;;

            --ref-16s) require_value "$1" "${2:-}"; REF_16S="$2"; shift 2 ;;
            --ref-16s=*) REF_16S="${1#*=}"; shift ;;
            --ref-18s) require_value "$1" "${2:-}"; REF_18S="$2"; shift 2 ;;
            --ref-18s=*) REF_18S="${1#*=}"; shift ;;
            --ref-its) require_value "$1" "${2:-}"; REF_ITS="$2"; shift 2 ;;
            --ref-its=*) REF_ITS="${1#*=}"; shift ;;
            --index-16s) require_value "$1" "${2:-}"; INDEX_16S="$2"; shift 2 ;;
            --index-16s=*) INDEX_16S="${1#*=}"; shift ;;
            --index-18s) require_value "$1" "${2:-}"; INDEX_18S="$2"; shift 2 ;;
            --index-18s=*) INDEX_18S="${1#*=}"; shift ;;
            --index-its) require_value "$1" "${2:-}"; INDEX_ITS="$2"; shift 2 ;;
            --index-its=*) INDEX_ITS="${1#*=}"; shift ;;
            --taxonomy) require_value "$1" "${2:-}"; TAXONOMY="$2"; shift 2 ;;
            --taxonomy=*) TAXONOMY="${1#*=}"; shift ;;

            -j|--jobs) require_value "$1" "${2:-}"; JOBS="$2"; shift 2 ;;
            --jobs=*) JOBS="${1#*=}"; shift ;;
            -p|--threads-per-sample) require_value "$1" "${2:-}"; THREADS_PER_SAMPLE="$2"; shift 2 ;;
            --threads-per-sample=*) THREADS_PER_SAMPLE="${1#*=}"; shift ;;
            --seqkit-threads) require_value "$1" "${2:-}"; SEQKIT_THREADS="$2"; shift 2 ;;
            --seqkit-threads=*) SEQKIT_THREADS="${1#*=}"; shift ;;
            --bbmap-mem) require_value "$1" "${2:-}"; BBMAP_MEM="$2"; shift 2 ;;
            --bbmap-mem=*) BBMAP_MEM="${1#*=}"; shift ;;

            --k-16s) require_value "$1" "${2:-}"; K_16S="$2"; shift 2 ;;
            --k-16s=*) K_16S="${1#*=}"; shift ;;
            --hdist-16s) require_value "$1" "${2:-}"; HDIST_16S="$2"; shift 2 ;;
            --hdist-16s=*) HDIST_16S="${1#*=}"; shift ;;
            --mkh-16s) require_value "$1" "${2:-}"; MKH_16S="$2"; shift 2 ;;
            --mkh-16s=*) MKH_16S="${1#*=}"; shift ;;
            --mink-16s) require_value "$1" "${2:-}"; MINK_16S="$2"; shift 2 ;;
            --mink-16s=*) MINK_16S="${1#*=}"; shift ;;
            --k-18s) require_value "$1" "${2:-}"; K_18S="$2"; shift 2 ;;
            --k-18s=*) K_18S="${1#*=}"; shift ;;
            --hdist-18s) require_value "$1" "${2:-}"; HDIST_18S="$2"; shift 2 ;;
            --hdist-18s=*) HDIST_18S="${1#*=}"; shift ;;
            --mkh-18s) require_value "$1" "${2:-}"; MKH_18S="$2"; shift 2 ;;
            --mkh-18s=*) MKH_18S="${1#*=}"; shift ;;
            --mink-18s) require_value "$1" "${2:-}"; MINK_18S="$2"; shift 2 ;;
            --mink-18s=*) MINK_18S="${1#*=}"; shift ;;
            --k-its) require_value "$1" "${2:-}"; K_ITS="$2"; shift 2 ;;
            --k-its=*) K_ITS="${1#*=}"; shift ;;
            --hdist-its) require_value "$1" "${2:-}"; HDIST_ITS="$2"; shift 2 ;;
            --hdist-its=*) HDIST_ITS="${1#*=}"; shift ;;
            --mkh-its) require_value "$1" "${2:-}"; MKH_ITS="$2"; shift 2 ;;
            --mkh-its=*) MKH_ITS="${1#*=}"; shift ;;
            --mink-its) require_value "$1" "${2:-}"; MINK_ITS="$2"; shift 2 ;;
            --mink-its=*) MINK_ITS="${1#*=}"; shift ;;

            --minimap2-preset) require_value "$1" "${2:-}"; MINIMAP2_PRESET="$2"; shift 2 ;;
            --minimap2-preset=*) MINIMAP2_PRESET="${1#*=}"; shift ;;
            --minimap2-N) require_value "$1" "${2:-}"; MINIMAP2_N="$2"; shift 2 ;;
            --minimap2-N=*) MINIMAP2_N="${1#*=}"; shift ;;

            --rank) require_value "$1" "${2:-}"; RANK="$2"; shift 2 ;;
            --rank=*) RANK="${1#*=}"; shift ;;
            --min-identity-16s) require_value "$1" "${2:-}"; MIN_IDENTITY_16S="$2"; shift 2 ;;
            --min-identity-16s=*) MIN_IDENTITY_16S="${1#*=}"; shift ;;
            --min-identity-18s) require_value "$1" "${2:-}"; MIN_IDENTITY_18S="$2"; shift 2 ;;
            --min-identity-18s=*) MIN_IDENTITY_18S="${1#*=}"; shift ;;
            --min-identity-its) require_value "$1" "${2:-}"; MIN_IDENTITY_ITS="$2"; shift 2 ;;
            --min-identity-its=*) MIN_IDENTITY_ITS="${1#*=}"; shift ;;
            --min-aln-len-16s) require_value "$1" "${2:-}"; MIN_ALN_LEN_16S="$2"; shift 2 ;;
            --min-aln-len-16s=*) MIN_ALN_LEN_16S="${1#*=}"; shift ;;
            --min-aln-len-18s) require_value "$1" "${2:-}"; MIN_ALN_LEN_18S="$2"; shift 2 ;;
            --min-aln-len-18s=*) MIN_ALN_LEN_18S="${1#*=}"; shift ;;
            --min-aln-len-its) require_value "$1" "${2:-}"; MIN_ALN_LEN_ITS="$2"; shift 2 ;;
            --min-aln-len-its=*) MIN_ALN_LEN_ITS="${1#*=}"; shift ;;
            --min-qcov-16s) require_value "$1" "${2:-}"; MIN_QCOV_16S="$2"; shift 2 ;;
            --min-qcov-16s=*) MIN_QCOV_16S="${1#*=}"; shift ;;
            --min-qcov-18s) require_value "$1" "${2:-}"; MIN_QCOV_18S="$2"; shift 2 ;;
            --min-qcov-18s=*) MIN_QCOV_18S="${1#*=}"; shift ;;
            --min-qcov-its) require_value "$1" "${2:-}"; MIN_QCOV_ITS="$2"; shift 2 ;;
            --min-qcov-its=*) MIN_QCOV_ITS="${1#*=}"; shift ;;
            --min-mapq) require_value "$1" "${2:-}"; MIN_MAPQ="$2"; shift 2 ;;
            --min-mapq=*) MIN_MAPQ="${1#*=}"; shift ;;
            --top-identity-diff) require_value "$1" "${2:-}"; TOP_IDENTITY_DIFF="$2"; shift 2 ;;
            --top-identity-diff=*) TOP_IDENTITY_DIFF="${1#*=}"; shift ;;
            --top-aln-len-diff) require_value "$1" "${2:-}"; TOP_ALN_LEN_DIFF="$2"; shift 2 ;;
            --top-aln-len-diff=*) TOP_ALN_LEN_DIFF="${1#*=}"; shift ;;

            --force) FORCE=1; shift ;;
            --no-resume) RESUME=0; shift ;;
            --clean-tmp) CLEAN_TMP=1; shift ;;
            --dry-run) DRY_RUN=1; shift ;;
            --check-deps) CHECK_DEPS_ONLY=1; shift ;;
            --print-format) show_input_format; exit 0 ;;
            --print-defaults) show_defaults; exit 0 ;;
            -h|--help) show_help; exit 0 ;;
            --) shift; break ;;
            -*) die "Unknown option: $1. Run '${PROGRAM} --help'." ;;
            *) die "Unexpected positional argument: $1. Run '${PROGRAM} --help'." ;;
        esac
    done
}

# ==============================================================================
# Dependency checks
# ==============================================================================
need_tool() {
    local tool="$1" path
    if path="$(command -v "${tool}" 2>/dev/null)"; then
        printf '%b[%s] [%bOK%b] %-12s %s%b\n' "${DIM}" "$(ts)" "${GREEN}" "${RESET}${DIM}" "${tool}" "${path}" "${RESET}" >&2
        _log_plain OK "${tool}: ${path}"
        return 0
    else
        printf '%b[%s] [%bMISSING%b] %-12s not found in PATH%b\n' "${DIM}" "$(ts)" "${RED}" "${RESET}${DIM}" "${tool}" "${RESET}" >&2
        _log_plain MISSING "${tool} not found in PATH"
        return 1
    fi
}

check_dependencies() {
    local missing=0
    log_step "Checking external dependencies"

    # Core POSIX/GNU tools used by the driver and abundance parser.
    for t in awk sort find xargs; do
        need_tool "${t}" || missing=1
    done

    if should_run_step reads_count; then
        need_tool seqkit || missing=1
    fi
    if should_run_step extract; then
        need_tool bbduk.sh || missing=1
    fi
    if should_run_step align; then
        need_tool minimap2 || missing=1
    fi

    if [[ "${missing}" -ne 0 ]]; then
        die "Dependency check failed. Load/install the missing tools before running this pipeline."
    fi

    log "Dependency check passed"
}

# ==============================================================================
# Validation and initialization
# ==============================================================================
init_paths() {
    MANIFEST="${OUTDIR}/00_manifest/samples.tsv"
    LOG_DIR="${OUTDIR}/00_logs"
    STATUS_DIR="${OUTDIR}/99_checkpoints"
    CLEAN_DIR="${OUTDIR}/01_reads_count"
    CLEAN_OUT="${CLEAN_DIR}/clean_reads.tsv"
    CLEAN_TMP_DIR="${CLEAN_DIR}/tmp_counts"
    MARKER_DIR="${OUTDIR}/02_marker_reads"
    ALIGN_DIR="${OUTDIR}/03_align"
    ABUND_DIR="${OUTDIR}/04_abundance"
    SCRIPT_DIR="${OUTDIR}/00_scripts"

    mkdir -p \
        "${OUTDIR}" "${LOG_DIR}" "${STATUS_DIR}" "${CLEAN_DIR}" "${CLEAN_TMP_DIR}" \
        "${MARKER_DIR}/16S" "${MARKER_DIR}/18S" "${MARKER_DIR}/ITS" "${MARKER_DIR}/stats" "${MARKER_DIR}/logs" "${MARKER_DIR}/tmp" \
        "${ALIGN_DIR}/16S" "${ALIGN_DIR}/18S" "${ALIGN_DIR}/ITS" "${ALIGN_DIR}/logs" "${ALIGN_DIR}/tmp" \
        "${ABUND_DIR}" "${SCRIPT_DIR}" \
        "${STATUS_DIR}/reads_count" "${STATUS_DIR}/extract" "${STATUS_DIR}/align" "${STATUS_DIR}/abundance"

    MAIN_LOG="${LOG_DIR}/marker_count.$(date '+%Y%m%d_%H%M%S').log"
    touch "${MAIN_LOG}"
}

validate_args() {
    is_positive_int "${JOBS}" || die "--jobs must be a positive integer: ${JOBS}"
    is_positive_int "${THREADS_PER_SAMPLE}" || die "--threads-per-sample must be a positive integer: ${THREADS_PER_SAMPLE}"
    is_positive_int "${SEQKIT_THREADS}" || die "--seqkit-threads must be a positive integer: ${SEQKIT_THREADS}"

    for v in K_16S K_18S K_ITS MKH_16S MKH_18S MKH_ITS MINIMAP2_N MIN_ALN_LEN_16S MIN_ALN_LEN_18S MIN_ALN_LEN_ITS TOP_ALN_LEN_DIFF; do
        is_positive_int "${!v}" || die "--${v,,} must be a positive integer: ${!v}"
    done
    for v in HDIST_16S HDIST_18S HDIST_ITS MINK_16S MINK_18S MINK_ITS MIN_MAPQ; do
        is_nonnegative_int "${!v}" || die "--${v,,} must be a non-negative integer: ${!v}"
    done

    STEPS="$(normalize_steps "${STEPS}")"

    MARKERS="$(normalize_markers "${MARKERS}")"

    [[ "${RANK}" =~ ^(domain|phylum|class|order|family|genus|species|all)$ ]] || \
        die "--rank must be one of domain,phylum,class,order,family,genus,species,all. Got: ${RANK}"

    if [[ -n "${INPUT}" ]]; then
        [[ -s "${INPUT}" ]] || die "Input TSV not found or empty: ${INPUT}"
        [[ -z "${SAMPLE_ID}${R1_PATH}${R2_PATH}" ]] || die "Use either --input or --sample-id/--r1/--r2, not both."
    else
        [[ -n "${SAMPLE_ID}" && -n "${R1_PATH}" && -n "${R2_PATH}" ]] || \
            die "Missing input. Use --input data_path.tsv, or --sample-id with --r1 and --r2."
        [[ -s "${R1_PATH}" ]] || die "R1 not found or empty: ${R1_PATH}"
        [[ -s "${R2_PATH}" ]] || die "R2 not found or empty: ${R2_PATH}"
    fi

    if should_run_step reads_count; then
        command -v seqkit >/dev/null 2>&1 || die "seqkit not found in PATH."
    fi

    if should_run_step extract; then
        command -v bbduk.sh >/dev/null 2>&1 || die "bbduk.sh not found in PATH."
        if list_has_marker 16S; then [[ -s "${REF_16S}" ]] || die "Missing --ref-16s for 16S extraction."; fi
        if list_has_marker 18S; then [[ -s "${REF_18S}" ]] || die "Missing --ref-18s for 18S extraction."; fi
        if list_has_marker ITS; then [[ -s "${REF_ITS}" ]] || die "Missing --ref-its for ITS extraction."; fi
    fi

    if should_run_step align; then
        command -v minimap2 >/dev/null 2>&1 || die "minimap2 not found in PATH."
        if list_has_marker 16S; then [[ -s "${INDEX_16S}" ]] || die "Missing --index-16s for 16S alignment."; fi
        if list_has_marker 18S; then [[ -s "${INDEX_18S}" ]] || die "Missing --index-18s for 18S alignment."; fi
        if list_has_marker ITS; then [[ -s "${INDEX_ITS}" ]] || die "Missing --index-its for ITS alignment."; fi
    fi

    if should_run_step abundance; then
        command -v awk >/dev/null 2>&1 || die "awk not found in PATH."
        [[ -s "${TAXONOMY}" ]] || die "Missing --taxonomy for abundance step."
    fi
}

prepare_manifest() {
    mkdir -p "$(dirname "${MANIFEST}")"

    if [[ -n "${INPUT}" ]]; then
        awk -F '\t' -v OFS='\t' '
            function trim_cr(s) { sub(/\r$/, "", s); return s }
            function err(msg) { printf "[%s] [ERROR] %s\n", strftime("%Y-%m-%d %H:%M:%S"), msg > "/dev/stderr" }
            function value(name,    v) {
                if (!(name in col)) return "NA"
                v = $(col[name])
                v = trim_cr(v)
                return (v == "" ? "NA" : v)
            }
            function is_standard(name) {
                return (name == "sample_id" || name == "year" || name == "month" || name == "depth" || name == "r1_path" || name == "r2_path")
            }
            BEGIN {
                req[1] = "sample_id"
                req[2] = "r1_path"
                req[3] = "r2_path"
            }
            NR == 1 {
                for (i = 1; i <= NF; i++) {
                    name = trim_cr($i)
                    if (name == "") continue
                    if (name in col) {
                        err("Duplicate input column name in header: " name)
                        exit 2
                    }
                    col[name] = i
                    header[++header_n] = name
                }

                missing = ""
                for (i = 1; i <= 3; i++) {
                    if (!(req[i] in col)) {
                        missing = missing (missing == "" ? "" : ",") req[i]
                    }
                }
                if (missing != "") {
                    err("Input TSV must contain a header with required columns: sample_id,r1_path,r2_path")
                    err("Missing required column(s): " missing)
                    exit 2
                }

                extra_n = 0
                for (i = 1; i <= header_n; i++) {
                    name = header[i]
                    if (!is_standard(name)) extra[++extra_n] = name
                }

                printf "sample_id\tyear\tmonth\tdepth\tr1_path\tr2_path"
                for (i = 1; i <= extra_n; i++) printf "\t%s", extra[i]
                printf "\n"
                next
            }
            {
                trim_cr($0)
                if ($0 == "") next
                sid = value("sample_id")
                year = value("year")
                month = value("month")
                depth = value("depth")
                r1 = value("r1_path")
                r2 = value("r2_path")
                if (sid == "" || sid == "NA") {
                    err("Empty sample_id at input line " NR)
                    exit 2
                }
                if (sid ~ /\//) {
                    err("sample_id contains slash at input line " NR ": " sid)
                    exit 2
                }
                if (r1 == "" || r1 == "NA" || r2 == "" || r2 == "NA") {
                    err("Empty r1_path/r2_path at input line " NR)
                    exit 2
                }
                printf "%s\t%s\t%s\t%s\t%s\t%s", sid, year, month, depth, r1, r2
                for (i = 1; i <= extra_n; i++) printf "\t%s", value(extra[i])
                printf "\n"
            }
        ' "${INPUT}" > "${MANIFEST}.tmp"
        mv -f "${MANIFEST}.tmp" "${MANIFEST}"
        log "Input TSV columns were matched by header names. Required: sample_id, r1_path, r2_path. Optional metadata columns were retained."
    else
        printf 'sample_id\tyear\tmonth\tdepth\tr1_path\tr2_path\n' > "${MANIFEST}"
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
            "${SAMPLE_ID}" "${YEAR}" "${MONTH}" "${DEPTH}" "${R1_PATH}" "${R2_PATH}" >> "${MANIFEST}"
    fi

    local n
    n=$(awk -F '\t' 'NR>1 && NF>=6 {n++} END{print n+0}' "${MANIFEST}")
    [[ "${n}" -gt 0 ]] || die "No valid sample rows found in manifest: ${MANIFEST}"

    awk -F '\t' '
        function err(msg) { printf "[%s] [ERROR] %s\n", strftime("%Y-%m-%d %H:%M:%S"), msg > "/dev/stderr" }
        NR>1 {
            if ($1 == "") { err("Empty sample_id at manifest line " NR); exit 2 }
            if ($1 ~ /\//) { err("sample_id contains slash at manifest line " NR ": " $1); exit 2 }
            if ($5 == "" || $6 == "") { err("Empty R1/R2 path at manifest line " NR); exit 2 }
        }' "${MANIFEST}"

    log "Manifest written: ${MANIFEST}"
    log "Sample count: ${n}"
}

print_run_plan() {
    cat >&2 <<EOF_PLAN

${BOLD}Run plan${RESET}
  outdir             : ${OUTDIR}
  manifest           : ${MANIFEST}
  steps              : ${STEPS}
  markers            : ${MARKERS}
  rank               : ${RANK}
  jobs               : ${JOBS}
  threads/sample     : ${THREADS_PER_SAMPLE}
  seqkit threads     : ${SEQKIT_THREADS}
  resume             : ${RESUME}
  force              : ${FORCE}
  reads count table  : ${CLEAN_OUT}
  marker read dir    : ${MARKER_DIR}
  alignment dir      : ${ALIGN_DIR}
  abundance dir      : ${ABUND_DIR}
  main log           : ${MAIN_LOG}
EOF_PLAN
}

# ==============================================================================
# Step 1: reads count
# ==============================================================================
clean_count_one_sample() {
    local row_no="$1" sample_id="$2" year="$3" month="$4" depth="$5" r1_path="$6" r2_path="$7"
    local out="${CLEAN_TMP_DIR}/${row_no}.${sample_id}.clean_reads.tsv"
    local done="${STATUS_DIR}/reads_count/${sample_id}.done"
    local log_file="${LOG_DIR}/${sample_id}.reads_count.log"

    if [[ "${FORCE}" -eq 0 && "${RESUME}" -eq 1 && -s "${out}" && -s "${done}" ]]; then
        echo "[$(ts)] [INFO] reads_count skip: ${sample_id}" >> "${log_file}"
        return 0
    fi

    [[ -s "${r1_path}" ]] || { echo "[$(ts)] [ERROR] Missing R1: ${r1_path}" >> "${log_file}"; return 1; }
    [[ -s "${r2_path}" ]] || { echo "[$(ts)] [ERROR] Missing R2: ${r2_path}" >> "${log_file}"; return 1; }

    local read_pairs clean_reads_total
    read_pairs=$(seqkit stats -j "${SEQKIT_THREADS}" -T "${r1_path}" 2>>"${log_file}" | awk 'NR==2 {gsub(/,/, "", $4); print $4}')
    [[ "${read_pairs}" =~ ^[0-9]+$ ]] || { echo "[$(ts)] [ERROR] Failed to parse seqkit count for ${r1_path}" >> "${log_file}"; return 1; }
    clean_reads_total=$(( read_pairs * 2 ))

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${sample_id}" "${year}" "${month}" "${depth}" "${r1_path}" "${r2_path}" "${read_pairs}" "${clean_reads_total}" > "${out}.tmp"
    mv -f "${out}.tmp" "${out}"
    printf 'sample_id=%s\tread_pairs=%s\tclean_reads_total=%s\n' "${sample_id}" "${read_pairs}" "${clean_reads_total}" > "${done}"
    echo "[$(ts)] [INFO] reads_count done: ${sample_id} read_pairs=${read_pairs}" >> "${log_file}"
}

run_clean_counts() {
    local sample_count
    sample_count=$(awk -F '\t' 'NR>1 && NF>=6 {n++} END{print n+0}' "${MANIFEST}")

    if [[ "${FORCE}" -eq 0 && "${RESUME}" -eq 1 && -s "${CLEAN_OUT}" ]]; then
        local existing_rows
        existing_rows=$(awk -F '\t' 'NR>1 {n++} END{print n+0}' "${CLEAN_OUT}")
        if [[ "${existing_rows}" -eq "${sample_count}" ]]; then
            log "Reads-count table exists and row count matches; skipping reads_count step: ${CLEAN_OUT}"
            return 0
        fi
    fi

    log_step "Step reads_count: counting reads from clean FASTQ files with seqkit"
    rm -f "${CLEAN_TMP_DIR}"/*.clean_reads.tsv 2>/dev/null || true

    export CLEAN_TMP_DIR STATUS_DIR LOG_DIR SEQKIT_THREADS FORCE RESUME
    export -f ts clean_count_one_sample

    awk -F '\t' 'NR>1 && NF>=6 {
        row++
        printf "%s%c%s%c%s%c%s%c%s%c%s%c%s%c", row,0,$1,0,$2,0,$3,0,$4,0,$5,0,$6,0
    }' "${MANIFEST}" \
        | xargs -0 -r -P "${JOBS}" -n 7 bash -euo pipefail -c 'clean_count_one_sample "$@"' _

    printf 'sample_id\tyear\tmonth\tdepth\tr1_path\tr2_path\tread_pairs\tclean_reads_total\n' > "${CLEAN_OUT}.tmp"
    local i f
    for i in $(seq 1 "${sample_count}"); do
        f=$(find "${CLEAN_TMP_DIR}" -maxdepth 1 -type f -name "${i}.*.clean_reads.tsv" | head -n 1)
        [[ -s "${f}" ]] || die "Missing clean-count temporary result for row ${i}. Check logs in ${LOG_DIR}."
        cat "${f}" >> "${CLEAN_OUT}.tmp"
    done
    mv -f "${CLEAN_OUT}.tmp" "${CLEAN_OUT}"
    log "Written clean counts: ${CLEAN_OUT}"
}

# ==============================================================================
# Step 2: BBDuk marker extraction
# ==============================================================================
marker_ref() {
    case "$1" in
        16S) echo "${REF_16S}" ;;
        18S) echo "${REF_18S}" ;;
        ITS) echo "${REF_ITS}" ;;
    esac
}
marker_index() {
    case "$1" in
        16S) echo "${INDEX_16S}" ;;
        18S) echo "${INDEX_18S}" ;;
        ITS) echo "${INDEX_ITS}" ;;
    esac
}
marker_k() {
    case "$1" in
        16S) echo "${K_16S}" ;;
        18S) echo "${K_18S}" ;;
        ITS) echo "${K_ITS}" ;;
    esac
}
marker_hdist() {
    case "$1" in
        16S) echo "${HDIST_16S}" ;;
        18S) echo "${HDIST_18S}" ;;
        ITS) echo "${HDIST_ITS}" ;;
    esac
}
marker_mkh() {
    case "$1" in
        16S) echo "${MKH_16S}" ;;
        18S) echo "${MKH_18S}" ;;
        ITS) echo "${MKH_ITS}" ;;
    esac
}
marker_mink() {
    case "$1" in
        16S) echo "${MINK_16S}" ;;
        18S) echo "${MINK_18S}" ;;
        ITS) echo "${MINK_ITS}" ;;
    esac
}

extract_one_marker() {
    local sample_id="$1" r1_path="$2" r2_path="$3" marker="$4"
    local ref k hdist mkh mink
    ref="$(marker_ref "${marker}")"
    k="$(marker_k "${marker}")"
    hdist="$(marker_hdist "${marker}")"
    mkh="$(marker_mkh "${marker}")"
    mink="$(marker_mink "${marker}")"

    local out_r1="${MARKER_DIR}/${marker}/${sample_id}.${marker}.R1.fq.gz"
    local out_r2="${MARKER_DIR}/${marker}/${sample_id}.${marker}.R2.fq.gz"
    local stat="${MARKER_DIR}/stats/${sample_id}.${marker}.bbduk.stats.txt"
    local log_file="${MARKER_DIR}/logs/${sample_id}.${marker}.bbduk.log"
    local done="${STATUS_DIR}/extract/${sample_id}.${marker}.done"
    local tag="${sample_id}.${marker}.$$"
    local tmp_r1="${MARKER_DIR}/tmp/${tag}.R1.fq.gz.tmp"
    local tmp_r2="${MARKER_DIR}/tmp/${tag}.R2.fq.gz.tmp"
    local tmp_stat="${MARKER_DIR}/tmp/${tag}.stats.tmp"

    if [[ "${FORCE}" -eq 0 && "${RESUME}" -eq 1 && -s "${out_r1}" && -s "${out_r2}" && -s "${stat}" && -s "${done}" ]]; then
        echo "[$(ts)] [INFO] extract skip: ${sample_id} ${marker}" >> "${log_file}"
        return 0
    fi

    rm -f "${tmp_r1}" "${tmp_r2}" "${tmp_stat}"
    local cmd=(bbduk.sh "-Xmx${BBMAP_MEM}" "in1=${r1_path}" "in2=${r2_path}" "outm1=${tmp_r1}" "outm2=${tmp_r2}" "ref=${ref}" "k=${k}" "hdist=${hdist}" "mkh=${mkh}" "t=${THREADS_PER_SAMPLE}" "ordered=t" "stats=${tmp_stat}")
    if [[ "${mink}" -gt 0 ]]; then
        cmd+=("mink=${mink}")
    fi

    echo "[$(ts)] [INFO] extract start: ${sample_id} ${marker}" > "${log_file}"
    { printf '[%s] [CMD]' "$(ts)"; printf ' %q' "${cmd[@]}"; printf '\n'; } >> "${log_file}"

    "${cmd[@]}" >> "${log_file}" 2>&1

    mv -f "${tmp_r1}" "${out_r1}"
    mv -f "${tmp_r2}" "${out_r2}"
    mv -f "${tmp_stat}" "${stat}"
    printf 'sample_id=%s\tmarker=%s\tout_r1=%s\tout_r2=%s\n' "${sample_id}" "${marker}" "${out_r1}" "${out_r2}" > "${done}"
    echo "[$(ts)] [INFO] extract done: ${sample_id} ${marker}" >> "${log_file}"
}

extract_one_sample() {
    local sample_id="$1" year="$2" month="$3" depth="$4" r1_path="$5" r2_path="$6" marker
    [[ -s "${r1_path}" ]] || { echo "[$(ts)] [ERROR] Missing R1: ${r1_path}" >&2; return 1; }
    [[ -s "${r2_path}" ]] || { echo "[$(ts)] [ERROR] Missing R2: ${r2_path}" >&2; return 1; }
    for marker in $(for_each_marker); do
        extract_one_marker "${sample_id}" "${r1_path}" "${r2_path}" "${marker}"
    done
}

run_extract() {
    log_step "Step extract: extracting marker reads with BBDuk"
    export REF_16S REF_18S REF_ITS K_16S HDIST_16S MKH_16S MINK_16S K_18S HDIST_18S MKH_18S MINK_18S K_ITS HDIST_ITS MKH_ITS MINK_ITS
    export MARKER_DIR STATUS_DIR THREADS_PER_SAMPLE BBMAP_MEM FORCE RESUME MARKERS
    export -f ts for_each_marker marker_ref marker_k marker_hdist marker_mkh marker_mink extract_one_marker extract_one_sample

    awk -F '\t' 'NR>1 && NF>=6 { printf "%s%c%s%c%s%c%s%c%s%c%s%c", $1,0,$2,0,$3,0,$4,0,$5,0,$6,0 }' "${MANIFEST}" \
        | xargs -0 -r -P "${JOBS}" -n 6 bash -euo pipefail -c 'extract_one_sample "$@"' _
    log "Step extract finished"
}

# ==============================================================================
# Step 3: minimap2 alignment
# ==============================================================================
align_one_marker() {
    local sample_id="$1" marker="$2"
    local index
    index="$(marker_index "${marker}")"

    local in_r1="${MARKER_DIR}/${marker}/${sample_id}.${marker}.R1.fq.gz"
    local in_r2="${MARKER_DIR}/${marker}/${sample_id}.${marker}.R2.fq.gz"
    local out_r1="${ALIGN_DIR}/${marker}/${sample_id}.${marker}.R1.paf"
    local out_r2="${ALIGN_DIR}/${marker}/${sample_id}.${marker}.R2.paf"
    local done="${STATUS_DIR}/align/${sample_id}.${marker}.done"

    [[ -e "${in_r1}" ]] || { echo "[$(ts)] [ERROR] Missing marker R1: ${in_r1}" >&2; return 1; }
    [[ -e "${in_r2}" ]] || { echo "[$(ts)] [ERROR] Missing marker R2: ${in_r2}" >&2; return 1; }

    if [[ "${FORCE}" -eq 0 && "${RESUME}" -eq 1 && -e "${out_r1}" && -e "${out_r2}" && -s "${done}" ]]; then
        echo "[$(ts)] [INFO] align skip: ${sample_id} ${marker}" >> "${ALIGN_DIR}/logs/${sample_id}.${marker}.minimap2.log"
        return 0
    fi

    local fq mate out tmp log_file
    for mate in R1 R2; do
        if [[ "${mate}" == "R1" ]]; then fq="${in_r1}"; out="${out_r1}"; else fq="${in_r2}"; out="${out_r2}"; fi
        tmp="${ALIGN_DIR}/tmp/${sample_id}.${marker}.${mate}.$$.paf.tmp"
        log_file="${ALIGN_DIR}/logs/${sample_id}.${marker}.${mate}.minimap2.log"

        if [[ "${FORCE}" -eq 0 && "${RESUME}" -eq 1 && -e "${out}" ]]; then
            echo "[$(ts)] [INFO] align skip existing: ${out}" >> "${log_file}"
            continue
        fi

        echo "[$(ts)] [INFO] align start: ${sample_id} ${marker} ${mate}" > "${log_file}"
        echo "[$(ts)] [CMD] minimap2 -x ${MINIMAP2_PRESET} -N ${MINIMAP2_N} -t ${THREADS_PER_SAMPLE} ${index} ${fq}" >> "${log_file}"
        minimap2 -x "${MINIMAP2_PRESET}" -N "${MINIMAP2_N}" -t "${THREADS_PER_SAMPLE}" "${index}" "${fq}" > "${tmp}" 2>> "${log_file}"
        mv -f "${tmp}" "${out}"
        echo "[$(ts)] [INFO] align done: ${out}" >> "${log_file}"
    done

    printf 'sample_id=%s\tmarker=%s\tout_r1=%s\tout_r2=%s\n' "${sample_id}" "${marker}" "${out_r1}" "${out_r2}" > "${done}"
}

align_one_sample() {
    local sample_id="$1" year="$2" month="$3" depth="$4" r1_path="$5" r2_path="$6" marker
    for marker in $(for_each_marker); do
        align_one_marker "${sample_id}" "${marker}"
    done
}

run_align() {
    log_step "Step align: aligning marker reads with minimap2"
    export INDEX_16S INDEX_18S INDEX_ITS MARKER_DIR ALIGN_DIR STATUS_DIR THREADS_PER_SAMPLE FORCE RESUME MARKERS MINIMAP2_PRESET MINIMAP2_N
    export -f ts for_each_marker marker_index align_one_marker align_one_sample

    awk -F '\t' 'NR>1 && NF>=6 { printf "%s%c%s%c%s%c%s%c%s%c%s%c", $1,0,$2,0,$3,0,$4,0,$5,0,$6,0 }' "${MANIFEST}" \
        | xargs -0 -r -P "${JOBS}" -n 6 bash -euo pipefail -c 'align_one_sample "$@"' _
    log "Step align finished"
}

# ==============================================================================
# Step 4: abundance by shell/awk
# ==============================================================================
write_embedded_abundance_script() {
    local helper="${SCRIPT_DIR}/marker_count_abundance.awk.sh"
    cat > "${helper}" <<'AWK_HELPER'
#!/usr/bin/env bash
set -euo pipefail

DATA_PATH=""
CLEAN_COUNTS=""
TAXONOMY=""
ALIGN_DIR=""
MARKERS="16S,ITS"
RANK="genus"
OUTPUT_PREFIX=""
MIN_IDENTITY_16S="0.97"
MIN_IDENTITY_18S="0.97"
MIN_IDENTITY_ITS="0.95"
MIN_ALN_LEN_16S=80
MIN_ALN_LEN_18S=80
MIN_ALN_LEN_ITS=80
MIN_QCOV_16S="0.60"
MIN_QCOV_18S="0.60"
MIN_QCOV_ITS="0.60"
MIN_MAPQ=0
TOP_IDENTITY_DIFF="0.010"
TOP_ALN_LEN_DIFF=10

ats() { date '+%Y-%m-%d %H:%M:%S'; }
alog() {
    local level="$1"; shift
    printf '[%s] [%s] %s\n' "$(ats)" "${level}" "$*" >&2
}

require_value() {
    local opt="$1" value="${2:-}"
    [[ -n "${value}" && "${value}" != -* ]] || { alog ERROR "Option ${opt} requires a value."; exit 2; }
}

while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --data-path) require_value "$1" "${2:-}"; DATA_PATH="$2"; shift 2 ;;
        --clean-counts) require_value "$1" "${2:-}"; CLEAN_COUNTS="$2"; shift 2 ;;
        --taxonomy) require_value "$1" "${2:-}"; TAXONOMY="$2"; shift 2 ;;
        --align-dir) require_value "$1" "${2:-}"; ALIGN_DIR="$2"; shift 2 ;;
        --markers) require_value "$1" "${2:-}"; MARKERS="$2"; shift 2 ;;
        --rank) require_value "$1" "${2:-}"; RANK="$2"; shift 2 ;;
        --output-prefix) require_value "$1" "${2:-}"; OUTPUT_PREFIX="$2"; shift 2 ;;
        --min-identity-16s) require_value "$1" "${2:-}"; MIN_IDENTITY_16S="$2"; shift 2 ;;
        --min-identity-18s) require_value "$1" "${2:-}"; MIN_IDENTITY_18S="$2"; shift 2 ;;
        --min-identity-its) require_value "$1" "${2:-}"; MIN_IDENTITY_ITS="$2"; shift 2 ;;
        --min-aln-len-16s) require_value "$1" "${2:-}"; MIN_ALN_LEN_16S="$2"; shift 2 ;;
        --min-aln-len-18s) require_value "$1" "${2:-}"; MIN_ALN_LEN_18S="$2"; shift 2 ;;
        --min-aln-len-its) require_value "$1" "${2:-}"; MIN_ALN_LEN_ITS="$2"; shift 2 ;;
        --min-qcov-16s) require_value "$1" "${2:-}"; MIN_QCOV_16S="$2"; shift 2 ;;
        --min-qcov-18s) require_value "$1" "${2:-}"; MIN_QCOV_18S="$2"; shift 2 ;;
        --min-qcov-its) require_value "$1" "${2:-}"; MIN_QCOV_ITS="$2"; shift 2 ;;
        --min-mapq) require_value "$1" "${2:-}"; MIN_MAPQ="$2"; shift 2 ;;
        --top-identity-diff) require_value "$1" "${2:-}"; TOP_IDENTITY_DIFF="$2"; shift 2 ;;
        --top-aln-len-diff) require_value "$1" "${2:-}"; TOP_ALN_LEN_DIFF="$2"; shift 2 ;;
        *) alog ERROR "Unknown argument for abundance parser: $1"; exit 2 ;;
    esac
done

[[ -s "${DATA_PATH}" ]] || { alog ERROR "Missing data-path: ${DATA_PATH}"; exit 1; }
[[ -s "${CLEAN_COUNTS}" ]] || { alog ERROR "Missing clean-counts: ${CLEAN_COUNTS}"; exit 1; }
[[ -s "${TAXONOMY}" ]] || { alog ERROR "Missing taxonomy: ${TAXONOMY}"; exit 1; }
[[ -n "${OUTPUT_PREFIX}" ]] || { alog ERROR "Missing --output-prefix"; exit 1; }
mkdir -p "$(dirname "${OUTPUT_PREFIX}")"

run_rank() {
    local rank="$1"
    local long_out="${OUTPUT_PREFIX}.${rank}.long.tsv"
    local matrix_out="${OUTPUT_PREFIX}.${rank}.matrix.tsv"
    local domain_out="${OUTPUT_PREFIX}.${rank}.domain_total.tsv"
    local stats_out="${OUTPUT_PREFIX}.${rank}.assignment_stats.tsv"

    awk \
        -v DATA_PATH="${DATA_PATH}" \
        -v CLEAN_COUNTS="${CLEAN_COUNTS}" \
        -v TAXONOMY="${TAXONOMY}" \
        -v ALIGN_DIR="${ALIGN_DIR}" \
        -v MARKERS="${MARKERS}" \
        -v RANK="${rank}" \
        -v LONG_OUT="${long_out}" \
        -v MATRIX_OUT="${matrix_out}" \
        -v DOMAIN_OUT="${domain_out}" \
        -v STATS_OUT="${stats_out}" \
        -v MIN_IDENTITY_16S="${MIN_IDENTITY_16S}" \
        -v MIN_IDENTITY_18S="${MIN_IDENTITY_18S}" \
        -v MIN_IDENTITY_ITS="${MIN_IDENTITY_ITS}" \
        -v MIN_ALN_LEN_16S="${MIN_ALN_LEN_16S}" \
        -v MIN_ALN_LEN_18S="${MIN_ALN_LEN_18S}" \
        -v MIN_ALN_LEN_ITS="${MIN_ALN_LEN_ITS}" \
        -v MIN_QCOV_16S="${MIN_QCOV_16S}" \
        -v MIN_QCOV_18S="${MIN_QCOV_18S}" \
        -v MIN_QCOV_ITS="${MIN_QCOV_ITS}" \
        -v MIN_MAPQ="${MIN_MAPQ}" \
        -v TOP_IDENTITY_DIFF="${TOP_IDENTITY_DIFF}" \
        -v TOP_ALN_LEN_DIFF="${TOP_ALN_LEN_DIFF}" '
function trim_cr(s) { sub(/\r$/, "", s); return s }
function clean_value(s) {
    if (s == "" || s == "NA" || s == "nan" || s == "NaN") return "Unclassified"
    return s
}
function marker_min_identity(m) { return (m == "16S" ? MIN_IDENTITY_16S : (m == "18S" ? MIN_IDENTITY_18S : MIN_IDENTITY_ITS)) + 0 }
function marker_min_aln_len(m) { return (m == "16S" ? MIN_ALN_LEN_16S : (m == "18S" ? MIN_ALN_LEN_18S : MIN_ALN_LEN_ITS)) + 0 }
function marker_min_qcov(m) { return (m == "16S" ? MIN_QCOV_16S : (m == "18S" ? MIN_QCOV_18S : MIN_QCOV_ITS)) + 0 }
function read_tsv_header(path, col,    line,n,i,a) {
    delete last_header_name
    last_header_n = 0
    if ((getline line < path) <= 0) { printf "[ERROR] Cannot read header: %s\n", path > "/dev/stderr"; exit 1 }
    trim_cr(line)
    n = split(line, a, "\t")
    for (i=1; i<=n; i++) {
        a[i] = trim_cr(a[i])
        col[a[i]] = i
        last_header_name[i] = a[i]
    }
    last_header_n = n
    close(path)
}
function is_non_output_path_col(name) {
    return (name == "r1_path" || name == "r2_path")
}
function read_manifest(    line,n,a,path,col,i,name,sid,v) {
    path = DATA_PATH
    read_tsv_header(path, col)

    meta_n = 0
    for (i=1; i<=last_header_n; i++) {
        name = last_header_name[i]
        if (name == "" || name == "sample_id" || is_non_output_path_col(name)) continue
        meta_n++
        meta_name[meta_n] = name
    }

    getline line < path
    while ((getline line < path) > 0) {
        trim_cr(line)
        if (line == "") continue
        n = split(line, a, "\t")
        sample_n++
        sid = a[col["sample_id"]]
        sample_id[sample_n] = sid
        for (i=1; i<=meta_n; i++) {
            name = meta_name[i]
            v = ((name in col) ? a[col[name]] : "NA")
            v = trim_cr(v)
            if (v == "") v = "NA"
            sample_meta[sid SUBSEP name] = v
        }
    }
    close(path)
}
function read_clean(    line,n,a,path,col,sid) {
    path = CLEAN_COUNTS
    read_tsv_header(path, col)
    getline line < path
    while ((getline line < path) > 0) {
        trim_cr(line)
        if (line == "") continue
        n = split(line, a, "\t")
        sid = a[col["sample_id"]]
        clean_total[sid] = a[col["clean_reads_total"]] + 0
    }
    close(path)
}
function print_meta_header(out,    i) {
    printf "sample_id" > out
    for (i=1; i<=meta_n; i++) printf "\t%s", meta_name[i] >> out
}
function append_sample_meta(out, sid,    i,name,v) {
    printf "%s", sid >> out
    for (i=1; i<=meta_n; i++) {
        name = meta_name[i]
        v = sample_meta[sid SUBSEP name]
        if (v == "") v = "NA"
        printf "\t%s", v >> out
    }
}
function read_taxonomy(    line,n,a,path,col,ref) {
    path = TAXONOMY
    read_tsv_header(path, col)
    getline line < path
    if (!("ref_id" in col) || !("marker" in col) || !("domain" in col)) {
        printf "[ERROR] Taxonomy table requires ref_id, marker, domain columns.\n" > "/dev/stderr"; exit 1
    }
    while ((getline line < path) > 0) {
        trim_cr(line)
        if (line == "") continue
        n = split(line, a, "\t")
        ref = a[col["ref_id"]]
        tax_marker[ref] = a[col["marker"]]
        tax_domain[ref] = clean_value(a[col["domain"]])
        tax_rank[ref] = ((RANK in col) ? clean_value(a[col[RANK]]) : "Unclassified")
    }
    close(path)
}
function parse_markers(    n,i,a,m) {
    n = split(MARKERS, a, ",")
    for (i=1; i<=n; i++) {
        m = a[i]
        if (m != "") { marker_n++; marker_list[marker_n] = m }
    }
}
function paf_file(sid, marker, mate) {
    return ALIGN_DIR "/" marker "/" sid "." marker "." mate ".paf"
}
function reset_sample_marker(    k) {
    for (k in hit_count) delete hit_count[k]
    for (k in hit_id) delete hit_id[k]
    for (k in hit_aln) delete hit_aln[k]
    for (k in hit_matches) delete hit_matches[k]
    for (k in hit_mapq) delete hit_mapq[k]
    for (k in hit_domain) delete hit_domain[k]
    for (k in hit_taxon) delete hit_taxon[k]
    for (k in read_seen) delete read_seen[k]
    total_alignments = 0
    passed_alignments = 0
    reads_with_passed_hits = 0
    assigned_reads = 0
    ambiguous_reads = 0
}
function collect_paf(path, marker, mate,    line,a,n,qname,qlen,qstart,qend,tname,matches,aln_len,mapq,id,qcov,idx,key) {
    while ((getline line < path) > 0) {
        trim_cr(line)
        if (line == "") continue
        n = split(line, a, "\t")
        if (n < 12) continue
        total_alignments++
        qname = mate "::" a[1]
        sub(/ .*/, "", qname)
        qlen = a[2] + 0
        qstart = a[3] + 0
        qend = a[4] + 0
        tname = a[6]
        sub(/ .*/, "", tname)
        matches = a[10] + 0
        aln_len = a[11] + 0
        mapq = a[12] + 0
        if (!(tname in tax_marker)) continue
        if (tax_marker[tname] != marker) continue
        if (aln_len <= 0 || qlen <= 0) continue
        id = matches / aln_len
        qcov = (qend - qstart) / qlen
        if (id < marker_min_identity(marker)) continue
        if (aln_len < marker_min_aln_len(marker)) continue
        if (qcov < marker_min_qcov(marker)) continue
        if (mapq < MIN_MAPQ) continue
        passed_alignments++
        if (!(qname in read_seen)) { read_seen[qname] = 1; reads_with_passed_hits++ }
        hit_count[qname]++
        idx = hit_count[qname]
        key = qname SUBSEP idx
        hit_id[key] = id
        hit_aln[key] = aln_len
        hit_matches[key] = matches
        hit_mapq[key] = mapq
        hit_domain[key] = tax_domain[tname]
        hit_taxon[key] = tax_rank[tname]
    }
    close(path)
}
function assign_reads(marker, sid,    q,i,key,bid,baln,bmatches,bmapq,bestset,taxset,domset,ntax,ndom,one_tax,one_dom,ck) {
    for (q in hit_count) {
        bid = -1; baln = -1; bmatches = -1; bmapq = -1
        for (i=1; i<=hit_count[q]; i++) {
            key = q SUBSEP i
            if (hit_id[key] > bid || \
                (hit_id[key] == bid && hit_aln[key] > baln) || \
                (hit_id[key] == bid && hit_aln[key] == baln && hit_matches[key] > bmatches) || \
                (hit_id[key] == bid && hit_aln[key] == baln && hit_matches[key] == bmatches && hit_mapq[key] > bmapq)) {
                bid = hit_id[key]; baln = hit_aln[key]; bmatches = hit_matches[key]; bmapq = hit_mapq[key]
            }
        }
        delete taxset; delete domset; ntax = 0; ndom = 0
        for (i=1; i<=hit_count[q]; i++) {
            key = q SUBSEP i
            if ((bid - hit_id[key]) <= TOP_IDENTITY_DIFF && (baln - hit_aln[key]) <= TOP_ALN_LEN_DIFF) {
                if (!(hit_taxon[key] in taxset)) { taxset[hit_taxon[key]] = 1; ntax++; one_tax = hit_taxon[key] }
                if (!(hit_domain[key] in domset)) { domset[hit_domain[key]] = 1; ndom++; one_dom = hit_domain[key] }
            }
        }
        if (ntax == 1 && ndom == 1) {
            ck = one_dom SUBSEP marker SUBSEP one_tax
            count_reads[ck]++
            assigned_reads++
        } else {
            ambiguous_reads++
        }
    }
}
function process_sample_marker(sid, marker,    r1,r2,ck,a,domain,taxon,nreads,rpm,fid,dkey) {
    reset_sample_marker()
    r1 = paf_file(sid, marker, "R1")
    r2 = paf_file(sid, marker, "R2")
    collect_paf(r1, marker, "R1")
    collect_paf(r2, marker, "R2")
    assign_reads(marker, sid)

    append_sample_meta(STATS_OUT, sid)
    printf "\t%s\t%d\t%d\t%d\t%d\t%d\n", marker, total_alignments, passed_alignments, reads_with_passed_hits, assigned_reads, ambiguous_reads >> STATS_OUT

    for (ck in count_reads) {
        split(ck, a, SUBSEP)
        domain = a[1]; marker = a[2]; taxon = a[3]
        nreads = count_reads[ck]
        rpm = (clean_total[sid] > 0 ? nreads / clean_total[sid] * 1000000 : 0)
        append_sample_meta(LONG_OUT, sid)
        printf "\t%s\t%s\t%s\t%s\t%d\t%d\t%.10g\n", marker, domain, RANK, taxon, nreads, clean_total[sid], rpm >> LONG_OUT
        fid = domain "|" marker "|" RANK "__" taxon
        if (!(fid in feature_seen)) { feature_seen[fid] = 1; feature_n++; feature_list[feature_n] = fid }
        matrix_value[sid SUBSEP fid] += rpm
        dkey = sid SUBSEP domain SUBSEP marker
        domain_reads[dkey] += nreads
        domain_rpm[dkey] += rpm
        domain_clean[dkey] = clean_total[sid]
        delete count_reads[ck]
    }
}
BEGIN {
    OFS="\t"
    read_manifest()
    read_clean()
    read_taxonomy()
    parse_markers()

    print_meta_header(LONG_OUT)
    print "\tmarker\tdomain\trank\ttaxon\ttaxon_marker_reads\tclean_reads_total\tmarker_rpm" >> LONG_OUT
    print_meta_header(STATS_OUT)
    print "\tmarker\ttotal_alignments\tpassed_alignments\treads_with_passed_hits\tassigned_reads\tambiguous_reads" >> STATS_OUT

    for (si=1; si<=sample_n; si++) {
        sid = sample_id[si]
        if (!(sid in clean_total)) {
            printf "[ERROR] Sample %s missing in reads-count table.\n", sid > "/dev/stderr"; exit 1
        }
        for (mi=1; mi<=marker_n; mi++) process_sample_marker(sid, marker_list[mi])
    }

    print_meta_header(MATRIX_OUT)
    for (fi=1; fi<=feature_n; fi++) printf "\t%s", feature_list[fi] >> MATRIX_OUT
    printf "\n" >> MATRIX_OUT
    for (si=1; si<=sample_n; si++) {
        sid = sample_id[si]
        append_sample_meta(MATRIX_OUT, sid)
        for (fi=1; fi<=feature_n; fi++) {
            fid = feature_list[fi]
            printf "\t%.10g", matrix_value[sid SUBSEP fid] + 0 >> MATRIX_OUT
        }
        printf "\n" >> MATRIX_OUT
    }

    print_meta_header(DOMAIN_OUT)
    print "\tdomain\tmarker\tdomain_marker_reads\tclean_reads_total\tdomain_marker_rpm" >> DOMAIN_OUT
    for (dkey in domain_reads) {
        split(dkey, da, SUBSEP)
        sid = da[1]
        append_sample_meta(DOMAIN_OUT, sid)
        printf "\t%s\t%s\t%d\t%d\t%.10g\n", da[2], da[3], domain_reads[dkey], domain_clean[dkey], domain_rpm[dkey] >> DOMAIN_OUT
    }
}
' </dev/null

    alog INFO "Written: ${long_out}"
    alog INFO "Written: ${matrix_out}"
    alog INFO "Written: ${domain_out}"
    alog INFO "Written: ${stats_out}"
}

if [[ "${RANK}" == "all" ]]; then
    for r in domain phylum class order family genus species; do
        run_rank "${r}"
    done
else
    run_rank "${RANK}"
fi
AWK_HELPER
    chmod +x "${helper}"
    echo "${helper}"
}

run_abundance() {
    log_step "Step abundance: parsing PAF and calculating marker RPM"
    [[ -s "${CLEAN_OUT}" ]] || die "Reads-count table not found: ${CLEAN_OUT}. Run reads_count step first or provide existing file at this path."

    local helper
    helper=$(write_embedded_abundance_script)

    local out_prefix="${ABUND_DIR}/marker_rpm"
    local log_file="${ABUND_DIR}/abundance.log"
    local expected="${out_prefix}.${RANK}.long.tsv"
    if [[ "${RANK}" == "all" ]]; then
        expected="${out_prefix}.genus.long.tsv"
    fi

    if [[ "${FORCE}" -eq 0 && "${RESUME}" -eq 1 && -s "${expected}" ]]; then
        log "Abundance output exists; skipping abundance step: ${expected}"
        return 0
    fi

    {
        printf '[%s] [CMD] %q' "$(ts)" "${helper}"
        printf ' --data-path %q' "${MANIFEST}"
        printf ' --clean-counts %q' "${CLEAN_OUT}"
        printf ' --taxonomy %q' "${TAXONOMY}"
        printf ' --align-dir %q' "${ALIGN_DIR}"
        printf ' --markers %q' "${MARKERS}"
        printf ' --rank %q' "${RANK}"
        printf ' --output-prefix %q' "${out_prefix}"
        printf '\n'
    } > "${log_file}"

    "${helper}" \
        --data-path "${MANIFEST}" \
        --clean-counts "${CLEAN_OUT}" \
        --taxonomy "${TAXONOMY}" \
        --align-dir "${ALIGN_DIR}" \
        --markers "${MARKERS}" \
        --rank "${RANK}" \
        --output-prefix "${out_prefix}" \
        --min-identity-16s "${MIN_IDENTITY_16S}" \
        --min-identity-18s "${MIN_IDENTITY_18S}" \
        --min-identity-its "${MIN_IDENTITY_ITS}" \
        --min-aln-len-16s "${MIN_ALN_LEN_16S}" \
        --min-aln-len-18s "${MIN_ALN_LEN_18S}" \
        --min-aln-len-its "${MIN_ALN_LEN_ITS}" \
        --min-qcov-16s "${MIN_QCOV_16S}" \
        --min-qcov-18s "${MIN_QCOV_18S}" \
        --min-qcov-its "${MIN_QCOV_ITS}" \
        --min-mapq "${MIN_MAPQ}" \
        --top-identity-diff "${TOP_IDENTITY_DIFF}" \
        --top-aln-len-diff "${TOP_ALN_LEN_DIFF}" \
        >> "${log_file}" 2>&1

    cat "${log_file}" >&2
    printf 'rank=%s\toutput_prefix=%s\n' "${RANK}" "${out_prefix}" > "${STATUS_DIR}/abundance/abundance.done"
    log "Step abundance finished"
}

# ==============================================================================
# Summary and cleanup
# ==============================================================================
write_run_config() {
    local config="${LOG_DIR}/run_config.tsv"
    {
        printf 'key\tvalue\n'
        printf 'program\t%s\n' "${PROGRAM}"
        printf 'outdir\t%s\n' "${OUTDIR}"
        printf 'steps\t%s\n' "${STEPS}"
        printf 'markers\t%s\n' "${MARKERS}"
        printf 'jobs\t%s\n' "${JOBS}"
        printf 'threads_per_sample\t%s\n' "${THREADS_PER_SAMPLE}"
        printf 'seqkit_threads\t%s\n' "${SEQKIT_THREADS}"
        printf 'bbmap_mem\t%s\n' "${BBMAP_MEM}"
        printf 'ref_16s\t%s\n' "${REF_16S}"
        printf 'ref_18s\t%s\n' "${REF_18S}"
        printf 'ref_its\t%s\n' "${REF_ITS}"
        printf 'index_16s\t%s\n' "${INDEX_16S}"
        printf 'index_18s\t%s\n' "${INDEX_18S}"
        printf 'index_its\t%s\n' "${INDEX_ITS}"
        printf 'taxonomy\t%s\n' "${TAXONOMY}"
        printf 'rank\t%s\n' "${RANK}"
        printf 'force\t%s\n' "${FORCE}"
        printf 'resume\t%s\n' "${RESUME}"
    } > "${config}"
    log "Run config written: ${config}"
}

cleanup_tmp() {
    if [[ "${CLEAN_TMP}" -eq 1 ]]; then
        rm -rf "${CLEAN_TMP_DIR}" "${MARKER_DIR}/tmp" "${ALIGN_DIR}/tmp"
        log "Temporary directories removed"
    fi
}

main() {
    init_colors
    parse_args "$@"
    STEPS="$(normalize_steps "${STEPS}")"
    MARKERS="$(normalize_markers "${MARKERS}")"
    init_paths
    check_dependencies
    if [[ "${CHECK_DEPS_ONLY}" -eq 1 ]]; then
        exit 0
    fi
    validate_args
    prepare_manifest
    write_run_config
    print_run_plan

    if [[ "${DRY_RUN}" -eq 1 ]]; then
        log "Dry-run finished. No analysis was executed."
        exit 0
    fi

    if should_run_step reads_count; then run_clean_counts; fi
    if should_run_step extract; then run_extract; fi
    if should_run_step align; then run_align; fi
    if should_run_step abundance; then run_abundance; fi

    cleanup_tmp
    log "All requested steps finished successfully."
}

main "$@"

