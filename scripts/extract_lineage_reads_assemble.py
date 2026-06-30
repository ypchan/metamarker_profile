#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Extract marker reads assigned to a target lineage and assemble them with miniasm.

The selection logic mirrors metamarker_profile abundance assignment:
  1. parse mate-labeled PAF query IDs ending in /1 or /2;
  2. keep the best valid hit per pair_id + mate;
  3. arbitrate R1/R2 conflicts at pair level;
  4. extract reads assigned to the requested rank/lineage.

Outputs include mate-specific FASTQ, merged_segments.fq.gz, selected-read tables,
minimap2 all-vs-all overlaps, miniasm GFA, and FASTA contigs.
"""

from __future__ import annotations

import argparse
import ast
import csv
import gzip
import re
import shlex
import shutil
import subprocess
import sys
import time
from contextlib import nullcontext
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from textwrap import dedent
from typing import Callable, Dict, Iterator, List, Optional, Sequence, Tuple

from rich.console import Console
from rich.panel import Panel
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    SpinnerColumn,
    TaskProgressColumn,
    TextColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
)
from rich.table import Table
from rich_argparse import RawDescriptionRichHelpFormatter


RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
MARKER_ORDER = ["16S", "18S", "ITS"]
NA_GROUP = "NA"
RANK_PREFIX = {
    "domain": "d__",
    "phylum": "p__",
    "class": "c__",
    "order": "o__",
    "family": "f__",
    "genus": "g__",
    "species": "s__",
}
DEFAULT_MIN_IDENTITY = {"16S": 0.97, "18S": 0.97, "ITS": 0.95}
DEFAULT_MIN_ALN_LEN = {"16S": 80, "18S": 80, "ITS": 80}
DEFAULT_MIN_QCOV = {"16S": 0.60, "18S": 0.60, "ITS": 0.60}
DEFAULT_MIN_MAPQ = 0
DEFAULT_OUTDIR = "metamarker_profile_out"
DEFAULT_RANK = "species"
DEFAULT_MODE = "assigned"
DEFAULT_MATCH = "exact"
DEFAULT_ASSEMBLER = "miniasm"
DEFAULT_ASSEMBLY_THREADS = 4
DEFAULT_OVERLAP_PRESET = "ava-ont"
DEFAULT_JOBS = 4
DEFAULT_LOG_WIDTH_FRACTION = 0.5

RawDescriptionRichHelpFormatter.styles["argparse.prog"] = "bold magenta"
RawDescriptionRichHelpFormatter.styles["argparse.groups"] = "bold green"
RawDescriptionRichHelpFormatter.styles["argparse.args"] = "bold cyan"
RawDescriptionRichHelpFormatter.styles["argparse.metavar"] = "bold yellow"
RawDescriptionRichHelpFormatter.styles["argparse.help"] = "white"
RawDescriptionRichHelpFormatter.styles["argparse.text"] = "white"


@dataclass
class Hit:
    sample_id: str
    marker: str
    query_id: str
    pair_id: str
    mate: str
    strand: str
    target_id: str
    target_len: int
    matches: int
    aln_len: int
    mapq: int
    identity: float
    qcov: float
    taxonomy: Dict[str, str]


@dataclass
class SelectedRead:
    sample_id: str
    marker: str
    pair_id: str
    mate: str
    strand: str
    target_id: str
    identity: float
    aln_len: int
    mapq: int
    mate_euk_group: str
    mate_lineage: str
    mate_full_lineage: str
    assigned_euk_group: str
    assigned_lineage: str
    assigned_full_lineage: str
    assigned_target_id: str
    read_weight: int


@dataclass
class SelectResult:
    index: int
    sample_id: str
    marker: str
    paf_path: Path
    best_hits: int
    selected_reads: int
    elapsed: float
    selected: List[SelectedRead]


def open_text_auto(path: Path, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode, encoding=None if "b" in mode else "utf-8")
    return open(path, mode, encoding=None if "b" in mode else "utf-8")


def open_gzip_text(path: Path, mode: str = "wt"):
    return gzip.open(path, mode, encoding="utf-8")


def terminal_width(default: int = 120) -> int:
    return shutil.get_terminal_size((default, 24)).columns


def dynamic_log_width() -> int:
    return max(60, int(terminal_width() * DEFAULT_LOG_WIDTH_FRACTION))


def make_console() -> Console:
    return Console(stderr=True, width=dynamic_log_width())


def now() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def format_seconds(seconds: float) -> str:
    if seconds < 60:
        return f"{seconds:.2f}s"
    minutes, sec = divmod(seconds, 60)
    if minutes < 60:
        return f"{int(minutes)}m{sec:04.1f}s"
    hours, minutes = divmod(minutes, 60)
    return f"{int(hours)}h{int(minutes):02d}m{sec:04.1f}s"


def rate_text(count: int | float, seconds: float, unit: str) -> str:
    if seconds <= 0:
        return f"{count:g} {unit}/s"
    return f"{float(count) / seconds:.2f} {unit}/s"


def append_run_log(log_file: Path, level: str, message: str) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    with open(log_file, "a", encoding="utf-8") as fh:
        fh.write(f"[{now()}]\t{level}\t{message}\n")


def console_log(console: Console, log_file: Optional[Path], level: str, message: str, style: str = "white") -> None:
    console.print(f"[dim]{now()}[/dim] [{style}]{level}[/{style}] {message}")
    if log_file is not None:
        append_run_log(log_file, level, message)


def progress_bar(console: Console, label: str, total: int) -> Progress:
    return Progress(
        SpinnerColumn(),
        TextColumn(f"[bold cyan]{label}"),
        BarColumn(bar_width=None),
        TaskProgressColumn(),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        console=console,
    )


def outputs_ready(paths: Sequence[Path]) -> bool:
    return all(p.exists() for p in paths)


def stage_completed(log_file: Path, stage: str) -> bool:
    if not log_file.exists():
        return False
    needle = f"\tSTAGE_DONE\tstage={stage}"
    with open(log_file, "r", encoding="utf-8", errors="replace") as fh:
        return any(needle in line for line in fh)


def mark_stage_done(log_file: Path, stage: str, message: str = "") -> None:
    suffix = f"\t{message}" if message else ""
    append_run_log(log_file, "STAGE_DONE", f"stage={stage}{suffix}")


def resume_ready(log_file: Path, stage: str, outputs: Sequence[Path], force: bool, no_resume: bool) -> bool:
    return (not force) and (not no_resume) and stage_completed(log_file, stage) and outputs_ready(outputs)


def extraction_outputs(outdir: Path) -> List[Path]:
    return [
        outdir / "selected_reads.tsv",
        outdir / "summary.tsv",
        outdir / "reads.R1.fq.gz",
        outdir / "reads.R2.fq.gz",
        outdir / "merged_segments.fq.gz",
        outdir / "merged_segments.fa",
    ]


def assembly_outputs(outdir: Path) -> List[Path]:
    return [
        outdir / "miniasm.overlaps.paf",
        outdir / "miniasm.gfa",
        outdir / "miniasm_contigs.fa",
        outdir / "assembly_summary.tsv",
    ]


def write_status(path: Path, rows: Sequence[Tuple[str, str]]) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("key\tvalue\n")
        for key, value in rows:
            fh.write(f"{key}\t{value}\n")


def write_empty_extraction_outputs(outdir: Path) -> None:
    for name in ["reads.R1.fq.gz", "reads.R2.fq.gz", "merged_segments.fq.gz"]:
        with open_gzip_text(outdir / name, "wt"):
            pass
    with open(outdir / "merged_segments.fa", "w", encoding="utf-8"):
        pass


def validate_args(args: argparse.Namespace) -> None:
    if args.jobs < 1:
        raise ValueError("--jobs must be positive.")
    if args.assembly_threads < 1:
        raise ValueError("--assembly-threads must be positive.")
    if args.min_contig_len < 0:
        raise ValueError("--min-contig-len must be >= 0.")


def print_run_plan(
    console: Console,
    args: argparse.Namespace,
    output_dir: Path,
    log_file: Path,
    mode_name: str,
    markers: Optional[Sequence[str]] = None,
    sample_count: Optional[int] = None,
    segments_fastq: Optional[Path] = None,
) -> None:
    table = Table(show_header=False, box=None, pad_edge=False)
    table.add_column("parameter", style="bold cyan", no_wrap=True)
    table.add_column("value", style="white")
    table.add_row("mode", f"[green]{mode_name}[/green]")
    if segments_fastq is not None:
        table.add_row("segments_fastq", str(segments_fastq))
    if args.lineage:
        table.add_row("lineage", f"[yellow]{args.lineage}[/yellow]")
        table.add_row("rank/match", f"{args.rank} / {args.match}")
        table.add_row("selection", args.mode)
    if markers is not None:
        table.add_row("markers", ",".join(markers))
    if sample_count is not None:
        table.add_row("samples", str(sample_count))
    table.add_row("output_dir", str(output_dir))
    table.add_row("assembler", args.assembler)
    table.add_row("jobs", str(args.jobs))
    table.add_row("assembly_threads", str(args.assembly_threads))
    table.add_row("resume", "off" if args.no_resume else ("force recompute" if args.force else "on"))
    table.add_row("log", str(log_file))
    console.print(Panel(table, title="[bold green]Run parameters[/bold green]", border_style="green"))


def print_final_summary(console: Console, title: str, rows: Sequence[Tuple[str, object]]) -> None:
    table = Table(show_header=False, box=None, pad_edge=False)
    table.add_column("metric", style="bold cyan", no_wrap=True)
    table.add_column("value", style="white")
    for key, value in rows:
        table.add_row(key, str(value))
    console.print(Panel(table, title=f"[bold green]{title}[/bold green]", border_style="green"))


def parse_list(raw: Optional[str]) -> Optional[List[str]]:
    if raw is None or raw.strip() == "":
        return None
    return [x.strip() for x in raw.split(",") if x.strip()]


def parse_run_config(path: Path) -> Dict[str, object]:
    cfg: Dict[str, object] = {}
    if not path.exists():
        return cfg
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or "\t" not in line:
                continue
            key, value = line.split("\t", 1)
            try:
                cfg[key] = ast.literal_eval(value)
            except Exception:
                cfg[key] = value
    return cfg


def load_manifest(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Sample manifest not found: {path}")
    with open(path, "r", encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        raise ValueError(f"Sample manifest is empty: {path}")
    return rows


def load_taxonomy(path: Path) -> Dict[str, Dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Taxonomy table not found: {path}")
    tax: Dict[str, Dict[str, str]] = {}
    with open(path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"ref_id", "marker", "euk_group", *RANKS}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Taxonomy table missing columns: {','.join(sorted(missing))}")
        for row in reader:
            ref_id = (row.get("ref_id") or "").strip()
            if not ref_id:
                continue
            rec = {k: (row.get(k) or "").strip() for k in ["marker", *RANKS]}
            rec["euk_group"] = (row.get("euk_group") or "").strip()
            if not rec["euk_group"]:
                raise ValueError(f"Taxonomy table has empty euk_group for ref_id={ref_id}. Rebuild refs with scripts/build_ref_files.py.")
            if not rec["domain"]:
                rec["domain"] = "NA"
            for rank in RANKS:
                if rank != "domain" and not rec[rank]:
                    rec[rank] = "unidentified"
            tax[ref_id] = rec
    return tax


def get_threshold_map(cfg: Dict[str, object], key: str, default: Dict[str, object]) -> Dict[str, object]:
    value = cfg.get(key)
    if isinstance(value, dict):
        return {m: value.get(m, default[m]) for m in MARKER_ORDER}
    return dict(default)


def get_thresholds(cfg: Dict[str, object]) -> Tuple[Dict[str, float], Dict[str, int], Dict[str, float], int]:
    min_identity = {k: float(v) for k, v in get_threshold_map(cfg, "min_identity", DEFAULT_MIN_IDENTITY).items()}
    min_aln_len = {k: int(v) for k, v in get_threshold_map(cfg, "min_aln_len", DEFAULT_MIN_ALN_LEN).items()}
    min_qcov = {k: float(v) for k, v in get_threshold_map(cfg, "min_qcov", DEFAULT_MIN_QCOV).items()}
    min_mapq = int(cfg.get("min_mapq", DEFAULT_MIN_MAPQ))
    return min_identity, min_aln_len, min_qcov, min_mapq


def normalize_read_id(header_or_query: str) -> str:
    name = header_or_query.strip()
    if name.startswith("@"):
        name = name[1:]
    name = name.split()[0]
    name = re.sub(r"^R[12]::", "", name)
    name = re.sub(r"/[12]$", "", name)
    return name


def mate_from_query(query_id: str) -> str:
    if query_id.endswith("/1"):
        return "R1"
    if query_id.endswith("/2"):
        return "R2"
    return "unknown"


def better_hit(candidate: Hit, current: Optional[Hit]) -> bool:
    if current is None:
        return True
    for attr in ["identity", "aln_len", "matches", "mapq", "qcov", "target_len"]:
        a = getattr(candidate, attr)
        b = getattr(current, attr)
        if a > b:
            return True
        if a < b:
            return False
    return candidate.target_id < current.target_id


def better_winner(candidate: Hit, current: Optional[Hit]) -> bool:
    if current is None:
        return True
    for attr in ["identity", "aln_len", "matches", "mapq", "qcov", "target_len"]:
        a = getattr(candidate, attr)
        b = getattr(current, attr)
        if a > b:
            return True
        if a < b:
            return False
    if candidate.mate != current.mate:
        return candidate.mate < current.mate
    return candidate.target_id < current.target_id


def parse_paf_best_hits(
    paf_path: Path,
    sample_id: str,
    marker: str,
    taxonomy: Dict[str, Dict[str, str]],
    min_identity: Dict[str, float],
    min_aln_len: Dict[str, int],
    min_qcov: Dict[str, float],
    min_mapq: int,
) -> Dict[Tuple[str, str], Hit]:
    best: Dict[Tuple[str, str], Hit] = {}
    if not paf_path.exists() or paf_path.stat().st_size == 0:
        return best

    with open(paf_path, "r", encoding="utf-8", errors="replace") as fh:
        for line_no, line in enumerate(fh, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            try:
                query_id = fields[0]
                query_len = int(fields[1])
                query_start = int(fields[2])
                query_end = int(fields[3])
                strand = fields[4]
                target_id = fields[5]
                target_len = int(fields[6])
                matches = int(fields[9])
                aln_len = int(fields[10])
                mapq = int(fields[11])
            except ValueError:
                continue
            if query_len <= 0 or aln_len <= 0:
                continue

            rec = taxonomy.get(target_id)
            if not rec or rec.get("marker") != marker:
                continue

            identity = matches / float(aln_len)
            qcov = (query_end - query_start) / float(query_len)
            if identity < min_identity[marker]:
                continue
            if aln_len < min_aln_len[marker]:
                continue
            if qcov < min_qcov[marker]:
                continue
            if mapq < min_mapq:
                continue

            mate = mate_from_query(query_id)
            if mate == "unknown":
                raise ValueError(f"PAF query ID lacks /1 or /2 mate label at {paf_path}:{line_no}: {query_id}")
            pair_id = normalize_read_id(query_id)
            hit = Hit(
                sample_id=sample_id,
                marker=marker,
                query_id=query_id,
                pair_id=pair_id,
                mate=mate,
                strand=strand,
                target_id=target_id,
                target_len=target_len,
                matches=matches,
                aln_len=aln_len,
                mapq=mapq,
                identity=identity,
                qcov=qcov,
                taxonomy=rec,
            )
            key = (pair_id, mate)
            if better_hit(hit, best.get(key)):
                best[key] = hit
    return best


def build_matcher(lineage: str, match_mode: str, ignore_case: bool) -> Callable[[str], bool]:
    target = lineage.lower() if ignore_case else lineage
    if match_mode == "regex":
        flags = re.IGNORECASE if ignore_case else 0
        pattern = re.compile(lineage, flags=flags)
        return lambda value: bool(pattern.search(value or ""))
    if match_mode == "contains":
        return lambda value: target in ((value or "").lower() if ignore_case else (value or ""))
    return lambda value: ((value or "").lower() if ignore_case else (value or "")) == target


def lineage_ranks(rank: str) -> List[str]:
    return RANKS[: RANKS.index(rank) + 1]


def full_lineage(taxonomy: Dict[str, str], rank: str) -> str:
    parts = []
    for lineage_rank in lineage_ranks(rank):
        default_value = "NA" if lineage_rank == "domain" else "unidentified"
        parts.append(f"{RANK_PREFIX[lineage_rank]}{taxonomy.get(lineage_rank) or default_value}")
    return ";".join(parts)


def matches_lineage(hit: Hit, rank: str, lineage_match: Callable[[str], bool]) -> bool:
    return lineage_match(hit.taxonomy.get(rank, "unidentified")) or lineage_match(full_lineage(hit.taxonomy, rank))


def select_reads(
    hits: Dict[Tuple[str, str], Hit],
    rank: str,
    lineage_match: Callable[[str], bool],
    mode: str,
) -> List[SelectedRead]:
    by_pair: Dict[str, List[Hit]] = {}
    for hit in hits.values():
        by_pair.setdefault(hit.pair_id, []).append(hit)

    selected: List[SelectedRead] = []
    if mode == "mate-hit":
        for hit in hits.values():
            mate_lineage = hit.taxonomy.get(rank, "unidentified")
            mate_full_lineage = full_lineage(hit.taxonomy, rank)
            if not matches_lineage(hit, rank, lineage_match):
                continue
            selected.append(
                SelectedRead(
                    sample_id=hit.sample_id,
                    marker=hit.marker,
                    pair_id=hit.pair_id,
                    mate=hit.mate,
                    strand=hit.strand,
                    target_id=hit.target_id,
                    identity=hit.identity,
                    aln_len=hit.aln_len,
                    mapq=hit.mapq,
                    mate_euk_group=hit.taxonomy.get("euk_group", NA_GROUP),
                    mate_lineage=mate_lineage,
                    mate_full_lineage=mate_full_lineage,
                    assigned_euk_group=hit.taxonomy.get("euk_group", NA_GROUP),
                    assigned_lineage=mate_lineage,
                    assigned_full_lineage=mate_full_lineage,
                    assigned_target_id=hit.target_id,
                    read_weight=1,
                )
            )
        return selected

    for pair_id, pair_hits in by_pair.items():
        winner: Optional[Hit] = None
        for hit in pair_hits:
            if better_winner(hit, winner):
                winner = hit
        if winner is None:
            continue
        assigned_lineage = winner.taxonomy.get(rank, "unidentified")
        assigned_full_lineage = full_lineage(winner.taxonomy, rank)
        if not matches_lineage(winner, rank, lineage_match):
            continue
        read_weight = len(pair_hits)
        for hit in pair_hits:
            mate_full_lineage = full_lineage(hit.taxonomy, rank)
            selected.append(
                SelectedRead(
                    sample_id=hit.sample_id,
                    marker=hit.marker,
                    pair_id=pair_id,
                    mate=hit.mate,
                    strand=hit.strand,
                    target_id=hit.target_id,
                    identity=hit.identity,
                    aln_len=hit.aln_len,
                    mapq=hit.mapq,
                    mate_euk_group=hit.taxonomy.get("euk_group", NA_GROUP),
                    mate_lineage=hit.taxonomy.get(rank, "unidentified"),
                    mate_full_lineage=mate_full_lineage,
                    assigned_euk_group=winner.taxonomy.get("euk_group", NA_GROUP),
                    assigned_lineage=assigned_lineage,
                    assigned_full_lineage=assigned_full_lineage,
                    assigned_target_id=winner.target_id,
                    read_weight=read_weight,
                )
            )
    return selected


def select_sample_marker(
    index: int,
    sample_id: str,
    marker: str,
    pipeline_outdir: Path,
    taxonomy: Dict[str, Dict[str, str]],
    min_identity: Dict[str, float],
    min_aln_len: Dict[str, int],
    min_qcov: Dict[str, float],
    min_mapq: int,
    rank: str,
    lineage_match: Callable[[str], bool],
    mode: str,
) -> SelectResult:
    t0 = time.perf_counter()
    paf = pipeline_outdir / "03_align" / marker / f"{sample_id}.{marker}.pe.paf"
    hits = parse_paf_best_hits(
        paf,
        sample_id,
        marker,
        taxonomy,
        min_identity,
        min_aln_len,
        min_qcov,
        min_mapq,
    )
    selected = select_reads(hits, rank, lineage_match, mode)
    elapsed = time.perf_counter() - t0
    return SelectResult(
        index=index,
        sample_id=sample_id,
        marker=marker,
        paf_path=paf,
        best_hits=len(hits),
        selected_reads=len(selected),
        elapsed=elapsed,
        selected=selected,
    )


def select_lineage_reads_parallel(
    samples: Sequence[Dict[str, str]],
    markers: Sequence[str],
    pipeline_outdir: Path,
    taxonomy: Dict[str, Dict[str, str]],
    min_identity: Dict[str, float],
    min_aln_len: Dict[str, int],
    min_qcov: Dict[str, float],
    min_mapq: int,
    rank: str,
    lineage_match: Callable[[str], bool],
    mode: str,
    jobs: int,
    console: Optional[Console],
    quiet: bool,
) -> Tuple[List[SelectedRead], List[SelectResult]]:
    tasks: List[Tuple[int, str, str]] = []
    for sample in samples:
        sample_id = sample["sample_id"]
        for marker in markers:
            tasks.append((len(tasks), sample_id, marker))

    results: List[SelectResult] = []
    max_workers = max(1, min(jobs, len(tasks) or 1))
    progress_ctx = progress_bar(console, "parse_paf", len(tasks)) if console and not quiet and tasks else nullcontext(None)
    with progress_ctx as progress:
        task_id = progress.add_task("sample-marker PAF", total=len(tasks)) if progress else None
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futs = {
                executor.submit(
                    select_sample_marker,
                    index,
                    sample_id,
                    marker,
                    pipeline_outdir,
                    taxonomy,
                    min_identity,
                    min_aln_len,
                    min_qcov,
                    min_mapq,
                    rank,
                    lineage_match,
                    mode,
                ): (sample_id, marker)
                for index, sample_id, marker in tasks
            }
            for fut in as_completed(futs):
                result = fut.result()
                results.append(result)
                if progress and task_id is not None:
                    progress.update(
                        task_id,
                        advance=1,
                        description=f"{result.sample_id}:{result.marker}",
                    )

    selected: List[SelectedRead] = []
    for result in sorted(results, key=lambda x: x.index):
        selected.extend(result.selected)
    return selected, results


def iter_fastq(path: Path) -> Iterator[Tuple[str, str, str, str]]:
    with open_text_auto(path, "rt") as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            if not qual:
                raise ValueError(f"Truncated FASTQ record in {path}: {header.rstrip()}")
            yield header.rstrip("\r\n"), seq.rstrip("\r\n"), plus.rstrip("\r\n"), qual.rstrip("\r\n")


def revcomp(seq: str) -> str:
    table = str.maketrans("ACGTRYKMSWBDHVNacgtrykmswbdhvn", "TGCAYRMKSWVHDBNtgcayrmkswvhdbn")
    return seq.translate(table)[::-1]


def wrap_fasta(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def safe_name(value: str, max_len: int = 80) -> str:
    out = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    out = out.strip("._-") or "lineage"
    return out[:max_len]


def write_selected_table(path: Path, selected: Sequence[SelectedRead]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(
            [
                "sample_id",
                "marker",
                "pair_id",
                "mate",
                "strand",
                "target_id",
                "identity",
                "aln_len",
                "mapq",
                "mate_euk_group",
                "mate_lineage",
                "mate_full_lineage",
                "assigned_euk_group",
                "assigned_lineage",
                "assigned_full_lineage",
                "assigned_target_id",
                "read_weight",
            ]
        )
        for r in selected:
            writer.writerow(
                [
                    r.sample_id,
                    r.marker,
                    r.pair_id,
                    r.mate,
                    r.strand,
                    r.target_id,
                    f"{r.identity:.6f}",
                    r.aln_len,
                    r.mapq,
                    r.mate_euk_group,
                    r.mate_lineage,
                    r.mate_full_lineage,
                    r.assigned_euk_group,
                    r.assigned_lineage,
                    r.assigned_full_lineage,
                    r.assigned_target_id,
                    r.read_weight,
                ]
            )


def write_summary(path: Path, rows: Sequence[Dict[str, object]]) -> None:
    fields = ["sample_id", "marker", "selected_pairs", "selected_reads", "r1_reads", "r2_reads"]
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def read_summary_totals(path: Path) -> Dict[str, int]:
    totals = {"rows": 0, "selected_pairs": 0, "selected_reads": 0, "r1_reads": 0, "r2_reads": 0}
    if not path.exists():
        return totals
    with open(path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            totals["rows"] += 1
            for key in ["selected_pairs", "selected_reads", "r1_reads", "r2_reads"]:
                totals[key] += int(row.get(key) or 0)
    return totals


def read_metric_tsv(path: Path) -> Dict[str, str]:
    metrics: Dict[str, str] = {}
    if not path.exists():
        return metrics
    with open(path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            key = row.get("metric")
            if key:
                metrics[key] = row.get("value", "")
    return metrics


def extract_fastq_records(
    outdir: Path,
    pipeline_outdir: Path,
    selected: Sequence[SelectedRead],
    orient: bool,
    console: Optional[Console] = None,
    log_file: Optional[Path] = None,
    quiet: bool = False,
) -> List[Dict[str, object]]:
    selected_lookup: Dict[Tuple[str, str, str, str], SelectedRead] = {
        (r.sample_id, r.marker, r.mate, r.pair_id): r for r in selected
    }
    by_sample_marker: Dict[Tuple[str, str], List[SelectedRead]] = {}
    for r in selected:
        by_sample_marker.setdefault((r.sample_id, r.marker), []).append(r)

    r1_out = outdir / "reads.R1.fq.gz"
    r2_out = outdir / "reads.R2.fq.gz"
    segments_out = outdir / "merged_segments.fq.gz"
    fasta_out = outdir / "merged_segments.fa"

    summary: List[Dict[str, object]] = []
    written = set()
    total_units = sum(
        1
        for reads in by_sample_marker.values()
        for mate in ["R1", "R2"]
        if any(r.mate == mate for r in reads)
    )
    progress_ctx = progress_bar(console, "extract_fastq", total_units) if console and not quiet and total_units else nullcontext(None)
    with progress_ctx as progress:
        task_id = progress.add_task("mate FASTQ", total=total_units) if progress else None
        with open_gzip_text(r1_out, "wt") as r1_fh, open_gzip_text(r2_out, "wt") as r2_fh, open_gzip_text(
            segments_out, "wt"
        ) as segments_fh, open(fasta_out, "w", encoding="utf-8") as fa_fh:
            for (sample_id, marker), reads in sorted(by_sample_marker.items()):
                pair_ids = {r.pair_id for r in reads}
                mate_counts = {"R1": 0, "R2": 0}
                for mate in ["R1", "R2"]:
                    fastq = pipeline_outdir / "02_marker_reads" / marker / f"{sample_id}.{marker}.{mate}.fq.gz"
                    if not fastq.exists():
                        raise FileNotFoundError(f"Marker FASTQ not found: {fastq}")
                    wanted = {r.pair_id for r in reads if r.mate == mate}
                    if not wanted:
                        continue
                    mate_fh = r1_fh if mate == "R1" else r2_fh
                    for header, seq, plus, qual in iter_fastq(fastq):
                        pair_id = normalize_read_id(header)
                        key = (sample_id, marker, mate, pair_id)
                        selected_read = selected_lookup.get(key)
                        if pair_id not in wanted or selected_read is None or key in written:
                            continue
                        written.add(key)
                        new_id = f"{sample_id}|{marker}|{mate}|{pair_id}"
                        record = f"@{new_id}\n{seq}\n+\n{qual}\n"
                        mate_fh.write(record)
                        fa_seq = revcomp(seq) if orient and selected_read.strand == "-" else seq
                        fa_qual = qual[::-1] if orient and selected_read.strand == "-" else qual
                        segments_fh.write(f"@{new_id}\n{fa_seq}\n+\n{fa_qual}\n")
                        fa_fh.write(f">{new_id} target={selected_read.target_id} strand={selected_read.strand}\n")
                        fa_fh.write(wrap_fasta(fa_seq) + "\n")
                        mate_counts[mate] += 1
                    if progress and task_id is not None:
                        progress.update(task_id, advance=1, description=f"{sample_id}:{marker}:{mate}")
                summary.append(
                    {
                        "sample_id": sample_id,
                        "marker": marker,
                        "selected_pairs": len(pair_ids),
                        "selected_reads": mate_counts["R1"] + mate_counts["R2"],
                        "r1_reads": mate_counts["R1"],
                        "r2_reads": mate_counts["R2"],
                    }
                )

    missing = len(selected) - len(written)
    if missing:
        message = f"{missing} selected reads were not found in marker FASTQ files."
        if console and not quiet:
            console_log(console, log_file, "WARN", message, "yellow")
        else:
            print(f"WARNING: {message}", file=sys.stderr)
    return summary


def require_tool(name: str) -> str:
    path = shutil.which(name)
    if not path:
        raise RuntimeError(f"Missing required external tool in PATH: {name}")
    return path


def run_command(cmd: Sequence[str], log_file: Path, stdout_file: Optional[Path] = None) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    with open(log_file, "w", encoding="utf-8") as log:
        log.write(f"[{now()}]\tSTART\tcmd={shlex.join([str(x) for x in cmd])}\n")
        if stdout_file is not None:
            log.write(f"[{now()}]\tSTDOUT\tpath={stdout_file}\n")
        if stdout_file is None:
            proc = subprocess.run([str(x) for x in cmd], stdout=log, stderr=subprocess.STDOUT, text=True)
        else:
            with open(stdout_file, "w", encoding="utf-8") as out:
                proc = subprocess.run([str(x) for x in cmd], stdout=out, stderr=log, text=True)
        elapsed = time.perf_counter() - t0
        log.write(f"[{now()}]\tDONE\treturncode={proc.returncode}\telapsed={format_seconds(elapsed)}\n")
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed; log={log_file}")


def count_fastq_records(path: Path) -> int:
    count = 0
    for count, _record in enumerate(iter_fastq(path), start=1):
        pass
    return count


def count_lines(path: Path) -> int:
    if not path.exists():
        return 0
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        return sum(1 for _ in fh)


def gfa_to_fasta(gfa_path: Path, fasta_path: Path, min_contig_len: int) -> Dict[str, int]:
    contigs = 0
    total_bases = 0
    longest = 0
    with open(gfa_path, "r", encoding="utf-8", errors="replace") as gfa, open(fasta_path, "w", encoding="utf-8") as out:
        for line in gfa:
            if not line.startswith("S\t"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            name, seq = fields[1], fields[2]
            if seq == "*" or len(seq) < min_contig_len:
                continue
            contigs += 1
            total_bases += len(seq)
            longest = max(longest, len(seq))
            out.write(f">{name} len={len(seq)}\n{wrap_fasta(seq)}\n")
    return {"contigs": contigs, "total_bases": total_bases, "longest_contig": longest}


def run_miniasm_assembly(
    segments_fastq: Path,
    outdir: Path,
    threads: int,
    overlap_preset: str,
    minimap2_opts: Sequence[str],
    miniasm_opts: Sequence[str],
    min_contig_len: int,
    console: Optional[Console] = None,
    log_file: Optional[Path] = None,
    quiet: bool = False,
) -> Dict[str, int | str]:
    if not segments_fastq.exists():
        raise FileNotFoundError(f"merged segments FASTQ not found: {segments_fastq}")
    if threads < 1:
        raise ValueError("--assembly-threads must be positive.")

    minimap2 = require_tool("minimap2")
    miniasm = require_tool("miniasm")
    t0 = time.perf_counter()
    input_reads = count_fastq_records(segments_fastq)

    overlaps_paf = outdir / "miniasm.overlaps.paf"
    gfa_out = outdir / "miniasm.gfa"
    fasta_out = outdir / "miniasm_contigs.fa"
    summary_out = outdir / "assembly_summary.tsv"

    overlap_cmd = [
        minimap2,
        "-x",
        overlap_preset,
        "-t",
        str(threads),
        *minimap2_opts,
        str(segments_fastq),
        str(segments_fastq),
    ]
    assembly_cmd = [miniasm, *miniasm_opts, "-f", str(segments_fastq), str(overlaps_paf)]
    if console and not quiet:
        console_log(
            console,
            log_file,
            "ASSEMBLE",
            f"miniasm input={segments_fastq} reads={input_reads} threads={threads}",
            "cyan",
        )

    progress_ctx = progress_bar(console, "assemble", 3) if console and not quiet else nullcontext(None)
    with progress_ctx as progress:
        task_id = progress.add_task("minimap2/miniasm", total=3) if progress else None
        run_command(overlap_cmd, outdir / "minimap2.overlap.log", stdout_file=overlaps_paf)
        if progress and task_id is not None:
            progress.update(task_id, advance=1, description="minimap2 overlaps")

        run_command(assembly_cmd, outdir / "miniasm.log", stdout_file=gfa_out)
        if progress and task_id is not None:
            progress.update(task_id, advance=1, description="miniasm graph")

        contig_stats = gfa_to_fasta(gfa_out, fasta_out, min_contig_len)
        if progress and task_id is not None:
            progress.update(task_id, advance=1, description="export fasta")

    elapsed = time.perf_counter() - t0
    overlaps = count_lines(overlaps_paf)
    summary: Dict[str, int | str] = {
        "segments_fastq": str(segments_fastq),
        "input_reads": input_reads,
        "overlaps": overlaps,
        "contigs": contig_stats["contigs"],
        "total_bases": contig_stats["total_bases"],
        "longest_contig": contig_stats["longest_contig"],
        "overlap_preset": overlap_preset,
        "minimap2_opts": " ".join(minimap2_opts),
        "miniasm_opts": " ".join(miniasm_opts),
        "elapsed": format_seconds(elapsed),
        "overlap_rate": rate_text(overlaps, elapsed, "overlaps"),
    }
    with open(summary_out, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["metric", "value"])
        for key, value in summary.items():
            writer.writerow([key, value])
    return summary


def choose_taxonomy(args: argparse.Namespace, cfg: Dict[str, object]) -> Path:
    if args.taxonomy:
        return Path(args.taxonomy)
    value = cfg.get("taxonomy")
    if isinstance(value, str) and value:
        return Path(value)
    raise ValueError("Taxonomy path not found. Provide --taxonomy or use an outdir with run_config.tsv.")


def choose_markers(args: argparse.Namespace, cfg: Dict[str, object]) -> List[str]:
    markers = parse_list(args.markers)
    if markers:
        return markers
    cfg_markers = cfg.get("markers")
    if isinstance(cfg_markers, list) and cfg_markers:
        return [str(x) for x in cfg_markers]
    return MARKER_ORDER


def input_help_text() -> str:
    return dedent("""
    [bold green]Input modes[/bold green]

    1. Lineage extraction mode

       Required or assumed:
         --lineage LINEAGE
         --outdir OUTDIR  (default: metamarker_profile_out)

       OUTDIR must be a completed metamarker_profile output directory containing:

         sample_manifest.tsv
         run_config.tsv
         02_marker_reads/<MARKER>/<SAMPLE>.<MARKER>.R1.fq.gz
         02_marker_reads/<MARKER>/<SAMPLE>.<MARKER>.R2.fq.gz
         03_align/<MARKER>/<SAMPLE>.<MARKER>.pe.paf

       The taxonomy TSV is read from run_config.tsv unless --taxonomy is set.
       Required taxonomy columns:

         ref_id marker domain euk_group phylum class order family genus species

       euk_group is metadata for eukaryotic group separation, not a rank.

       PAF query IDs must end in /1 or /2. The script strips that mate suffix
       and pairs mates by normalized read ID.

    2. Assembly-only mode

       Required:
         --segments-fastq merged_segments.fq.gz

       The FASTQ must be gzip-compressed and contain standard four-line FASTQ
       records. This mode skips lineage extraction and runs minimap2/miniasm
       directly on the supplied merged_segments.fq.gz.

    [bold green]Lineage matching[/bold green]

      --rank controls which taxonomy column is matched.
      --lineage is matched exactly by default.
      It can be either the short value at --rank, such as Vibrio,
      or the complete lineage path through --rank, such as:
        d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Vibrionales;f__Vibrionaceae;g__Vibrio
      --match contains or --match regex can be used for broader matching.

    [bold green]Selection modes[/bold green]

      assigned  Matches the abundance table after R1/R2 arbitration.
      mate-hit  Keeps only mates whose own best hit matches the lineage.
    """).strip()


def default_help_text() -> str:
    return dedent(f"""
    [bold green]Default parameters[/bold green]

      --outdir             {DEFAULT_OUTDIR}
      --rank               {DEFAULT_RANK}
      --mode               {DEFAULT_MODE}
      --match              {DEFAULT_MATCH}
      --assembler          {DEFAULT_ASSEMBLER}
      --assembly-threads   {DEFAULT_ASSEMBLY_THREADS}
      --overlap-preset     {DEFAULT_OVERLAP_PRESET}
      --jobs               {DEFAULT_JOBS}
      --min-contig-len     0
      --resume             enabled by default
      --log width          terminal width × {DEFAULT_LOG_WIDTH_FRACTION:g}

    [bold green]Thresholds inherited from run_config.tsv[/bold green]

      min_identity:
        16S={DEFAULT_MIN_IDENTITY["16S"]} 18S={DEFAULT_MIN_IDENTITY["18S"]} ITS={DEFAULT_MIN_IDENTITY["ITS"]}
      min_aln_len:
        16S={DEFAULT_MIN_ALN_LEN["16S"]} 18S={DEFAULT_MIN_ALN_LEN["18S"]} ITS={DEFAULT_MIN_ALN_LEN["ITS"]}
      min_qcov:
        16S={DEFAULT_MIN_QCOV["16S"]} 18S={DEFAULT_MIN_QCOV["18S"]} ITS={DEFAULT_MIN_QCOV["ITS"]}
      min_mapq:
        {DEFAULT_MIN_MAPQ}

    [bold green]Assembly defaults[/bold green]

      minimap2 overlap command:
        minimap2 -x {DEFAULT_OVERLAP_PRESET} -t {DEFAULT_ASSEMBLY_THREADS} merged_segments.fq.gz merged_segments.fq.gz

      miniasm command:
        miniasm -f merged_segments.fq.gz miniasm.overlaps.paf
    """).strip()


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="extract_lineage_reads_assemble.py",
        formatter_class=RawDescriptionRichHelpFormatter,
        description=dedent("""
        [bold]Extract reads assigned to a target lineage and assemble marker segments.[/bold]

        The script first reproduces metamarker_profile read assignment from PAF,
        writes a flat set of selected FASTQ/TSV files, then assembles
        merged_segments.fq.gz using minimap2 all-vs-all overlaps and miniasm.
        """),
        epilog=dedent("""
        [bold green]Examples[/bold green]

          python3 scripts/extract_lineage_reads_assemble.py \\
            --outdir metamarker_profile_out \\
            --rank genus \\
            --lineage "Targetus" \\
            --markers 16S \\
            --assembler miniasm

          python3 scripts/extract_lineage_reads_assemble.py \\
            --segments-fastq metamarker_profile_out/lineage_extract/genus_Targetus/merged_segments.fq.gz \\
            --assembler miniasm \\
            --miniasm-opts "-m 40 -s 40"
        """),
        add_help=False,
    )

    req = p.add_argument_group("Required arguments")
    req.add_argument("--lineage", metavar="STR", help="Lineage value to extract at --rank. Required unless --segments-fastq is used.")
    req.add_argument("--segments-fastq", metavar="FILE", help="Existing merged_segments.fq.gz for assembly-only mode.")

    inp = p.add_argument_group("Optional arguments - input and matching")
    inp.add_argument("--outdir", default=DEFAULT_OUTDIR, metavar="DIR", help=f"metamarker_profile output directory. Default: {DEFAULT_OUTDIR}.")
    inp.add_argument("--taxonomy", metavar="FILE", help="Override taxonomy TSV. Default: path from OUTDIR/run_config.tsv.")
    inp.add_argument("--rank", default=DEFAULT_RANK, choices=RANKS, help=f"Taxonomic rank to match. Default: {DEFAULT_RANK}.")
    inp.add_argument("--markers", metavar="LIST", help="Comma-separated markers. Default: markers from run_config.tsv, otherwise 16S,18S,ITS.")
    inp.add_argument("--sample-id", metavar="LIST", help="Comma-separated sample IDs. Default: all samples in sample_manifest.tsv.")
    inp.add_argument(
        "--mode",
        choices=["assigned", "mate-hit"],
        default=DEFAULT_MODE,
        help="assigned matches abundance after pair arbitration; mate-hit keeps only matching mates. Default: assigned.",
    )
    inp.add_argument("--match", choices=["exact", "contains", "regex"], default=DEFAULT_MATCH, help=f"Lineage match mode. Default: {DEFAULT_MATCH}.")
    inp.add_argument("--ignore-case", action="store_true", help="Case-insensitive lineage matching. Default: false.")
    inp.add_argument("--output-dir", metavar="DIR", help="Output directory. Default: OUTDIR/lineage_extract/<rank>_<lineage>.")

    asm = p.add_argument_group("Optional arguments - assembly")
    asm.add_argument("--no-orient", action="store_true", help="Do not reverse-complement reads whose best PAF hit is on the minus strand. Default: false.")
    asm.add_argument(
        "--assembler",
        choices=["miniasm", "none"],
        default=DEFAULT_ASSEMBLER,
        help=f"Assembly backend. Default: {DEFAULT_ASSEMBLER}.",
    )
    asm.add_argument("--assembly-threads", type=int, default=DEFAULT_ASSEMBLY_THREADS, metavar="INT", help=f"Threads for minimap2 overlap search. Default: {DEFAULT_ASSEMBLY_THREADS}.")
    asm.add_argument("--overlap-preset", default=DEFAULT_OVERLAP_PRESET, metavar="STR", help=f"minimap2 all-vs-all preset. Default: {DEFAULT_OVERLAP_PRESET}.")
    asm.add_argument("--minimap2-opts", default="", metavar="STR", help='Extra minimap2 overlap options, quoted as one string. Default: "".')
    asm.add_argument("--miniasm-opts", default="", metavar="STR", help='Extra miniasm options, quoted as one string. Default: "".')
    asm.add_argument("--min-contig-len", type=int, default=0, metavar="INT", help="Minimum contig length exported from GFA to FASTA. Default: 0.")

    run = p.add_argument_group("Optional arguments - runtime, resume and logging")
    run.add_argument("--jobs", type=int, default=DEFAULT_JOBS, metavar="INT", help=f"Parallel sample-marker PAF parsing jobs. Default: {DEFAULT_JOBS}.")
    run.add_argument("--force", action="store_true", help="Recompute outputs even when resume checks pass. Default: false.")
    run.add_argument("--no-resume", action="store_true", help="Disable checkpoint/output checks. Default: false.")
    run.add_argument("--quiet", action="store_true", help="Suppress Rich run-plan and progress details. Default: false.")

    helpg = p.add_argument_group("Help")
    helpg.add_argument("--help_input", "--help-input", action="store_true", help="Show detailed input file requirements and exit.")
    helpg.add_argument("--help_default", "--help-default", action="store_true", help="Show detailed default parameter notes and exit.")
    helpg.add_argument("-h", "--help", action="help", help="Show this help message and exit.")
    return p


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    console = make_console()

    if args.help_input:
        console.print(Panel(input_help_text(), title="[bold cyan]Input File Requirements[/bold cyan]", border_style="cyan"))
        return 0
    if args.help_default:
        console.print(Panel(default_help_text(), title="[bold cyan]Default Parameters[/bold cyan]", border_style="cyan"))
        return 0

    validate_args(args)
    total_t0 = time.perf_counter()

    if args.segments_fastq:
        segments_fastq = Path(args.segments_fastq)
        output_dir = Path(args.output_dir) if args.output_dir else segments_fastq.parent
        output_dir.mkdir(parents=True, exist_ok=True)
        run_log = output_dir / "extract_lineage_reads_assemble.log"
        append_run_log(run_log, "START", f"mode=assembly_only\tsegments_fastq={segments_fastq}")
        if not args.quiet:
            print_run_plan(console, args, output_dir, run_log, "assembly-only", segments_fastq=segments_fastq)

        assembly_summary: Dict[str, int | str] = {}
        if args.assembler == "none":
            console_log(console, run_log, "INFO", "No assembly requested because --assembler none was set.", "yellow")
        elif resume_ready(run_log, "assembly", assembly_outputs(output_dir), args.force, args.no_resume):
            console_log(console, run_log, "RESUME", "assembly outputs already completed; skipping miniasm.", "yellow")
            assembly_summary = read_metric_tsv(output_dir / "assembly_summary.tsv")
        else:
            assembly_t0 = time.perf_counter()
            summary = run_miniasm_assembly(
                segments_fastq,
                output_dir,
                args.assembly_threads,
                args.overlap_preset,
                shlex.split(args.minimap2_opts),
                shlex.split(args.miniasm_opts),
                args.min_contig_len,
                console=console,
                log_file=run_log,
                quiet=args.quiet,
            )
            assembly_summary = summary
            mark_stage_done(
                run_log,
                "assembly",
                (
                    f"elapsed={format_seconds(time.perf_counter() - assembly_t0)}"
                    f"\tinput_reads={summary['input_reads']}\toverlaps={summary['overlaps']}\tcontigs={summary['contigs']}"
                ),
            )

        total_elapsed = time.perf_counter() - total_t0
        write_status(
            output_dir / "run_status.tsv",
            [
                ("status", "completed"),
                ("mode", "assembly_only"),
                ("segments_fastq", str(segments_fastq)),
                ("assembler", args.assembler),
                ("elapsed", format_seconds(total_elapsed)),
            ],
        )
        final_rows: List[Tuple[str, object]] = [
            ("mode", "assembly_only"),
            ("assembler", args.assembler),
            ("elapsed", format_seconds(total_elapsed)),
            ("output_dir", output_dir),
        ]
        if assembly_summary:
            final_rows.extend(
                [
                    ("input_reads", assembly_summary.get("input_reads", "")),
                    ("overlaps", assembly_summary.get("overlaps", "")),
                    ("contigs", assembly_summary.get("contigs", "")),
                    ("longest_contig", assembly_summary.get("longest_contig", "")),
                ]
            )
        print_final_summary(console, "Completed", final_rows)
        return 0

    if not args.lineage:
        raise ValueError("--lineage is required unless --segments-fastq is used.")

    pipeline_outdir = Path(args.outdir)
    cfg = parse_run_config(pipeline_outdir / "run_config.tsv")
    taxonomy_path = choose_taxonomy(args, cfg)
    taxonomy = load_taxonomy(taxonomy_path)
    markers = choose_markers(args, cfg)
    invalid_markers = [m for m in markers if m not in MARKER_ORDER]
    if invalid_markers:
        raise ValueError(f"Unsupported marker(s): {','.join(invalid_markers)}")

    sample_filter = set(parse_list(args.sample_id) or [])
    samples = load_manifest(pipeline_outdir / "sample_manifest.tsv")
    if sample_filter:
        samples = [s for s in samples if s["sample_id"] in sample_filter]
    if not samples:
        raise ValueError("No samples selected.")

    min_identity, min_aln_len, min_qcov, min_mapq = get_thresholds(cfg)
    lineage_match = build_matcher(args.lineage, args.match, args.ignore_case)

    output_dir = Path(args.output_dir) if args.output_dir else pipeline_outdir / "lineage_extract" / (
        f"{args.rank}_{safe_name(args.lineage)}"
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    run_log = output_dir / "extract_lineage_reads_assemble.log"
    append_run_log(
        run_log,
        "START",
        f"mode=lineage_extract\tlineage={args.lineage}\trank={args.rank}\tmarkers={','.join(markers)}",
    )
    if not args.quiet:
        print_run_plan(console, args, output_dir, run_log, "lineage-extract", markers=markers, sample_count=len(samples))

    if resume_ready(run_log, "lineage_extract", extraction_outputs(output_dir), args.force, args.no_resume):
        console_log(console, run_log, "RESUME", "lineage extraction outputs already completed; skipping PAF/FASTQ scan.", "yellow")
    else:
        select_t0 = time.perf_counter()
        selected, select_results = select_lineage_reads_parallel(
            samples,
            markers,
            pipeline_outdir,
            taxonomy,
            min_identity,
            min_aln_len,
            min_qcov,
            min_mapq,
            args.rank,
            lineage_match,
            args.mode,
            args.jobs,
            console,
            args.quiet,
        )
        select_elapsed = time.perf_counter() - select_t0
        console_log(
            console,
            run_log,
            "SELECT",
            (
                f"tasks={len(select_results)} best_hits={sum(x.best_hits for x in select_results)} "
                f"selected_reads={len(selected)} elapsed={format_seconds(select_elapsed)} "
                f"rate={rate_text(len(select_results), select_elapsed, 'tasks')}"
            ),
            "cyan",
        )
        write_selected_table(output_dir / "selected_reads.tsv", selected)
        if not selected:
            write_summary(output_dir / "summary.tsv", [])
            write_empty_extraction_outputs(output_dir)
            mark_stage_done(run_log, "lineage_extract", "selected_reads=0")
            write_status(
                output_dir / "run_status.tsv",
                [
                    ("status", "completed"),
                    ("mode", "lineage_extract"),
                    ("lineage", args.lineage),
                    ("rank", args.rank),
                    ("selected_reads", "0"),
                    ("elapsed", format_seconds(time.perf_counter() - total_t0)),
                ],
            )
            console_log(
                console,
                run_log,
                "WARN",
                f"No reads matched lineage={args.lineage!r} rank={args.rank}; wrote empty outputs to {output_dir}",
                "yellow",
            )
            print_final_summary(
                console,
                "Completed",
                [
                    ("lineage", args.lineage),
                    ("selected_reads", 0),
                    ("elapsed", format_seconds(time.perf_counter() - total_t0)),
                    ("output_dir", output_dir),
                ],
            )
            return 0

        extract_t0 = time.perf_counter()
        summary = extract_fastq_records(
            output_dir,
            pipeline_outdir,
            selected,
            orient=not args.no_orient,
            console=console,
            log_file=run_log,
            quiet=args.quiet,
        )
        write_summary(output_dir / "summary.tsv", summary)
        extract_elapsed = time.perf_counter() - extract_t0
        total_reads = sum(int(row["selected_reads"]) for row in summary)
        total_pairs = sum(int(row["selected_pairs"]) for row in summary)
        console_log(
            console,
            run_log,
            "EXTRACT",
            (
                f"selected_reads={total_reads} selected_pairs={total_pairs} "
                f"elapsed={format_seconds(extract_elapsed)} rate={rate_text(total_reads, extract_elapsed, 'reads')}"
            ),
            "cyan",
        )
        mark_stage_done(
            run_log,
            "lineage_extract",
            f"selected_reads={total_reads}\tselected_pairs={total_pairs}\telapsed={format_seconds(extract_elapsed)}",
        )

    segments_fastq = output_dir / "merged_segments.fq.gz"
    totals = read_summary_totals(output_dir / "summary.tsv")
    assembly_summary: Dict[str, int | str] = {}
    if args.assembler == "none":
        console_log(console, run_log, "INFO", "No assembly requested because --assembler none was set.", "yellow")
    elif totals["selected_reads"] == 0:
        console_log(console, run_log, "INFO", "No selected reads; skipping miniasm assembly.", "yellow")
    elif resume_ready(run_log, "assembly", assembly_outputs(output_dir), args.force, args.no_resume):
        console_log(console, run_log, "RESUME", "assembly outputs already completed; skipping miniasm.", "yellow")
        assembly_summary = read_metric_tsv(output_dir / "assembly_summary.tsv")
    else:
        assembly_t0 = time.perf_counter()
        assembly_summary = run_miniasm_assembly(
            segments_fastq,
            output_dir,
            args.assembly_threads,
            args.overlap_preset,
            shlex.split(args.minimap2_opts),
            shlex.split(args.miniasm_opts),
            args.min_contig_len,
            console=console,
            log_file=run_log,
            quiet=args.quiet,
        )
        mark_stage_done(
            run_log,
            "assembly",
            (
                f"elapsed={format_seconds(time.perf_counter() - assembly_t0)}"
                f"\tinput_reads={assembly_summary['input_reads']}"
                f"\toverlaps={assembly_summary['overlaps']}\tcontigs={assembly_summary['contigs']}"
            ),
        )

    total_elapsed = time.perf_counter() - total_t0
    write_status(
        output_dir / "run_status.tsv",
        [
            ("status", "completed"),
            ("mode", "lineage_extract"),
            ("lineage", args.lineage),
            ("rank", args.rank),
            ("markers", ",".join(markers)),
            ("samples", str(len(samples))),
            ("selected_pairs", str(totals["selected_pairs"])),
            ("selected_reads", str(totals["selected_reads"])),
            ("r1_reads", str(totals["r1_reads"])),
            ("r2_reads", str(totals["r2_reads"])),
            ("assembler", args.assembler),
            ("elapsed", format_seconds(total_elapsed)),
        ],
    )

    final_rows = [
        ("lineage", args.lineage),
        ("rank", args.rank),
        ("selected_pairs", totals["selected_pairs"]),
        ("selected_reads", totals["selected_reads"]),
        ("R1/R2", f"{totals['r1_reads']}/{totals['r2_reads']}"),
        ("assembler", args.assembler),
        ("elapsed", format_seconds(total_elapsed)),
        ("output_dir", output_dir),
    ]
    if assembly_summary:
        final_rows.extend(
            [
                ("overlaps", assembly_summary.get("overlaps", "")),
                ("contigs", assembly_summary.get("contigs", "")),
                ("longest_contig", assembly_summary.get("longest_contig", "")),
            ]
        )
    print_final_summary(console, "Completed", final_rows)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except BrokenPipeError:
        raise SystemExit(0)
    except Exception as exc:
        make_console().print(f"[bold red]ERROR[/bold red] {exc}")
        raise SystemExit(1)
