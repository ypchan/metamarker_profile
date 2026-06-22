#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
metamarker_profile


Workflow
--------
1. Count clean paired-end reads.
2. Extract marker candidate reads with BBDuk.
3. Align marker candidate R1/R2 reads together with one minimap2 call and write paired PAF.
4. Parse PAF with Python + Polars and calculate marker RPM/RPKM.

Design notes
------------
- External align/extract tools are still BBDuk and minimap2.
- Polars is the abundance engine.
- Abundance is parallelized at sample level, while extract/align are scheduled at sample-marker level.
- R1/R2 are streamed into one paired minimap2 process per sample-marker task.
  Query names are normalized to /1 and /2 mate suffixes so minimap2 keeps PE pairing.
- R1/R2 taxonomic conflicts are arbitrated at pair level:
  the better mate hit wins and both mates are assigned to that winner.

Dependencies
------------
Python:
  rich
  rich-argparse
  polars

External:
  bbduk.sh
  minimap2
  seqkit  optional; used by default for fast read counting when available

Example
-------
metamarker_profile \
  --input data_path.tsv \
  --outdir metamarker_profile_out \
  --markers 16S,ITS \
  --jobs 8 \
  --threads-per-sample 4
"""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
import os
import re
import shutil
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from dataclasses import dataclass, asdict
from datetime import datetime
from multiprocessing import get_context
from pathlib import Path
from textwrap import dedent
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from rich.console import Console
from rich.logging import RichHandler
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


PROGRAM = Path(sys.argv[0]).name

SILVA_REF_PREFIX = "SILVA_138.2_SSURef_NR99_tax_silva"
UNITE_REF_PREFIX = "UNITE_public_19.02.2025"

DEFAULT_OUTDIR = "metamarker_profile_out"
DEFAULT_STEPS = "all"
DEFAULT_MARKERS = "16S,ITS"
DEFAULT_RANK = "genus"

DEFAULT_JOBS = 4
DEFAULT_THREADS_PER_SAMPLE = 4
DEFAULT_SEQKIT_THREADS = "auto"
DEFAULT_SEQKIT_MAX_THREADS = 4
DEFAULT_BBDUK_THREADS = "auto"
DEFAULT_BBDUK_MAX_THREADS = 12
DEFAULT_MINIMAP2_THREADS = "auto"
DEFAULT_MINIMAP2_MAX_THREADS = 12
DEFAULT_ABUNDANCE_JOBS = "auto"
DEFAULT_POLARS_THREADS = 2
DEFAULT_BBMAP_MEM = "64G"

DEFAULT_K = {"16S": 31, "18S": 31, "ITS": 25}
DEFAULT_HDIST = {"16S": 0, "18S": 0, "ITS": 0}
DEFAULT_MKH = {"16S": 1, "18S": 1, "ITS": 1}
DEFAULT_MINK = {"16S": 0, "18S": 0, "ITS": 0}

DEFAULT_MINIMAP2_PRESET = "sr"
DEFAULT_MINIMAP2_N = 10

DEFAULT_MIN_IDENTITY = {"16S": 0.97, "18S": 0.97, "ITS": 0.95}
DEFAULT_MIN_ALN_LEN = {"16S": 80, "18S": 80, "ITS": 80}
DEFAULT_MIN_QCOV = {"16S": 0.60, "18S": 0.60, "ITS": 0.60}
DEFAULT_MIN_MAPQ = 0

RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
MARKER_ORDER = ["16S", "18S", "ITS"]

# Rich argparse style
RawDescriptionRichHelpFormatter.styles["argparse.prog"] = "bold magenta"
RawDescriptionRichHelpFormatter.styles["argparse.groups"] = "bold green"
RawDescriptionRichHelpFormatter.styles["argparse.args"] = "bold cyan"
RawDescriptionRichHelpFormatter.styles["argparse.metavar"] = "bold yellow"
RawDescriptionRichHelpFormatter.styles["argparse.help"] = "white"
RawDescriptionRichHelpFormatter.styles["argparse.text"] = "white"
RawDescriptionRichHelpFormatter.styles["argparse.syntax"] = "bold"


@dataclass
class Sample:
    row_no: int
    sample_id: str
    year: str
    month: str
    depth: str
    r1_path: str
    r2_path: str
    extra: Dict[str, str]


@dataclass
class Paths:
    outdir: str
    manifest: str
    clean_out: str
    marker_dir: str
    align_dir: str
    abund_dir: str
    tmp_dir: str
    clean_tmp_dir: str
    abundance_tmp_dir: str
    status_dir: str
    task_log_dir: str
    command_dir: str
    main_log: str


@dataclass
class Config:
    markers: List[str]
    rank: str
    steps: List[str]
    ref_dir: str
    ref_16s: str
    ref_18s: str
    ref_its: str
    index_16s: str
    index_18s: str
    index_its: str
    taxonomy: str

    jobs: int
    threads_per_sample: int
    total_thread_budget: int
    seqkit_threads: int
    bbduk_threads: int
    minimap2_threads: int
    abundance_jobs: int
    polars_threads: int
    bbmap_mem: str

    k: Dict[str, int]
    hdist: Dict[str, int]
    mkh: Dict[str, int]
    mink: Dict[str, int]
    minimap2_preset: str
    minimap2_n: int

    min_identity: Dict[str, float]
    min_aln_len: Dict[str, int]
    min_qcov: Dict[str, float]
    min_mapq: int

    force: bool
    resume: bool
    clean_tmp: bool
    dry_run: bool
    keep_task_logs: bool
    read_count_method: str


def now() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def terminal_width(default: int = 120) -> int:
    return shutil.get_terminal_size((default, 24)).columns


def make_console() -> Console:
    width = max(80, min(terminal_width(), 160))
    return Console(stderr=True, width=width, highlight=True, soft_wrap=True)


def format_seconds(seconds: float) -> str:
    seconds = max(0.0, seconds)
    if seconds < 60:
        return f"{seconds:.1f}s"
    minutes, sec = divmod(int(seconds), 60)
    hours, minutes = divmod(minutes, 60)
    if hours:
        return f"{hours}h{minutes:02d}m{sec:02d}s"
    return f"{minutes}m{sec:02d}s"


def setup_logger(console: Console, log_file: str, verbose: bool = False) -> logging.Logger:
    logger = logging.getLogger("metamarker_profile")
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    logger.handlers.clear()
    logger.propagate = False

    rich_handler = RichHandler(
        console=console,
        show_time=True,
        show_level=True,
        show_path=False,
        markup=True,
        rich_tracebacks=True,
    )
    rich_handler.setFormatter(logging.Formatter("%(message)s"))
    logger.addHandler(rich_handler)

    file_handler = logging.FileHandler(log_file, encoding="utf-8")
    file_handler.setFormatter(logging.Formatter("[%(asctime)s] [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S"))
    file_handler.setLevel(logging.DEBUG)
    logger.addHandler(file_handler)

    return logger


def default_ref_base() -> Path:
    return Path.cwd()


def read_default_ref_dir(base_dir: Path) -> str:
    if os.environ.get("METAMARKER_PROFILE_REF_DIR"):
        return os.environ["METAMARKER_PROFILE_REF_DIR"]

    config_file = Path(
        os.environ.get(
            "METAMARKER_PROFILE_REF_CONFIG",
            Path.home() / ".config" / "metamarker_profile" / "ref_dir",
        )
    )
    if config_file.exists() and config_file.stat().st_size > 0:
        for raw in config_file.read_text(encoding="utf-8").splitlines():
            line = raw.split("#", 1)[0].strip()
            if not line:
                continue
            p = Path(line)
            if p.is_absolute():
                return str(p)
            return str(base_dir / p)

    return str(base_dir / "refs")


def normalize_markers(raw: str) -> List[str]:
    tokens = [x.strip().upper() for x in raw.replace(";", ",").split(",") if x.strip()]
    out: List[str] = []
    for t in tokens:
        if t not in MARKER_ORDER:
            raise ValueError(f"Unsupported marker: {t}. Allowed: 16S,18S,ITS")
        if t not in out:
            out.append(t)
    if not out:
        raise ValueError("--markers is empty.")
    return out


def normalize_steps(raw: str) -> List[str]:
    raw = raw.strip().lower().replace(";", ",").replace(" ", "")
    if raw == "all":
        return ["reads_count", "extract", "align", "abundance"]

    aliases = {
        "read_count": "reads_count",
        "read-count": "reads_count",
        "reads-count": "reads_count",
        "count": "reads_count",
        "counts": "reads_count",
        "clean": "reads_count",
    }
    out: List[str] = []
    for t in raw.split(","):
        if not t:
            continue
        step = aliases.get(t, t)
        if step not in {"reads_count", "extract", "align", "abundance"}:
            raise ValueError(f"Unsupported step: {t}. Allowed: all or reads_count,extract,align,abundance")
        if step not in out:
            out.append(step)
    if not out:
        raise ValueError("--steps is empty.")
    return out


def is_pos_int(x: Any) -> bool:
    try:
        return int(x) > 0
    except Exception:
        return False


def resolve_tool_threads(label: str, requested: str, cap: int, threads_per_sample: int, logger: logging.Logger) -> int:
    if requested in {"auto", "0", 0}:
        req = threads_per_sample
    else:
        req = int(requested)
        if req < 1:
            raise ValueError(f"{label} threads must be positive, auto, or 0.")
    if req > cap:
        logger.warning(f"[yellow]{label} threads capped at {cap}; requested {req}.[/yellow]")
    return max(1, min(req, cap))


def calc_step_jobs(sample_count: int, total_budget: int, tool_threads: int) -> int:
    return max(1, min(sample_count, total_budget // max(1, tool_threads)))


def marker_ref(cfg: Config, marker: str) -> str:
    return {"16S": cfg.ref_16s, "18S": cfg.ref_18s, "ITS": cfg.ref_its}[marker]


def marker_index(cfg: Config, marker: str) -> str:
    return {"16S": cfg.index_16s, "18S": cfg.index_18s, "ITS": cfg.index_its}[marker]


def ensure_dir(path: str | Path) -> None:
    Path(path).mkdir(parents=True, exist_ok=True)


def open_text(path: str | Path, mode: str = "rt"):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode, encoding=None if "b" in mode else "utf-8")


def shlex_join(cmd: Sequence[str]) -> str:
    import shlex
    return " ".join(shlex.quote(str(x)) for x in cmd)


def safe_fraction(numerator: int | float, denominator: int | float) -> float:
    denominator = float(denominator)
    if denominator <= 0:
        return 0.0
    return float(numerator) / denominator


def fmt_float(value: int | float) -> str:
    return f"{float(value):.10g}"


def write_command(paths: Paths, stage: str, sample_id: str, cmd: Sequence[str]) -> None:
    d = Path(paths.command_dir) / stage
    ensure_dir(d)
    f = d / f"{sample_id}.commands.txt"
    with open(f, "a", encoding="utf-8") as out:
        out.write(f"# [{now()}] stage={stage} sample_id={sample_id}\n")
        out.write(f"cd {Path.cwd()}\n")
        out.write(shlex_join([str(x) for x in cmd]) + "\n\n")


def run_cmd(
    cmd: Sequence[str],
    log_file: str,
    stage: str,
    sample_id: str,
    paths: Paths,
    stdout_file: Optional[str] = None,
) -> None:
    write_command(paths, stage, sample_id, cmd)
    ensure_dir(Path(log_file).parent)
    t0 = time.time()
    with open(log_file, "a", encoding="utf-8") as log:
        log.write(f"[{now()}] [START] stage={stage} sample={sample_id}\n")
        log.write(f"[{now()}] [CMD] {shlex_join(cmd)}\n")
        if stdout_file:
            ensure_dir(Path(stdout_file).parent)
            with open(stdout_file, "w", encoding="utf-8") as out:
                proc = subprocess.run([str(x) for x in cmd], stdout=out, stderr=log, text=True)
        else:
            proc = subprocess.run([str(x) for x in cmd], stdout=log, stderr=subprocess.STDOUT, text=True)
        elapsed = format_seconds(time.time() - t0)
        level = "DONE" if proc.returncode == 0 else "ERROR"
        log.write(f"[{now()}] [{level}] stage={stage} sample={sample_id} elapsed={elapsed} returncode={proc.returncode}\n")
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed [{stage}] sample={sample_id}; log={log_file}")


def command_exists(tool: str) -> bool:
    return shutil.which(tool) is not None


def require_file(label: str, path: str, opt: str) -> None:
    if not Path(path).is_file() or Path(path).stat().st_size == 0:
        raise FileNotFoundError(f"Missing {label}: {path}. Provide {opt} or --ref-dir.")


def init_paths(outdir: str) -> Paths:
    out = Path(outdir)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    paths = Paths(
        outdir=str(out),
        manifest=str(out / "sample_manifest.tsv"),
        clean_out=str(out / "reads_stat.tsv"),
        marker_dir=str(out / "02_marker_reads"),
        align_dir=str(out / "03_align"),
        abund_dir=str(out),
        tmp_dir=str(out / ".tmp"),
        clean_tmp_dir=str(out / ".tmp" / "01_reads_count"),
        abundance_tmp_dir=str(out / ".tmp" / "abundance"),
        status_dir=str(out / ".checkpoints"),
        task_log_dir=str(out / ".tmp" / "task_logs"),
        command_dir=str(out / "commands"),
        main_log=str(out / f"metamarker_profile.{timestamp}.log"),
    )

    for d in [
        paths.outdir,
        paths.marker_dir,
        paths.align_dir,
        paths.tmp_dir,
        paths.clean_tmp_dir,
        paths.abundance_tmp_dir,
        paths.task_log_dir,
        paths.command_dir,
        Path(paths.command_dir) / "reads_count",
        Path(paths.command_dir) / "extract",
        Path(paths.command_dir) / "align",
        Path(paths.command_dir) / "abundance",
        Path(paths.status_dir) / "reads_count",
        Path(paths.status_dir) / "extract",
        Path(paths.status_dir) / "align",
        Path(paths.status_dir) / "abundance",
    ]:
        ensure_dir(d)

    for marker in MARKER_ORDER:
        ensure_dir(Path(paths.marker_dir) / marker)
        ensure_dir(Path(paths.align_dir) / marker)
    ensure_dir(Path(paths.marker_dir) / "stats")
    ensure_dir(Path(paths.marker_dir) / "tmp")
    ensure_dir(Path(paths.align_dir) / "tmp")
    Path(paths.main_log).touch()
    return paths


def read_samples_from_input(path: str) -> Tuple[List[Sample], List[str]]:
    samples: List[Sample] = []
    with open(path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("Input TSV has no header.")
        required = {"sample_id", "r1_path", "r2_path"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"Input TSV missing required columns: {','.join(sorted(missing))}")
        standard = {"sample_id", "year", "month", "depth", "r1_path", "r2_path"}
        extra_cols = [c for c in reader.fieldnames if c not in standard]
        for row_no, row in enumerate(reader, start=1):
            sid = (row.get("sample_id") or "").strip()
            if not sid:
                raise ValueError(f"Empty sample_id at input row {row_no}")
            if "/" in sid:
                raise ValueError(f"sample_id contains slash at input row {row_no}: {sid}")
            r1 = (row.get("r1_path") or "").strip()
            r2 = (row.get("r2_path") or "").strip()
            if not r1 or not r2:
                raise ValueError(f"Empty r1_path/r2_path at input row {row_no}")
            extra = {c: (row.get(c) or "NA").strip() or "NA" for c in extra_cols}
            samples.append(
                Sample(
                    row_no=row_no,
                    sample_id=sid,
                    year=(row.get("year") or "NA").strip() or "NA",
                    month=(row.get("month") or "NA").strip() or "NA",
                    depth=(row.get("depth") or "NA").strip() or "NA",
                    r1_path=r1,
                    r2_path=r2,
                    extra=extra,
                )
            )
    if not samples:
        raise ValueError("No sample rows found.")
    return samples, extra_cols


def build_single_sample(args: argparse.Namespace) -> Tuple[List[Sample], List[str]]:
    sid = args.sample_id
    if "/" in sid:
        raise ValueError(f"sample_id contains slash: {sid}")
    return [
        Sample(
            row_no=1,
            sample_id=sid,
            year="NA",
            month="NA",
            depth="NA",
            r1_path=args.r1,
            r2_path=args.r2,
            extra={},
        )
    ], []


def write_manifest(samples: List[Sample], extra_cols: List[str], paths: Paths) -> None:
    with open(paths.manifest, "w", encoding="utf-8", newline="") as out:
        fields = ["sample_id", "year", "month", "depth", "r1_path", "r2_path"] + extra_cols
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for s in samples:
            row = {
                "sample_id": s.sample_id,
                "year": s.year,
                "month": s.month,
                "depth": s.depth,
                "r1_path": s.r1_path,
                "r2_path": s.r2_path,
            }
            row.update({c: s.extra.get(c, "NA") for c in extra_cols})
            writer.writerow(row)


def check_inputs(samples: List[Sample]) -> None:
    for s in samples:
        if not Path(s.r1_path).is_file() or Path(s.r1_path).stat().st_size == 0:
            raise FileNotFoundError(f"Missing R1 for {s.sample_id}: {s.r1_path}")
        if not Path(s.r2_path).is_file() or Path(s.r2_path).stat().st_size == 0:
            raise FileNotFoundError(f"Missing R2 for {s.sample_id}: {s.r2_path}")


def read_config(args: argparse.Namespace, paths: Paths, logger: logging.Logger) -> Config:
    base_dir = default_ref_base()
    ref_dir = args.ref_dir or read_default_ref_dir(base_dir)

    ref_16s = args.ref_16s or str(Path(ref_dir) / f"{SILVA_REF_PREFIX}.dna.arc_bac.shortid.fasta")
    ref_18s = args.ref_18s or str(Path(ref_dir) / f"{SILVA_REF_PREFIX}.dna.euk.shortid.fasta")
    ref_its = args.ref_its or str(Path(ref_dir) / f"{UNITE_REF_PREFIX}.shortid.fasta")

    index_16s = args.index_16s or f"{ref_16s}.mmi"
    index_18s = args.index_18s or f"{ref_18s}.mmi"
    index_its = args.index_its or f"{ref_its}.mmi"
    taxonomy = args.taxonomy or str(Path(ref_dir) / "ref_taxonomy.tsv")

    markers = normalize_markers(args.markers)
    steps = normalize_steps(args.steps)

    total_budget = args.jobs * args.threads_per_sample
    seqkit_threads = resolve_tool_threads("seqkit", str(args.seqkit_threads), args.seqkit_max_threads, args.threads_per_sample, logger)
    bbduk_threads = resolve_tool_threads("BBDuk", str(args.bbduk_threads), args.bbduk_max_threads, args.threads_per_sample, logger)
    minimap2_threads = resolve_tool_threads("minimap2", str(args.minimap2_threads), args.minimap2_max_threads, args.threads_per_sample, logger)

    # abundance is sample-level parallel + Polars internal threads.
    if args.abundance_jobs == "auto":
        abundance_jobs = max(1, total_budget // max(1, args.polars_threads))
    else:
        abundance_jobs = int(args.abundance_jobs)

    return Config(
        markers=markers,
        rank=args.rank,
        steps=steps,
        ref_dir=ref_dir,
        ref_16s=ref_16s,
        ref_18s=ref_18s,
        ref_its=ref_its,
        index_16s=index_16s,
        index_18s=index_18s,
        index_its=index_its,
        taxonomy=taxonomy,
        jobs=args.jobs,
        threads_per_sample=args.threads_per_sample,
        total_thread_budget=total_budget,
        seqkit_threads=seqkit_threads,
        bbduk_threads=bbduk_threads,
        minimap2_threads=minimap2_threads,
        abundance_jobs=abundance_jobs,
        polars_threads=args.polars_threads,
        bbmap_mem=args.bbmap_mem,
        k={"16S": args.k_16s, "18S": args.k_18s, "ITS": args.k_its},
        hdist={"16S": args.hdist_16s, "18S": args.hdist_18s, "ITS": args.hdist_its},
        mkh={"16S": args.mkh_16s, "18S": args.mkh_18s, "ITS": args.mkh_its},
        mink={"16S": args.mink_16s, "18S": args.mink_18s, "ITS": args.mink_its},
        minimap2_preset=args.minimap2_preset,
        minimap2_n=args.minimap2_N,
        min_identity={"16S": args.min_identity_16s, "18S": args.min_identity_18s, "ITS": args.min_identity_its},
        min_aln_len={"16S": args.min_aln_len_16s, "18S": args.min_aln_len_18s, "ITS": args.min_aln_len_its},
        min_qcov={"16S": args.min_qcov_16s, "18S": args.min_qcov_18s, "ITS": args.min_qcov_its},
        min_mapq=args.min_mapq,
        force=args.force,
        resume=not args.no_resume,
        clean_tmp=args.clean_tmp,
        dry_run=args.dry_run,
        keep_task_logs=args.keep_task_logs,
        read_count_method=args.read_count_method,
    )


def check_dependencies(cfg: Config, logger: logging.Logger) -> None:
    missing = []
    if "extract" in cfg.steps and not command_exists("bbduk.sh"):
        missing.append("bbduk.sh")
    if "align" in cfg.steps and not command_exists("minimap2"):
        missing.append("minimap2")
    if "reads_count" in cfg.steps and cfg.read_count_method == "seqkit" and not command_exists("seqkit"):
        missing.append("seqkit")

    if "abundance" in cfg.steps:
        # Strong dependency; no try wrapper. Missing polars should fail clearly.
        import polars  # noqa: F401

    if missing:
        raise RuntimeError(f"Missing external tools in PATH: {', '.join(missing)}")

    if "extract" in cfg.steps:
        if "16S" in cfg.markers:
            require_file("16S reference FASTA", cfg.ref_16s, "--ref-16s")
        if "18S" in cfg.markers:
            require_file("18S reference FASTA", cfg.ref_18s, "--ref-18s")
        if "ITS" in cfg.markers:
            require_file("ITS reference FASTA", cfg.ref_its, "--ref-its")

    if "align" in cfg.steps:
        if "16S" in cfg.markers:
            require_file("16S minimap2 index", cfg.index_16s, "--index-16s")
        if "18S" in cfg.markers:
            require_file("18S minimap2 index", cfg.index_18s, "--index-18s")
        if "ITS" in cfg.markers:
            require_file("ITS minimap2 index", cfg.index_its, "--index-its")

    if "abundance" in cfg.steps:
        require_file("taxonomy table", cfg.taxonomy, "--taxonomy")

    logger.info("[green]Dependency and reference check passed[/green]")


def print_run_plan(console: Console, cfg: Config, paths: Paths, samples: List[Sample], extra_cols: List[str]) -> None:
    table = Table.grid(padding=(0, 1))
    table.add_column(style="bold cyan", justify="right")
    table.add_column(style="white")
    table.add_row("Samples", f"{len(samples):,}")
    table.add_row("Markers", ",".join(cfg.markers))
    table.add_row("Steps", ",".join(cfg.steps))
    table.add_row("Rank", cfg.rank)
    table.add_row("Jobs", str(cfg.jobs))
    table.add_row("Threads/sample", str(cfg.threads_per_sample))
    table.add_row("CPU budget", str(cfg.total_thread_budget))
    table.add_row("seqkit threads", str(cfg.seqkit_threads))
    table.add_row("BBDuk threads", str(cfg.bbduk_threads))
    table.add_row("minimap2 threads", str(cfg.minimap2_threads))
    table.add_row("minimap2 mode", "paired R1/R2 single process")
    table.add_row("abundance jobs", str(cfg.abundance_jobs))
    table.add_row("Polars threads/task", str(cfg.polars_threads))
    table.add_row("Outdir", paths.outdir)
    table.add_row("Manifest", paths.manifest)
    table.add_row("Taxonomy", cfg.taxonomy)
    table.add_row("Extra metadata", ",".join(extra_cols) if extra_cols else "none")
    console.print(Panel(table, title="[bold green]Run plan[/bold green]", border_style="green", expand=False))


def count_fastq_records_python(path: str) -> int:
    opener = gzip.open if path.endswith(".gz") else open
    n = 0
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for n, _ in enumerate(fh, start=1):
            pass
    if n % 4 != 0:
        raise ValueError(f"FASTQ line count is not divisible by 4: {path} lines={n}")
    return n // 4


def count_reads_one(sample: Sample, cfg: Config, paths: Paths) -> Dict[str, Any]:
    out = Path(paths.clean_tmp_dir) / f"{sample.row_no}.{sample.sample_id}.reads_stat.tsv"
    done = Path(paths.status_dir) / "reads_count" / f"{sample.sample_id}.done"
    log_file = Path(paths.task_log_dir) / f"reads_count.{sample.sample_id}.log"

    if cfg.resume and not cfg.force and out.exists() and out.stat().st_size > 0 and done.exists():
        return {"sample_id": sample.sample_id, "skipped": True}

    ensure_dir(out.parent)
    ensure_dir(done.parent)
    method = cfg.read_count_method
    if method == "auto":
        method = "seqkit" if command_exists("seqkit") else "python"

    if method == "seqkit":
        tmp = Path(paths.clean_tmp_dir) / f"{sample.row_no}.{sample.sample_id}.seqkit.stats.tmp"
        cmd = ["seqkit", "stats", "-j", str(cfg.seqkit_threads), "-T", sample.r1_path]
        run_cmd(cmd, str(log_file), "reads_count", sample.sample_id, paths, stdout_file=str(tmp))
        with open(tmp, "r", encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            row = next(reader)
            read_pairs = int(row.get("num_seqs", row.get("num_seqs ", "0")).replace(",", ""))
        tmp.unlink(missing_ok=True)
    else:
        with open(log_file, "w", encoding="utf-8") as log:
            log.write(f"[{now()}] [INFO] Native Python FASTQ count: {sample.r1_path}\n")
        read_pairs = count_fastq_records_python(sample.r1_path)

    clean_total = read_pairs * 2

    with open(out.with_suffix(out.suffix + ".tmp"), "w", encoding="utf-8") as tmp_out:
        tmp_out.write(
            "\t".join([
                sample.sample_id,
                sample.year,
                sample.month,
                sample.depth,
                sample.r1_path,
                sample.r2_path,
                str(read_pairs),
                str(clean_total),
            ]) + "\n"
        )
    out.with_suffix(out.suffix + ".tmp").replace(out)
    done.write_text(f"sample_id={sample.sample_id}\tread_pairs={read_pairs}\tclean_reads_total={clean_total}\n", encoding="utf-8")

    if not cfg.keep_task_logs:
        Path(log_file).unlink(missing_ok=True)

    return {"sample_id": sample.sample_id, "read_pairs": read_pairs, "clean_reads_total": clean_total}


def run_reads_count(samples: List[Sample], cfg: Config, paths: Paths, logger: logging.Logger) -> None:
    if cfg.resume and not cfg.force and Path(paths.clean_out).exists():
        with open(paths.clean_out, "r", encoding="utf-8") as fh:
            n = sum(1 for _ in fh) - 1
        if n == len(samples):
            logger.info(f"[green]Reads-stat table exists; skipping:[/green] {paths.clean_out}")
            return

    logger.info("[bold cyan]Step 01 reads_count: clean FASTQ read totals[/bold cyan]")
    jobs = calc_step_jobs(len(samples), cfg.total_thread_budget, cfg.seqkit_threads)
    logger.info(f"reads_count parallelism: jobs={jobs}, seqkit_threads={cfg.seqkit_threads}, cpu_budget={cfg.total_thread_budget}")

    if cfg.force:
        for p in Path(paths.clean_tmp_dir).glob("*.reads_stat.tsv"):
            p.unlink(missing_ok=True)

    progress_console = Console(stderr=True, width=terminal_width())
    with Progress(
        SpinnerColumn(),
        TextColumn("[bold cyan]reads_count"),
        BarColumn(bar_width=None),
        TaskProgressColumn(),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        console=progress_console,
    ) as progress:
        task = progress.add_task("reads_count", total=len(samples))
        with ThreadPoolExecutor(max_workers=jobs) as ex:
            futs = {ex.submit(count_reads_one, s, cfg, paths): s.sample_id for s in samples}
            for fut in as_completed(futs):
                sample_id = futs[fut]
                try:
                    fut.result()
                except Exception:
                    logger.exception(f"reads_count failed: sample={sample_id}")
                    raise
                progress.update(task, advance=1)

    with open(paths.clean_out + ".tmp", "w", encoding="utf-8") as out:
        out.write("sample_id\tyear\tmonth\tdepth\tr1_path\tr2_path\tread_pairs\tclean_reads_total\n")
        for s in samples:
            f = Path(paths.clean_tmp_dir) / f"{s.row_no}.{s.sample_id}.reads_stat.tsv"
            if not f.exists():
                raise FileNotFoundError(f"Missing reads_count temporary result: {f}")
            out.write(f.read_text(encoding="utf-8"))
    Path(paths.clean_out + ".tmp").replace(paths.clean_out)
    logger.info(f"[green]Written reads stats:[/green] {paths.clean_out}")


def extract_one_marker(sample: Sample, marker: str, cfg: Config, paths: Paths) -> None:
    out_r1 = Path(paths.marker_dir) / marker / f"{sample.sample_id}.{marker}.R1.fq.gz"
    out_r2 = Path(paths.marker_dir) / marker / f"{sample.sample_id}.{marker}.R2.fq.gz"
    stat = Path(paths.marker_dir) / "stats" / f"{sample.sample_id}.{marker}.bbduk.stats.txt"
    done = Path(paths.status_dir) / "extract" / f"{sample.sample_id}.{marker}.done"
    log_file = Path(paths.task_log_dir) / f"extract.{sample.sample_id}.{marker}.bbduk.log"
    tag = f"{sample.sample_id}.{marker}.{os.getpid()}"
    tmp_r1 = Path(paths.marker_dir) / "tmp" / f"{tag}.R1.fq.gz.tmp"
    tmp_r2 = Path(paths.marker_dir) / "tmp" / f"{tag}.R2.fq.gz.tmp"
    tmp_stat = Path(paths.marker_dir) / "tmp" / f"{tag}.stats.tmp"

    if cfg.resume and not cfg.force and out_r1.exists() and out_r2.exists() and stat.exists() and done.exists():
        return

    for p in [tmp_r1, tmp_r2, tmp_stat]:
        p.unlink(missing_ok=True)

    cmd = [
        "bbduk.sh",
        f"-Xmx{cfg.bbmap_mem}",
        f"in1={sample.r1_path}",
        f"in2={sample.r2_path}",
        f"outm1={tmp_r1}",
        f"outm2={tmp_r2}",
        f"ref={marker_ref(cfg, marker)}",
        f"k={cfg.k[marker]}",
        f"hdist={cfg.hdist[marker]}",
        f"mkh={cfg.mkh[marker]}",
        f"t={cfg.bbduk_threads}",
        "ordered=t",
        f"stats={tmp_stat}",
    ]
    if cfg.mink[marker] > 0:
        cmd.append(f"mink={cfg.mink[marker]}")

    run_cmd(cmd, str(log_file), "extract", sample.sample_id, paths)
    tmp_r1.replace(out_r1)
    tmp_r2.replace(out_r2)
    tmp_stat.replace(stat)
    done.write_text(f"sample_id={sample.sample_id}\tmarker={marker}\tout_r1={out_r1}\tout_r2={out_r2}\n", encoding="utf-8")

    if not cfg.keep_task_logs:
        log_file.unlink(missing_ok=True)


def run_extract(samples: List[Sample], cfg: Config, paths: Paths, logger: logging.Logger) -> None:
    logger.info("[bold cyan]Step 02 extract: marker candidate reads with BBDuk[/bold cyan]")
    tasks = [(s, marker) for s in samples for marker in cfg.markers]
    jobs = calc_step_jobs(len(tasks), cfg.total_thread_budget, cfg.bbduk_threads)
    logger.info(
        f"extract parallelism: tasks={len(tasks)}, jobs={jobs}, "
        f"bbduk_threads={cfg.bbduk_threads}, cpu_budget={cfg.total_thread_budget}"
    )

    progress_console = Console(stderr=True, width=terminal_width())
    with Progress(
        SpinnerColumn(),
        TextColumn("[bold cyan]extract"),
        BarColumn(bar_width=None),
        TaskProgressColumn(),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        console=progress_console,
    ) as progress:
        task = progress.add_task("extract", total=len(tasks))
        with ThreadPoolExecutor(max_workers=jobs) as ex:
            futs = {ex.submit(extract_one_marker, s, marker, cfg, paths): (s.sample_id, marker) for s, marker in tasks}
            for fut in as_completed(futs):
                sample_id, marker = futs[fut]
                try:
                    fut.result()
                except Exception:
                    logger.exception(f"extract failed: sample={sample_id} marker={marker}")
                    raise
                progress.update(task, advance=1)

    logger.info("[green]Step extract finished[/green]")


def write_prefixed_fastq(src: Path, dst: Path, mate: str) -> None:
    line_no = 0
    mate_no = "1" if mate == "R1" else "2"
    opener = gzip.open if str(src).endswith(".gz") else open
    with open(dst, "w", encoding="utf-8", buffering=1024 * 1024) as out, opener(
        src, "rt", encoding="utf-8", errors="replace"
    ) as inp:
        for line_no, line in enumerate(inp, start=1):
            if line_no % 4 == 1:
                if not line.startswith("@"):
                    raise ValueError(f"Invalid FASTQ header in {src} at line {line_no}")
                header = line[1:].rstrip("\r\n")
                name, sep, comment = header.partition(" ")
                name = re.sub(r"^R[12]::", "", name)
                name = re.sub(r"/[12]$", "", name)
                labeled_name = f"{name}/{mate_no}"
                out.write(f"@{labeled_name}{sep}{comment}\n")
            else:
                out.write(line)
    if line_no % 4 != 0:
        raise ValueError(f"FASTQ line count is not divisible by 4: {src} lines={line_no}")


def stream_prefixed_fastq(src: Path, fifo: Path, mate: str, errors: List[BaseException]) -> None:
    try:
        write_prefixed_fastq(src, fifo, mate)
    except BrokenPipeError as exc:
        errors.append(exc)
    except BaseException as exc:
        errors.append(exc)


def run_paired_minimap2(
    cmd: Sequence[str],
    r1: Path,
    r2: Path,
    out_paf: Path,
    log_file: Path,
    stage: str,
    sample_id: str,
    marker: str,
    paths: Paths,
) -> None:
    ensure_dir(log_file.parent)
    ensure_dir(out_paf.parent)

    fifo_dir = Path(paths.align_dir) / "tmp"
    ensure_dir(fifo_dir)
    tag = f"{sample_id}.{marker}.{os.getpid()}.{time.time_ns()}"
    query_r1 = fifo_dir / f"{tag}.R1.fq"
    query_r2 = fifo_dir / f"{tag}.R2.fq"
    tmp_paf = out_paf.with_suffix(out_paf.suffix + ".tmp")

    for p in [query_r1, query_r2, tmp_paf]:
        p.unlink(missing_ok=True)

    use_fifo = True
    fifo_error = ""
    try:
        os.mkfifo(query_r1)
        os.mkfifo(query_r2)
    except OSError as exc:
        use_fifo = False
        fifo_error = str(exc)
        for p in [query_r1, query_r2]:
            p.unlink(missing_ok=True)
        try:
            write_prefixed_fastq(r1, query_r1, "R1")
            write_prefixed_fastq(r2, query_r2, "R2")
        except Exception:
            for p in [query_r1, query_r2, tmp_paf]:
                p.unlink(missing_ok=True)
            raise

    runtime_cmd = [str(x) for x in cmd]
    runtime_cmd[-2:] = [str(query_r1), str(query_r2)]
    write_command(paths, stage, sample_id, runtime_cmd)
    writer_errors: List[BaseException] = []
    t0 = time.time()

    try:
        with open(log_file, "w", encoding="utf-8") as log, open(tmp_paf, "w", encoding="utf-8") as paf:
            log.write(f"[{now()}] [START] stage={stage} sample={sample_id} marker={marker}\n")
            log.write(f"[{now()}] [INFO] R1 source={r1} streamed with query suffix /1\n")
            log.write(f"[{now()}] [INFO] R2 source={r2} streamed with query suffix /2\n")
            if not use_fifo:
                log.write(f"[{now()}] [WARN] FIFO unavailable; using temporary prefixed FASTQ files: {fifo_error}\n")
            log.write(f"[{now()}] [CMD] {shlex_join(runtime_cmd)}\n")
            if use_fifo:
                proc = subprocess.Popen(runtime_cmd, stdout=paf, stderr=log, text=True)
                writers = [
                    threading.Thread(target=stream_prefixed_fastq, args=(r1, query_r1, "R1", writer_errors), daemon=True),
                    threading.Thread(target=stream_prefixed_fastq, args=(r2, query_r2, "R2", writer_errors), daemon=True),
                ]
                for th in writers:
                    th.start()
                returncode = proc.wait()
                for th in writers:
                    th.join()
            else:
                proc = subprocess.run(runtime_cmd, stdout=paf, stderr=log, text=True)
                returncode = proc.returncode
            elapsed = format_seconds(time.time() - t0)
            level = "DONE" if returncode == 0 else "ERROR"
            log.write(f"[{now()}] [{level}] stage={stage} sample={sample_id} marker={marker} elapsed={elapsed} returncode={returncode}\n")

        if returncode != 0:
            raise RuntimeError(f"minimap2 failed: sample={sample_id} marker={marker}; log={log_file}")

        real_errors = [e for e in writer_errors if not isinstance(e, BrokenPipeError)]
        if real_errors:
            msg = "; ".join(str(e) for e in real_errors[:3])
            raise RuntimeError(f"FASTQ streaming failed: sample={sample_id} marker={marker}; {msg}; log={log_file}")

        tmp_paf.replace(out_paf)
    finally:
        for p in [query_r1, query_r2, tmp_paf]:
            p.unlink(missing_ok=True)


def align_one_marker(sample: Sample, marker: str, cfg: Config, paths: Paths) -> None:
    index = marker_index(cfg, marker)
    done = Path(paths.status_dir) / "align" / f"{sample.sample_id}.{marker}.done"
    out_paf = Path(paths.align_dir) / marker / f"{sample.sample_id}.{marker}.pe.paf"
    in_r1 = Path(paths.marker_dir) / marker / f"{sample.sample_id}.{marker}.R1.fq.gz"
    in_r2 = Path(paths.marker_dir) / marker / f"{sample.sample_id}.{marker}.R2.fq.gz"
    log_file = Path(paths.task_log_dir) / f"align.{sample.sample_id}.{marker}.minimap2.log"

    if cfg.resume and not cfg.force and out_paf.exists() and done.exists():
        return

    if not in_r1.exists():
        raise FileNotFoundError(f"Missing marker R1 reads: {in_r1}")
    if not in_r2.exists():
        raise FileNotFoundError(f"Missing marker R2 reads: {in_r2}")

    cmd = [
        "minimap2",
        "-x", cfg.minimap2_preset,
        "-N", str(cfg.minimap2_n),
        "-t", str(cfg.minimap2_threads),
        index,
        "R1_FIFO",
        "R2_FIFO",
    ]
    run_paired_minimap2(cmd, in_r1, in_r2, out_paf, log_file, "align", sample.sample_id, marker, paths)

    done.write_text(f"sample_id={sample.sample_id}\tmarker={marker}\tout_paf={out_paf}\n", encoding="utf-8")
    if not cfg.keep_task_logs:
        log_file.unlink(missing_ok=True)


def run_align(samples: List[Sample], cfg: Config, paths: Paths, logger: logging.Logger) -> None:
    logger.info("[bold cyan]Step 03 align: paired marker reads with minimap2[/bold cyan]")
    tasks = [(s, marker) for s in samples for marker in cfg.markers]
    jobs = calc_step_jobs(len(tasks), cfg.total_thread_budget, cfg.minimap2_threads)
    logger.info(
        f"align parallelism: tasks={len(tasks)}, jobs={jobs}, "
        f"minimap2_threads={cfg.minimap2_threads}, cpu_budget={cfg.total_thread_budget}, mode=paired"
    )

    progress_console = Console(stderr=True, width=terminal_width())
    with Progress(
        SpinnerColumn(),
        TextColumn("[bold cyan]align"),
        BarColumn(bar_width=None),
        TaskProgressColumn(),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        console=progress_console,
    ) as progress:
        task = progress.add_task("align", total=len(tasks))
        with ThreadPoolExecutor(max_workers=jobs) as ex:
            futs = {ex.submit(align_one_marker, s, marker, cfg, paths): (s.sample_id, marker) for s, marker in tasks}
            for fut in as_completed(futs):
                sample_id, marker = futs[fut]
                try:
                    fut.result()
                except Exception:
                    logger.exception(f"align failed: sample={sample_id} marker={marker}")
                    raise
                progress.update(task, advance=1)

    logger.info("[green]Step align finished[/green]")


# ---------- Polars abundance worker globals ----------
_TAX = None
_RANKS_TO_OUTPUT: List[str] = []


def _abundance_worker_init(taxonomy_path: str, ranks_to_output: List[str], polars_threads: int) -> None:
    os.environ["POLARS_MAX_THREADS"] = str(polars_threads)
    import polars as pl

    global _TAX, _RANKS_TO_OUTPUT
    _RANKS_TO_OUTPUT = ranks_to_output

    cols = ["ref_id", "marker"] + RANKS
    tax = pl.read_csv(
        taxonomy_path,
        separator="\t",
        has_header=True,
        infer_schema_length=1000,
        ignore_errors=True,
    )
    missing = [c for c in cols if c not in tax.columns]
    if missing:
        raise ValueError(f"Taxonomy table missing columns: {','.join(missing)}")
    tax = tax.select(cols).with_columns([pl.col(c).cast(pl.Utf8).str.strip_chars().alias(c) for c in cols])
    tax = tax.with_columns(
        [
            pl.when(pl.col("domain").is_null() | (pl.col("domain") == ""))
            .then(pl.lit("NA"))
            .otherwise(pl.col("domain"))
            .alias("domain"),
            *[
                pl.when(pl.col(rank).is_null() | (pl.col(rank) == ""))
                .then(pl.lit("unidentified"))
                .otherwise(pl.col(rank))
                .alias(rank)
                for rank in RANKS
                if rank != "domain"
            ],
        ]
    )
    _TAX = tax


def strip_pair_id_expr(colname: str):
    import polars as pl
    return (
        pl.col(colname)
        .cast(pl.Utf8)
        .str.replace(r"/[12]$", "")
    )


def mate_from_query_expr(colname: str):
    import polars as pl
    q = pl.col(colname).cast(pl.Utf8)
    return (
        pl.when(q.str.ends_with("/1"))
        .then(pl.lit("R1"))
        .when(q.str.ends_with("/2"))
        .then(pl.lit("R2"))
        .otherwise(pl.lit("unknown"))
    )


def read_paf_best_hits(paf_path: str, marker: str, cfg_dict: Dict[str, Any]) -> Tuple[Any, Dict[str, int]]:
    import polars as pl

    if not Path(paf_path).exists() or Path(paf_path).stat().st_size == 0:
        empty = pl.DataFrame(
            schema={
                "pair_id": pl.Utf8,
                "mate": pl.Utf8,
                "target_id": pl.Utf8,
                "target_len": pl.Int64,
                "matches": pl.Int64,
                "aln_len": pl.Int64,
                "mapq": pl.Int64,
                "identity": pl.Float64,
                "qcov": pl.Float64,
            }
        )
        return empty, {"total_alignments": 0, "passed_alignments": 0}

    names = ["query_id", "query_len", "query_start", "query_end", "target_id", "target_len", "matches", "aln_len", "mapq"]
    df = pl.read_csv(
        paf_path,
        separator="\t",
        has_header=False,
        columns=[0, 1, 2, 3, 5, 6, 9, 10, 11],
        new_columns=names,
        infer_schema_length=1000,
        ignore_errors=True,
        truncate_ragged_lines=True,
        schema_overrides={
            "query_id": pl.Utf8,
            "query_len": pl.Int64,
            "query_start": pl.Int64,
            "query_end": pl.Int64,
            "target_id": pl.Utf8,
            "target_len": pl.Int64,
            "matches": pl.Int64,
            "aln_len": pl.Int64,
            "mapq": pl.Int64,
        },
    )

    total_alignments = df.height
    if total_alignments == 0:
        return df, {"total_alignments": 0, "passed_alignments": 0}

    tax_marker = _TAX.filter(pl.col("marker") == marker)

    df = (
        df.with_columns([
            mate_from_query_expr("query_id").alias("mate"),
            strip_pair_id_expr("query_id").alias("pair_id"),
            (pl.col("matches") / pl.col("aln_len")).alias("identity"),
            ((pl.col("query_end") - pl.col("query_start")) / pl.col("query_len")).alias("qcov"),
        ])
        .join(tax_marker, left_on="target_id", right_on="ref_id", how="inner")
        .filter(
            (pl.col("identity") >= float(cfg_dict["min_identity"][marker])) &
            (pl.col("aln_len") >= int(cfg_dict["min_aln_len"][marker])) &
            (pl.col("qcov") >= float(cfg_dict["min_qcov"][marker])) &
            (pl.col("mapq") >= int(cfg_dict["min_mapq"]))
        )
    )

    unknown_mates = df.filter(pl.col("mate") == "unknown").height
    if unknown_mates:
        raise ValueError(
            f"PAF query IDs lack /1,/2 mate labels: {paf_path}. "
            "Rerun the align step with the current paired minimap2 workflow."
        )

    passed_alignments = df.height
    if passed_alignments == 0:
        return df, {"total_alignments": total_alignments, "passed_alignments": 0}

    best = (
        df.sort(
            ["pair_id", "mate", "identity", "aln_len", "matches", "mapq", "qcov", "target_len", "target_id"],
            descending=[False, False, True, True, True, True, True, True, False],
        )
        .unique(subset=["pair_id", "mate"], keep="first", maintain_order=True)
    )
    return best, {"total_alignments": total_alignments, "passed_alignments": passed_alignments}


def assign_pairs(hits):
    import polars as pl

    if hits.height == 0:
        return hits.with_columns(pl.lit(0).alias("read_weight")), {
            "reads_with_passed_hits": 0,
            "read_pairs_with_passed_hits": 0,
            "assigned_reads": 0,
            "discordant_pairs": 0,
            "reassigned_discordant_reads": 0,
            "tie_broken_reads": 0,
            "ambiguous_reads": 0,
        }

    mate_counts = hits.group_by("pair_id").agg([
        pl.len().alias("mate_count"),
        pl.col("target_id").n_unique().alias("target_n_unique"),
    ])

    winner = (
        hits.sort(
            ["pair_id", "identity", "aln_len", "matches", "mapq", "qcov", "target_len", "mate", "target_id"],
            descending=[False, True, True, True, True, True, True, False, False],
        )
        .unique(subset=["pair_id"], keep="first", maintain_order=True)
        .join(mate_counts, on="pair_id", how="left")
        .with_columns(pl.col("mate_count").cast(pl.Int64).alias("read_weight"))
    )

    discordant_pairs = mate_counts.filter((pl.col("mate_count") == 2) & (pl.col("target_n_unique") > 1)).height
    assigned_reads = int(winner.get_column("read_weight").sum()) if winner.height else 0

    stats = {
        "reads_with_passed_hits": hits.select(["pair_id", "mate"]).unique().height,
        "read_pairs_with_passed_hits": hits.select("pair_id").unique().height,
        "assigned_reads": assigned_reads,
        "discordant_pairs": discordant_pairs,
        "reassigned_discordant_reads": discordant_pairs * 2,
        "tie_broken_reads": 0,
        "ambiguous_reads": 0,
    }
    return winner, stats


def write_abundance_rows(
    assigned,
    sample: Dict[str, Any],
    marker: str,
    clean_total: int,
    ranks_to_output: List[str],
    out_files: Dict[str, str],
    extra_cols: List[str],
) -> None:
    import polars as pl

    if assigned.height == 0:
        return

    for rank in ranks_to_output:
        if clean_total > 0:
            rpm_expr = (pl.col("taxon_marker_reads") / float(clean_total) * 1_000_000.0).alias("marker_rpm")
            rpkm_expr = (
                pl.col("taxon_marker_reads")
                * 1_000_000_000.0
                / float(clean_total)
                / (pl.col("target_len_sum") / pl.col("taxon_marker_reads"))
            ).alias("marker_rpkm")
        else:
            rpm_expr = pl.lit(0.0).alias("marker_rpm")
            rpkm_expr = pl.lit(0.0).alias("marker_rpkm")

        df = (
            assigned
            .with_columns([
                pl.col(rank).fill_null("unidentified").alias("lineage"),
                (pl.col("target_len") * pl.col("read_weight")).alias("weighted_tlen"),
            ])
            .group_by(["domain", "lineage"])
            .agg([
                pl.col("read_weight").sum().alias("taxon_marker_reads"),
                pl.col("weighted_tlen").sum().alias("target_len_sum"),
            ])
            .with_columns([
                (pl.col("target_len_sum") / pl.col("taxon_marker_reads")).alias("mean_target_length_bp"),
                rpm_expr,
                rpkm_expr,
            ])
            .select(["domain", "lineage", "taxon_marker_reads", "mean_target_length_bp", "marker_rpm", "marker_rpkm"])
            .sort(["domain", "taxon_marker_reads", "lineage"], descending=[False, True, False])
        )

        with open(out_files[rank], "a", encoding="utf-8") as out:
            for row in df.iter_rows(named=True):
                base = [
                    sample["sample_id"],
                    marker,
                    row["domain"] or "NA",
                    rank,
                    row["lineage"] or "unidentified",
                    str(int(row["taxon_marker_reads"])),
                    str(clean_total),
                    fmt_float(row["marker_rpm"]),
                    fmt_float(row["marker_rpkm"]),
                    f"{float(row['mean_target_length_bp']):.2f}",
                ]
                extras = [sample.get(c, "NA") for c in extra_cols]
                out.write("\t".join(base + extras) + "\n")


def abundance_process_one_sample(payload: Dict[str, Any]) -> Dict[str, Any]:
    sample = payload["sample"]
    cfg_dict = payload["cfg"]
    paths = payload["paths"]
    clean_total = int(payload["clean_total"])
    extra_cols = payload["extra_cols"]
    ranks_to_output = payload["ranks_to_output"]

    sample_id = sample["sample_id"]
    tmp_dir = Path(paths["abundance_tmp_dir"])
    ensure_dir(tmp_dir)

    out_files = {r: str(tmp_dir / f"{sample['row_no']:06d}.{sample_id}.{r}.long.part.tsv") for r in ranks_to_output}
    stat_files = {r: str(tmp_dir / f"{sample['row_no']:06d}.{sample_id}.{r}.stats.part.tsv") for r in ranks_to_output}
    qc_file = str(tmp_dir / f"{sample['row_no']:06d}.{sample_id}.assignment_qc.part.tsv")
    for p in list(out_files.values()) + list(stat_files.values()) + [qc_file]:
        Path(p).unlink(missing_ok=True)
        Path(p).touch()

    for marker in cfg_dict["markers"]:
        paf = Path(paths["align_dir"]) / marker / f"{sample_id}.{marker}.pe.paf"

        hits, hit_stats = read_paf_best_hits(str(paf), marker, cfg_dict)
        assigned, astats = assign_pairs(hits)

        total_alignments = hit_stats["total_alignments"]
        passed_alignments = hit_stats["passed_alignments"]
        assigned_reads = astats["assigned_reads"]
        reads_with_passed_hits = astats["reads_with_passed_hits"]
        read_pairs_with_passed_hits = astats["read_pairs_with_passed_hits"]
        discordant_pairs = astats["discordant_pairs"]

        write_abundance_rows(assigned, sample, marker, clean_total, ranks_to_output, out_files, extra_cols)

        with open(qc_file, "a", encoding="utf-8") as out:
            base = [
                sample_id,
                sample.get("year", "NA"),
                sample.get("month", "NA"),
                sample.get("depth", "NA"),
                marker,
                str(clean_total),
                str(total_alignments),
                str(passed_alignments),
                fmt_float(safe_fraction(passed_alignments, total_alignments)),
                str(reads_with_passed_hits),
                str(read_pairs_with_passed_hits),
                str(assigned_reads),
                fmt_float(safe_fraction(assigned_reads, clean_total)),
                fmt_float(safe_fraction(assigned_reads, reads_with_passed_hits)),
                str(discordant_pairs),
                fmt_float(safe_fraction(discordant_pairs, read_pairs_with_passed_hits)),
                str(astats["reassigned_discordant_reads"]),
                str(astats["tie_broken_reads"]),
                str(astats["ambiguous_reads"]),
            ]
            extras = [sample.get(c, "NA") for c in extra_cols]
            out.write("\t".join(base + extras) + "\n")

        for rank in ranks_to_output:
            with open(stat_files[rank], "a", encoding="utf-8") as out:
                base = [
                    sample_id,
                    sample.get("year", "NA"),
                    sample.get("month", "NA"),
                    sample.get("depth", "NA"),
                    marker,
                    str(total_alignments),
                    str(passed_alignments),
                    str(astats["reads_with_passed_hits"]),
                    str(astats["read_pairs_with_passed_hits"]),
                    str(astats["assigned_reads"]),
                    str(astats["discordant_pairs"]),
                    str(astats["reassigned_discordant_reads"]),
                    str(astats["tie_broken_reads"]),
                    str(astats["ambiguous_reads"]),
                ]
                extras = [sample.get(c, "NA") for c in extra_cols]
                out.write("\t".join(base + extras) + "\n")

    return {"sample_id": sample_id, "long_files": out_files, "stat_files": stat_files, "qc_file": qc_file}


def read_clean_counts(clean_out: str) -> Dict[str, int]:
    out: Dict[str, int] = {}
    with open(clean_out, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            out[row["sample_id"]] = int(row["clean_reads_total"])
    return out


def samples_to_payload(samples: List[Sample], clean_counts: Dict[str, int], cfg: Config, paths: Paths, extra_cols: List[str], ranks: List[str]) -> List[Dict[str, Any]]:
    cfg_dict = asdict(cfg)
    paths_dict = asdict(paths)
    payloads: List[Dict[str, Any]] = []
    for s in samples:
        if s.sample_id not in clean_counts:
            raise ValueError(f"Sample missing in reads-count table: {s.sample_id}")
        row = {
            "row_no": s.row_no,
            "sample_id": s.sample_id,
            "year": s.year,
            "month": s.month,
            "depth": s.depth,
            **s.extra,
        }
        payloads.append({
            "sample": row,
            "clean_total": clean_counts[s.sample_id],
            "cfg": cfg_dict,
            "paths": paths_dict,
            "extra_cols": extra_cols,
            "ranks_to_output": ranks,
        })
    return payloads


def merge_abundance_outputs(samples: List[Sample], ranks: List[str], paths: Paths, extra_cols: List[str], logger: logging.Logger) -> None:
    long_header = [
        "sample_id", "marker", "domain", "rank", "lineage",
        "taxon_marker_reads", "clean_reads_total", "marker_rpm", "marker_rpkm", "mean_target_length_bp",
    ] + extra_cols
    stats_header = [
        "sample_id", "year", "month", "depth", "marker",
        "total_alignments", "passed_alignments", "reads_with_passed_hits", "read_pairs_with_passed_hits",
        "assigned_reads", "discordant_pairs", "reassigned_discordant_reads", "tie_broken_reads", "ambiguous_reads",
    ] + extra_cols
    qc_header = [
        "sample_id", "year", "month", "depth", "marker", "clean_reads_total",
        "total_alignments", "passed_alignments", "passed_alignment_rate",
        "reads_with_passed_hits", "read_pairs_with_passed_hits",
        "assigned_reads", "assigned_read_rate", "assigned_read_fraction_of_passed_reads",
        "discordant_pairs", "discordant_pair_rate",
        "reassigned_discordant_reads", "tie_broken_reads", "ambiguous_reads",
    ] + extra_cols

    qc_out = Path(paths.abund_dir) / "all.marker_rpm.assignment_qc.tsv"
    with open(str(qc_out) + ".tmp", "w", encoding="utf-8") as out:
        out.write("\t".join(qc_header) + "\n")
        for s in samples:
            part = Path(paths.abundance_tmp_dir) / f"{s.row_no:06d}.{s.sample_id}.assignment_qc.part.tsv"
            if part.exists():
                out.write(part.read_text(encoding="utf-8"))
    Path(str(qc_out) + ".tmp").replace(qc_out)
    logger.info(f"[green]Written:[/green] {qc_out}")

    for rank in ranks:
        long_out = Path(paths.abund_dir) / f"all.marker_rpm.{rank}.long.tsv"
        stats_out = Path(paths.abund_dir) / f"all.marker_rpm.{rank}.assignment_stats.tsv"
        with open(str(long_out) + ".tmp", "w", encoding="utf-8") as out:
            out.write("\t".join(long_header) + "\n")
            for s in samples:
                part = Path(paths.abundance_tmp_dir) / f"{s.row_no:06d}.{s.sample_id}.{rank}.long.part.tsv"
                if part.exists():
                    out.write(part.read_text(encoding="utf-8"))
        Path(str(long_out) + ".tmp").replace(long_out)

        with open(str(stats_out) + ".tmp", "w", encoding="utf-8") as out:
            out.write("\t".join(stats_header) + "\n")
            for s in samples:
                part = Path(paths.abundance_tmp_dir) / f"{s.row_no:06d}.{s.sample_id}.{rank}.stats.part.tsv"
                if part.exists():
                    out.write(part.read_text(encoding="utf-8"))
        Path(str(stats_out) + ".tmp").replace(stats_out)

        logger.info(f"[green]Written:[/green] {long_out}")
        logger.info(f"[green]Written:[/green] {stats_out}")


def run_abundance(samples: List[Sample], cfg: Config, paths: Paths, extra_cols: List[str], logger: logging.Logger) -> None:
    logger.info("[bold cyan]Step 04 abundance: Polars PAF parser and RPM/RPKM calculation[/bold cyan]")
    if not Path(paths.clean_out).exists():
        raise FileNotFoundError(f"Reads-stat table not found: {paths.clean_out}")

    ranks = RANKS if cfg.rank == "all" else [cfg.rank]
    expected = Path(paths.abund_dir) / f"all.marker_rpm.{ranks[-1]}.long.tsv"
    if cfg.resume and not cfg.force and expected.exists() and expected.stat().st_size > 0:
        logger.info(f"[green]Abundance output exists; skipping:[/green] {expected}")
        return

    os.environ["POLARS_MAX_THREADS"] = str(cfg.polars_threads)
    clean_counts = read_clean_counts(paths.clean_out)
    jobs = max(1, min(len(samples), cfg.abundance_jobs))
    payloads = samples_to_payload(samples, clean_counts, cfg, paths, extra_cols, ranks)

    logger.info(f"abundance parallelism: jobs={jobs}, polars_threads_per_job={cfg.polars_threads}, cpu_budget≈{jobs * cfg.polars_threads}")

    progress_console = Console(stderr=True, width=terminal_width())
    with Progress(
        SpinnerColumn(),
        TextColumn("[bold cyan]abundance"),
        BarColumn(bar_width=None),
        TaskProgressColumn(),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        console=progress_console,
    ) as progress:
        task = progress.add_task("abundance", total=len(payloads))
        ctx = get_context("spawn")
        with ProcessPoolExecutor(
            max_workers=jobs,
            mp_context=ctx,
            initializer=_abundance_worker_init,
            initargs=(cfg.taxonomy, ranks, cfg.polars_threads),
        ) as ex:
            futs = {ex.submit(abundance_process_one_sample, p): p["sample"]["sample_id"] for p in payloads}
            for fut in as_completed(futs):
                sample_id = futs[fut]
                try:
                    fut.result()
                except Exception:
                    logger.exception(f"abundance failed: sample={sample_id}")
                    raise
                progress.update(task, advance=1)

    merge_abundance_outputs(samples, ranks, paths, extra_cols, logger)
    done = Path(paths.status_dir) / "abundance" / "abundance.done"
    done.write_text(f"rank={cfg.rank}\tmarkers={','.join(cfg.markers)}\n", encoding="utf-8")


def write_run_config(cfg: Config, paths: Paths) -> None:
    out = Path(paths.outdir) / "run_config.tsv"
    with open(out, "w", encoding="utf-8") as fh:
        for k, v in asdict(cfg).items():
            fh.write(f"{k}\t{v}\n")


def build_parser() -> argparse.ArgumentParser:
    base_dir = default_ref_base()
    default_ref_dir = read_default_ref_dir(base_dir)

    epilog = dedent(f"""
    [bold green]Default parameters[/bold green]

      --outdir              {DEFAULT_OUTDIR}
      --steps               {DEFAULT_STEPS}
      --markers             {DEFAULT_MARKERS}
      --rank                {DEFAULT_RANK}
      --jobs                {DEFAULT_JOBS}
      --threads-per-sample  {DEFAULT_THREADS_PER_SAMPLE}
      --bbmap-mem           {DEFAULT_BBMAP_MEM}
      --minimap2-preset     {DEFAULT_MINIMAP2_PRESET}
      --minimap2-N          {DEFAULT_MINIMAP2_N}
      --abundance-jobs      auto
      --polars-threads      {DEFAULT_POLARS_THREADS}

    [bold green]BBDuk defaults[/bold green]

      16S: k=31, hdist=0, mkh=1, mink=0
      18S: k=31, hdist=0, mkh=1, mink=0
      ITS: k=25, hdist=0, mkh=1, mink=0

    [bold green]Abundance filters[/bold green]

      16S: identity>=0.97, aln_len>=80, qcov>=0.60
      18S: identity>=0.97, aln_len>=80, qcov>=0.60
      ITS: identity>=0.95, aln_len>=80, qcov>=0.60
      mapq>=0

    [bold green]R1/R2 conflict rule[/bold green]

      R1 and R2 are aligned together in one minimap2 process per marker.
      PAF query IDs are normalized to /1 and /2 during streaming.
      R1 and R2 are counted at read level but arbitrated at pair level.
      If R1/R2 disagree, the better mate hit wins:
      identity, aln_len, matches, MAPQ, qcov, target_len, mate, target_id.

    [bold green]Default reference paths[/bold green]

      ref-dir  {default_ref_dir}
      16S      {default_ref_dir}/{SILVA_REF_PREFIX}.dna.arc_bac.shortid.fasta
      18S      {default_ref_dir}/{SILVA_REF_PREFIX}.dna.euk.shortid.fasta
      ITS      {default_ref_dir}/{UNITE_REF_PREFIX}.shortid.fasta
      taxonomy {default_ref_dir}/ref_taxonomy.tsv
    """)

    p = argparse.ArgumentParser(
        prog="metamarker_profile",
        formatter_class=RawDescriptionRichHelpFormatter,
        description=dedent("""
        [bold]Integrated marker-read counting pipeline for paired-end clean reads.[/bold]

        Workflow components:
          seqkit/Python read counting
          BBDuk candidate-read extraction
          paired R1/R2 minimap2 PAF alignment
          Polars abundance parser
        """),
        epilog=epilog,
        add_help=False,
    )

    req = p.add_argument_group("Required arguments")
    req.add_argument("-i", "--input", metavar="FILE", help="Multi-sample TSV input with header.")
    req.add_argument("--sample-id", metavar="ID", help="Single-sample ID.")
    req.add_argument("--r1", metavar="FILE", help="Single-sample R1 FASTQ(.gz).")
    req.add_argument("--r2", metavar="FILE", help="Single-sample R2 FASTQ(.gz).")

    ref = p.add_argument_group("Reference data")
    ref.add_argument("--ref-dir", metavar="DIR", help=f"Reference directory. Default: {default_ref_dir}")
    ref.add_argument("--ref-16s", metavar="FILE", help="Override 16S reference FASTA for BBDuk.")
    ref.add_argument("--ref-18s", metavar="FILE", help="Override 18S reference FASTA for BBDuk.")
    ref.add_argument("--ref-its", metavar="FILE", help="Override ITS reference FASTA for BBDuk.")
    ref.add_argument("--index-16s", metavar="FILE", help="Override minimap2 index for 16S.")
    ref.add_argument("--index-18s", metavar="FILE", help="Override minimap2 index for 18S.")
    ref.add_argument("--index-its", metavar="FILE", help="Override minimap2 index for ITS.")
    ref.add_argument("--taxonomy", metavar="FILE", help="Override ref_id/marker/domain/ranks taxonomy TSV.")

    flow = p.add_argument_group("Optional arguments - workflow control")
    flow.add_argument("-o", "--outdir", default=DEFAULT_OUTDIR, metavar="DIR", help=f"Output directory. Default: {DEFAULT_OUTDIR}.")
    flow.add_argument("--steps", default=DEFAULT_STEPS, metavar="LIST", help="all or comma list: reads_count,extract,align,abundance. Default: all.")
    flow.add_argument("--markers", default=DEFAULT_MARKERS, metavar="LIST", help="Marker list: 16S,18S,ITS. Default: 16S,ITS.")
    flow.add_argument("--rank", default=DEFAULT_RANK, choices=RANKS + ["all"], help="Taxonomic rank. Default: genus.")

    par = p.add_argument_group("Optional arguments - parallelism and resources")
    par.add_argument("-j", "--jobs", type=int, default=DEFAULT_JOBS, metavar="INT", help=f"Number of samples processed in parallel. Default: {DEFAULT_JOBS}.")
    par.add_argument("-p", "--threads-per-sample", type=int, default=DEFAULT_THREADS_PER_SAMPLE, metavar="INT", help=f"Per-sample CPU budget. Default: {DEFAULT_THREADS_PER_SAMPLE}.")
    par.add_argument("--seqkit-threads", default=DEFAULT_SEQKIT_THREADS, metavar="INT|auto", help="seqkit threads per sample. Default: auto.")
    par.add_argument("--seqkit-max-threads", type=int, default=DEFAULT_SEQKIT_MAX_THREADS, metavar="INT", help=f"seqkit thread cap. Default: {DEFAULT_SEQKIT_MAX_THREADS}.")
    par.add_argument("--bbduk-threads", default=DEFAULT_BBDUK_THREADS, metavar="INT|auto", help="BBDuk threads per sample. Default: auto.")
    par.add_argument("--bbduk-max-threads", type=int, default=DEFAULT_BBDUK_MAX_THREADS, metavar="INT", help=f"BBDuk thread cap. Default: {DEFAULT_BBDUK_MAX_THREADS}.")
    par.add_argument("--minimap2-threads", default=DEFAULT_MINIMAP2_THREADS, metavar="INT|auto", help="minimap2 threads per sample. Default: auto.")
    par.add_argument("--minimap2-max-threads", type=int, default=DEFAULT_MINIMAP2_MAX_THREADS, metavar="INT", help=f"minimap2 thread cap. Default: {DEFAULT_MINIMAP2_MAX_THREADS}.")
    par.add_argument("--abundance-jobs", default=DEFAULT_ABUNDANCE_JOBS, metavar="INT|auto", help="Parallel sample-level Polars abundance jobs. Default: auto.")
    par.add_argument("--polars-threads", type=int, default=DEFAULT_POLARS_THREADS, metavar="INT", help=f"Polars threads per abundance job. Default: {DEFAULT_POLARS_THREADS}.")
    par.add_argument("--bbmap-mem", default=DEFAULT_BBMAP_MEM, metavar="STR", help=f"BBDuk Java heap. Default: {DEFAULT_BBMAP_MEM}.")

    bbduk = p.add_argument_group("Optional arguments - BBDuk extraction parameters")
    for m in ["16s", "18s", "its"]:
        upper = m.upper()
        bbduk.add_argument(f"--k-{m}", type=int, default=DEFAULT_K[upper], metavar="INT", help=f"Default: {DEFAULT_K[upper]}.")
        bbduk.add_argument(f"--hdist-{m}", type=int, default=DEFAULT_HDIST[upper], metavar="INT", help=f"Default: {DEFAULT_HDIST[upper]}.")
        bbduk.add_argument(f"--mkh-{m}", type=int, default=DEFAULT_MKH[upper], metavar="INT", help=f"Default: {DEFAULT_MKH[upper]}.")
        bbduk.add_argument(f"--mink-{m}", type=int, default=DEFAULT_MINK[upper], metavar="INT", help=f"0 disables. Default: {DEFAULT_MINK[upper]}.")

    mm = p.add_argument_group("Optional arguments - minimap2 parameters")
    mm.add_argument("--minimap2-preset", default=DEFAULT_MINIMAP2_PRESET, metavar="STR", help=f"Default: {DEFAULT_MINIMAP2_PRESET}.")
    mm.add_argument("--minimap2-N", type=int, default=DEFAULT_MINIMAP2_N, metavar="INT", help=f"Default: {DEFAULT_MINIMAP2_N}.")

    filt = p.add_argument_group("Optional arguments - abundance filters")
    filt.add_argument("--min-identity-16s", type=float, default=DEFAULT_MIN_IDENTITY["16S"], metavar="FLOAT", help="Default: 0.97.")
    filt.add_argument("--min-identity-18s", type=float, default=DEFAULT_MIN_IDENTITY["18S"], metavar="FLOAT", help="Default: 0.97.")
    filt.add_argument("--min-identity-its", type=float, default=DEFAULT_MIN_IDENTITY["ITS"], metavar="FLOAT", help="Default: 0.95.")
    filt.add_argument("--min-aln-len-16s", type=int, default=DEFAULT_MIN_ALN_LEN["16S"], metavar="INT", help="Default: 80.")
    filt.add_argument("--min-aln-len-18s", type=int, default=DEFAULT_MIN_ALN_LEN["18S"], metavar="INT", help="Default: 80.")
    filt.add_argument("--min-aln-len-its", type=int, default=DEFAULT_MIN_ALN_LEN["ITS"], metavar="INT", help="Default: 80.")
    filt.add_argument("--min-qcov-16s", type=float, default=DEFAULT_MIN_QCOV["16S"], metavar="FLOAT", help="Default: 0.60.")
    filt.add_argument("--min-qcov-18s", type=float, default=DEFAULT_MIN_QCOV["18S"], metavar="FLOAT", help="Default: 0.60.")
    filt.add_argument("--min-qcov-its", type=float, default=DEFAULT_MIN_QCOV["ITS"], metavar="FLOAT", help="Default: 0.60.")
    filt.add_argument("--min-mapq", type=int, default=DEFAULT_MIN_MAPQ, metavar="INT", help="Default: 0.")

    run = p.add_argument_group("Optional arguments - resume, overwrite and logging")
    run.add_argument("--read-count-method", choices=["auto", "seqkit", "python"], default="auto", help="Default: auto.")
    run.add_argument("--force", action="store_true", help="Re-run and overwrite outputs.")
    run.add_argument("--no-resume", action="store_true", help="Ignore checkpoints and output checks.")
    run.add_argument("--clean-tmp", action="store_true", help="Remove temporary files after successful run.")
    run.add_argument("--dry-run", action="store_true", help="Check inputs and print run plan only.")
    run.add_argument("--check-deps", action="store_true", help="Check dependencies and exit.")
    run.add_argument("--keep-task-logs", action="store_true", help="Keep per-sample logs even when tasks succeed.")
    run.add_argument("--quiet", action="store_true", help="Suppress banner and summary logs.")
    run.add_argument("--verbose", action="store_true", help="Show debug logs.")

    helpg = p.add_argument_group("Help")
    helpg.add_argument("-h", "--help", action="help", help="Show this help message and exit.")

    return p


def validate_cli(args: argparse.Namespace) -> None:
    dependency_check_only = args.check_deps and not (args.input or args.sample_id or args.r1 or args.r2)
    if dependency_check_only:
        return

    if args.input:
        if args.sample_id or args.r1 or args.r2:
            raise ValueError("Use either --input or --sample-id/--r1/--r2, not both.")
        if not Path(args.input).is_file():
            raise FileNotFoundError(f"Input TSV not found: {args.input}")
    else:
        if not (args.sample_id and args.r1 and args.r2):
            raise ValueError("Missing input. Use --input data_path.tsv or --sample-id with --r1 and --r2.")

    for name in ["jobs", "threads_per_sample", "seqkit_max_threads", "bbduk_max_threads", "minimap2_max_threads", "polars_threads"]:
        if int(getattr(args, name)) < 1:
            raise ValueError(f"--{name.replace('_', '-')} must be positive.")

    if args.abundance_jobs != "auto" and int(args.abundance_jobs) < 1:
        raise ValueError("--abundance-jobs must be positive or auto.")


def cleanup_tmp(paths: Paths, logger: logging.Logger) -> None:
    logger.info("[yellow]Cleaning temporary files[/yellow]")
    for d in [Path(paths.marker_dir) / "tmp", Path(paths.align_dir) / "tmp", Path(paths.abundance_tmp_dir)]:
        if d.exists():
            shutil.rmtree(d, ignore_errors=True)
            ensure_dir(d)


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    validate_cli(args)

    paths = init_paths(args.outdir)
    console = make_console()
    logger = setup_logger(console, paths.main_log, args.verbose)

    dependency_check_only = args.check_deps and not (args.input or args.sample_id or args.r1 or args.r2)
    if dependency_check_only:
        samples, extra_cols = [], []
    elif args.input:
        samples, extra_cols = read_samples_from_input(args.input)
    else:
        samples, extra_cols = build_single_sample(args)

    if samples:
        check_inputs(samples)
        write_manifest(samples, extra_cols, paths)

    cfg = read_config(args, paths, logger)
    check_dependencies(cfg, logger)
    write_run_config(cfg, paths)

    if not args.quiet:
        print_run_plan(console, cfg, paths, samples, extra_cols)

    if args.check_deps:
        logger.info("[green]Dependency check finished.[/green]")
        return

    if cfg.dry_run:
        logger.info("[yellow]Dry run finished. No workflow steps executed.[/yellow]")
        return

    t0 = time.time()

    if "reads_count" in cfg.steps:
        run_reads_count(samples, cfg, paths, logger)
    if "extract" in cfg.steps:
        run_extract(samples, cfg, paths, logger)
    if "align" in cfg.steps:
        run_align(samples, cfg, paths, logger)
    if "abundance" in cfg.steps:
        run_abundance(samples, cfg, paths, extra_cols, logger)

    if cfg.clean_tmp:
        cleanup_tmp(paths, logger)

    elapsed = time.time() - t0
    if not args.quiet:
        table = Table.grid(padding=(0, 1))
        table.add_column(style="bold cyan", justify="right")
        table.add_column(style="white")
        table.add_row("Samples", f"{len(samples):,}")
        table.add_row("Markers", ",".join(cfg.markers))
        table.add_row("Steps", ",".join(cfg.steps))
        table.add_row("Elapsed", f"{elapsed:.2f} sec")
        table.add_row("Outdir", paths.outdir)
        table.add_row("Main log", paths.main_log)
        console.print(Panel(table, title="[bold green]Finished[/bold green]", border_style="green", expand=False))

def run() -> int:
    try:
        main()
        return 0
    except BrokenPipeError:
        return 0
    except KeyboardInterrupt:
        Console(stderr=True).print("[yellow]Interrupted by user.[/yellow]")
        return 130
    except Exception as exc:
        logger = logging.getLogger("metamarker_profile")
        if logger.handlers:
            logger.exception("Workflow failed")
        Console(stderr=True).print(f"[bold red]ERROR:[/bold red] {exc}")
        return 1


if __name__ == "__main__":
    raise SystemExit(run())
