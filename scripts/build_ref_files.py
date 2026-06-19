#!/usr/bin/env python3
"""Build short-id marker references and taxonomy tables from SILVA/UNITE FASTA.

This script prepares reference files for marker-count workflows:

  * SILVA SSURef -> DNA FASTA split into arc/bac 16S and eukaryotic 18S
  * UNITE ITS    -> short-id ITS FASTA
  * taxonomy TSV -> ref_id, marker, domain, phylum, class, order, family, genus, species
  * minimap2 indexes, optionally built in parallel

Features:

  * Rich progress bars and colored logs
  * grouped, readable command-line help
  * timestamped log file
  * threaded record preparation
  * parallel minimap2 index building
  * explicit organelle filtering summary
"""

from __future__ import annotations

import argparse
import gzip
import logging
import os
import re
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Iterator, TextIO, NoReturn

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
from rich_argparse import RichHelpFormatter


PROGRAM = "build_ref_files.py"
MISSING = "Unclassified"
SOFTWARE_DIR = Path(__file__).resolve().parents[1]


def default_ref_dir() -> Path:
    env_ref = os.environ.get("META_MARKER_COUNT_REF_DIR")
    if env_ref:
        return Path(env_ref)

    config_path = Path(
        os.environ.get("META_MARKER_COUNT_REF_CONFIG", SOFTWARE_DIR / ".meta_marker_count_ref_dir")
    )
    if config_path.is_file():
        for raw_line in config_path.read_text(encoding="utf-8").splitlines():
            line = raw_line.split("#", 1)[0].strip()
            if not line:
                continue
            path = Path(line)
            return path if path.is_absolute() else SOFTWARE_DIR / path

    return SOFTWARE_DIR / "refs"


DEFAULT_REF_DIR = default_ref_dir()

TAX_HEADER = [
    "ref_id",
    "marker",
    "domain",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
]

RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]

PREFIX_RANK = {
    "d": "domain",
    "k": "domain",
    "p": "phylum",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species",
}

QIIME_RANK = {
    "D_0": "domain",
    "D_1": "phylum",
    "D_2": "class",
    "D_3": "order",
    "D_4": "family",
    "D_5": "genus",
    "D_6": "species",
}

NAME_RANK = {
    "domain": "domain",
    "superkingdom": "domain",
    "kingdom": "domain",
    "phylum": "phylum",
    "class": "class",
    "order": "order",
    "family": "family",
    "genus": "genus",
    "species": "species",
}


@dataclass(frozen=True)
class Taxonomy:
    domain: str = MISSING
    phylum: str = MISSING
    klass: str = MISSING
    order: str = MISSING
    family: str = MISSING
    genus: str = MISSING
    species: str = MISSING

    def values(self) -> list[str]:
        return [
            self.domain,
            self.phylum,
            self.klass,
            self.order,
            self.family,
            self.genus,
            self.species,
        ]


@dataclass(frozen=True)
class RefTaxonomy:
    ref_id: str
    marker: str
    taxonomy: Taxonomy

    def row(self) -> list[str]:
        return [self.ref_id, self.marker, *self.taxonomy.values()]


@dataclass
class ProcessStats:
    source: str
    total_seen: int = 0
    written: int = 0
    arc_bac: int = 0
    euk: int = 0
    unite: int = 0
    organelle_skipped: int = 0
    skipped_domain: int = 0
    duplicate_ids: int = 0
    elapsed_sec: float = 0.0
    output_files: list[Path] = field(default_factory=list)


@dataclass(frozen=True)
class PreparedRecord:
    target: str
    record: RefTaxonomy | None
    sequence: str
    skip_reason: str = ""


class IdRegistry:
    """Create unique short IDs while preserving the original short token when possible."""

    def __init__(self) -> None:
        self._counts: dict[str, int] = {}
        self.duplicate_count = 0

    def unique(self, raw_id: str) -> str:
        base = sanitize_ref_id(raw_id)
        if not base:
            base = "ref"
        n = self._counts.get(base, 0) + 1
        self._counts[base] = n
        if n == 1:
            return base
        self.duplicate_count += 1
        return f"{base}_{n}"


class PlainProgress:
    """Tiny fallback for environments without rich."""

    def __init__(self, logger: logging.Logger, log_every: int = 100000) -> None:
        self.logger = logger
        self.log_every = max(1, log_every)
        self.count = 0
        self.description = ""

    def __enter__(self) -> "PlainProgress":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:  # noqa: ANN001
        return None

    def add_task(self, description: str, total: int | None = None) -> int:
        self.description = description
        self.count = 0
        self.logger.info("%s started", description)
        return 1

    def update(self, task_id: int, advance: int = 0, description: str | None = None) -> None:
        del task_id
        if description:
            self.description = description
        self.count += advance
        if self.count and self.count % self.log_every == 0:
            self.logger.info("%s: %s records", self.description, f"{self.count:,}")


def now_stamp() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


def resolve_display_width(
    display_width: int | None = None,
    display_width_ratio: float = 0.5,
    min_width: int = 80,
) -> int:
    if display_width is not None and display_width > 0:
        return display_width

    terminal_width = shutil.get_terminal_size(fallback=(160, 24)).columns
    return max(min_width, int(terminal_width * display_width_ratio))


def make_console(
    no_rich: bool = False,
    display_width: int | None = None,
    display_width_ratio: float = 0.5,
):  # noqa: ANN201
    if no_rich:
        return None

    width = resolve_display_width(
        display_width=display_width,
        display_width_ratio=display_width_ratio,
    )
    return Console(stderr=True, width=width)


def setup_logging(outdir: Path, log_file: Path | None, quiet: bool, no_rich: bool) -> logging.Logger:
    outdir.mkdir(parents=True, exist_ok=True)
    if log_file is None:
        log_dir = outdir / "logs"
        log_dir.mkdir(parents=True, exist_ok=True)
        log_file = log_dir / f"build_ref_files.{time.strftime('%Y%m%d_%H%M%S')}.log"
    else:
        log_file.parent.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger("build_ref_files")
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()
    logger.propagate = False

    file_handler = logging.FileHandler(log_file, encoding="utf-8")
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(
        logging.Formatter(
            fmt="[%(asctime)s] [%(levelname)s] %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
    )
    logger.addHandler(file_handler)

    if not quiet:
        if not no_rich:
            console_handler = RichHandler(
                rich_tracebacks=True,
                show_time=True,
                show_level=True,
                show_path=False,
                markup=True,
            )
            console_handler.setLevel(logging.INFO)
            console_handler.setFormatter(logging.Formatter("%(message)s"))
        else:
            console_handler = logging.StreamHandler(sys.stderr)
            console_handler.setLevel(logging.INFO)
            console_handler.setFormatter(
                logging.Formatter(
                    fmt="[%(asctime)s] [%(levelname)s] %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S",
                )
            )
        logger.addHandler(console_handler)

    logger.info("Log file: %s", log_file)
    return logger


def create_progress(
    console,
    logger: logging.Logger,
    log_every: int,
    progress_bar_width: int | None = None,
):  # noqa: ANN001, ANN201
    if console is not None:
        if progress_bar_width is None:
            progress_bar_width = max(20, min(60, int(console.width * 0.35)))
        return Progress(
            SpinnerColumn(),
            TextColumn("[bold blue]{task.description}"),
            BarColumn(bar_width=progress_bar_width),
            TaskProgressColumn(),
            MofNCompleteColumn(),
            TimeElapsedColumn(),
            TimeRemainingColumn(),
            console=console,
            transient=False,
            expand=False,
        )
    return PlainProgress(logger=logger, log_every=log_every)


def open_text(path: Path) -> TextIO:
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def count_fasta_records(path: Path) -> int:
    count = 0
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
    return count


def read_fasta(path: Path) -> Iterator[tuple[str, str]]:
    header: str | None = None
    seq_parts: list[str] = []
    with open_text(path) as handle:
        for line in handle:
            line = line.rstrip("\n\r")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line.strip())
    if header is not None:
        yield header, "".join(seq_parts)


def batched(iterator: Iterator[tuple[str, str]], batch_size: int) -> Iterator[list[tuple[str, str]]]:
    batch: list[tuple[str, str]] = []
    for item in iterator:
        batch.append(item)
        if len(batch) >= batch_size:
            yield batch
            batch = []
    if batch:
        yield batch


def strip_fasta_suffix(path: Path) -> str:
    name = path.name
    for suffix in (
        ".fasta.gz",
        ".fa.gz",
        ".fna.gz",
        ".fas.gz",
        ".fasta",
        ".fa",
        ".fna",
        ".fas",
        ".gz",
    ):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def silva_output_prefix(path: Path) -> str:
    """Return a SILVA output prefix with exactly one trailing .dna."""
    base = strip_fasta_suffix(path)
    if base.endswith(".dna"):
        return base
    return f"{base}.dna"


def expected_output_paths(args: argparse.Namespace) -> list[Path]:
    """List output files that may be created by this run before any heavy work starts."""
    outdir = args.outdir
    paths: list[Path] = []

    if args.silva:
        base = silva_output_prefix(args.silva)
        silva_fastas = [
            outdir / f"{base}.arc_bac.shortid.fasta",
            outdir / f"{base}.euk.shortid.fasta",
        ]
        paths.extend(silva_fastas)
        paths.extend([
            outdir / f"{base}.arc_bac.taxonomy.tsv",
            outdir / f"{base}.euk.taxonomy.tsv",
        ])
        if not args.skip_index:
            paths.extend(Path(f"{fasta}.mmi") for fasta in silva_fastas)

    if args.unite:
        base = strip_fasta_suffix(args.unite)
        unite_fasta = outdir / f"{base}.shortid.fasta"
        paths.extend([
            unite_fasta,
            outdir / f"{base}.taxonomy.tsv",
        ])
        if not args.skip_index:
            paths.append(Path(f"{unite_fasta}.mmi"))

    if not args.skip_combined:
        paths.append(args.combined_taxonomy or (outdir / "ref_taxonomy.tsv"))

    return paths


def preflight_output_conflicts(args: argparse.Namespace) -> None:
    if args.force:
        return
    conflicts = [path for path in expected_output_paths(args) if path.exists()]
    if not conflicts:
        return

    shown = "\n".join(f"  - {path}" for path in conflicts[:20])
    more = "" if len(conflicts) <= 20 else f"\n  ... and {len(conflicts) - 20} more"
    die(
        "Output file(s) already exist. Use --force to overwrite, "
        "or choose a new --outdir.\n"
        f"{shown}{more}"
    )


def sanitize_ref_id(raw_id: str) -> str:
    raw_id = raw_id.strip()
    raw_id = raw_id.split()[0] if raw_id else ""
    raw_id = raw_id.strip(">|")
    raw_id = re.sub(r"[^A-Za-z0-9_.:-]+", "_", raw_id)
    raw_id = re.sub(r"_+", "_", raw_id)
    return raw_id.strip("_")


def short_id_from_header(header: str) -> str:
    token = header.split(None, 1)[0] if header else ""
    if "|" in token:
        for part in token.split("|"):
            part = part.strip()
            if part and not part.endswith("__") and ";" not in part:
                return part
    return token


def normalize_sequence(seq: str) -> str:
    return re.sub(r"\s+", "", seq).upper().replace("U", "T")


def write_fasta_record(handle: TextIO, ref_id: str, seq: str, width: int = 80) -> None:
    handle.write(f">{ref_id}\n")
    seq = normalize_sequence(seq)
    for i in range(0, len(seq), width):
        handle.write(seq[i : i + width] + "\n")


def split_header_description(header: str) -> str:
    parts = header.split(None, 1)
    return parts[1].strip() if len(parts) == 2 else ""


def taxonomy_candidate(header: str) -> str:
    candidates: list[str] = []
    desc = split_header_description(header)
    if desc:
        candidates.append(desc)
    if "|" in header:
        candidates.extend(part.strip() for part in header.split("|") if part.strip())

    scored: list[tuple[int, int, str]] = []
    for cand in candidates:
        score = cand.count(";")
        prefix_bonus = 1 if re.search(r"(^|[;\s|])(?:d|k|p|c|o|f|g|s)__", cand) else 0
        scored.append((score, prefix_bonus, cand))
    scored.sort(reverse=True)
    if scored and scored[0][0] > 0:
        return scored[0][2]
    return desc


def taxonomy_parts(header: str) -> list[str]:
    cand = taxonomy_candidate(header)
    return [p.strip() for p in cand.strip().strip(";").split(";") if p.strip()]


def detect_rank(raw: str) -> tuple[str | None, str]:
    value = raw.strip()

    m = re.match(r"^(D_[0-6])__(.*)$", value, flags=re.IGNORECASE)
    if m:
        return QIIME_RANK.get(m.group(1).upper()), m.group(2)

    m = re.match(r"^([dkpcofgs])__(.*)$", value, flags=re.IGNORECASE)
    if m:
        return PREFIX_RANK.get(m.group(1).lower()), m.group(2)

    m = re.match(r"^([A-Za-z_]+):(.*)$", value)
    if m:
        return NAME_RANK.get(m.group(1).lower()), m.group(2)

    return None, value


def clean_taxon(raw: str) -> str:
    _rank, value = detect_rank(raw)
    value = value.strip().strip(";")
    value = value.replace(" ", "_")
    value = re.sub(r"\s+", "_", value)
    value = re.sub(r"^_+|_+$", "", value)
    if not value or value in {"__", "_", "NA", "N/A", "nan", "NaN"}:
        return MISSING
    return value


def parse_taxonomy(header: str, default_domain: str = MISSING) -> Taxonomy:
    parts = taxonomy_parts(header)
    ranks: dict[str, str] = {}
    sequential: list[str] = []

    for part in parts:
        rank, _value = detect_rank(part)
        value = clean_taxon(part)
        if rank:
            ranks[rank] = value
        else:
            sequential.append(value)

    if ranks:
        values = {rank: ranks.get(rank, MISSING) for rank in RANKS}
    else:
        values = {rank: MISSING for rank in RANKS}
        for rank, value in zip(RANKS, sequential):
            values[rank] = value

    if values["domain"] == MISSING and default_domain != MISSING:
        values["domain"] = default_domain

    return Taxonomy(
        domain=values["domain"],
        phylum=values["phylum"],
        klass=values["class"],
        order=values["order"],
        family=values["family"],
        genus=values["genus"],
        species=values["species"],
    )


def taxonomy_from_values(values: dict[str, str]) -> Taxonomy:
    return Taxonomy(
        domain=values.get("domain", MISSING),
        phylum=values.get("phylum", MISSING),
        klass=values.get("class", MISSING),
        order=values.get("order", MISSING),
        family=values.get("family", MISSING),
        genus=values.get("genus", MISSING),
        species=values.get("species", MISSING),
    )


def domain_key(domain: str) -> str:
    return domain.lower().replace("_", "").replace(" ", "")


def parse_silva_taxonomy(header: str) -> Taxonomy:
    raw_parts = taxonomy_parts(header)
    parts = [clean_taxon(part) for part in raw_parts]
    if not parts:
        return parse_taxonomy(header)

    dkey = domain_key(parts[0])
    values = {rank: MISSING for rank in RANKS}
    values["domain"] = parts[0]

    if dkey in {"bacteria", "archaea"}:
        for rank, value in zip(RANKS, parts):
            values[rank] = value
        return taxonomy_from_values(values)

    if dkey == "eukaryota":
        lineage = parts[1:]
        upper = lineage[:-2] if len(lineage) >= 2 else lineage
        for rank, value in zip(["phylum", "class", "order", "family"], upper):
            values[rank] = value
        if len(lineage) >= 2:
            values["genus"] = lineage[-2]
            values["species"] = lineage[-1]
        elif len(lineage) == 1:
            values["species"] = lineage[0]
        return taxonomy_from_values(values)

    return parse_taxonomy(header)


def has_silva_organelle_lineage(header: str, keywords: tuple[str, ...]) -> bool:
    parts = [part.strip().lower() for part in taxonomy_parts(header)]
    text = ";".join(parts)
    return any(keyword.lower() in text for keyword in keywords)


def extract_unite_sh_id(header: str) -> str | None:
    match = re.search(r"\bSH\d+\.\d+[A-Za-z]*\b", header)
    return match.group(0) if match else None


def parse_rank_names(raw: str) -> tuple[str, ...]:
    rank_names = tuple(x.strip().lower() for x in re.split(r"[,;]", raw) if x.strip())
    allowed = set(RANKS)
    invalid = [rank for rank in rank_names if rank not in allowed]
    if invalid:
        raise argparse.ArgumentTypeError(
            "Unsupported rank(s): " + ", ".join(invalid) + ". Allowed: " + ", ".join(RANKS)
        )
    return rank_names


def is_unite_placeholder_taxon(value: str, rank: str) -> bool:
    text = value.strip().lower()
    compact = re.sub(r"[^a-z0-9]+", "_", text).strip("_")

    if not compact or compact in {"na", "n_a", "nan", "unclassified", "unidentified"}:
        return True
    if "incertae" in compact or "sedis" in compact:
        return True
    if "unclassified" in compact or "unidentified" in compact:
        return True
    if rank == "species" and (compact.endswith("_sp") or compact.endswith("_sp_") or compact == "sp"):
        return True
    return False


def add_unite_sh_suffix_to_taxonomy(
    tax: Taxonomy,
    sh_id: str | None,
    ranks_to_suffix: tuple[str, ...],
) -> Taxonomy:
    if not sh_id:
        return tax

    values = dict(zip(RANKS, tax.values(), strict=True))
    for rank in ranks_to_suffix:
        if rank == "domain":
            continue
        value = values.get(rank, MISSING)
        if "|SH" in value:
            continue
        if is_unite_placeholder_taxon(value, rank):
            values[rank] = f"{value}|{sh_id}"

    return taxonomy_from_values(values)


def write_tax_header(handle: TextIO) -> None:
    handle.write("\t".join(TAX_HEADER) + "\n")


def write_tax_row(handle: TextIO, record: RefTaxonomy) -> None:
    handle.write("\t".join(record.row()) + "\n")


def die(message: str) -> NoReturn:
    print(f"[ERROR] {message}", file=sys.stderr)
    raise SystemExit(1)


def ensure_can_write(path: Path, force: bool) -> None:
    if path.exists() and not force:
        die(f"Refusing to overwrite existing file without --force: {path}")


def prepare_unite_record(
    header: str,
    seq: str,
    ref_id: str,
    default_domain: str,
    add_sh_suffix: bool,
    sh_suffix_ranks: tuple[str, ...],
) -> PreparedRecord:
    tax = parse_taxonomy(header, default_domain=default_domain)
    if add_sh_suffix:
        tax = add_unite_sh_suffix_to_taxonomy(
            tax,
            extract_unite_sh_id(header),
            sh_suffix_ranks,
        )
    record = RefTaxonomy(ref_id=ref_id, marker="ITS", taxonomy=tax)
    return PreparedRecord(target="unite", record=record, sequence=normalize_sequence(seq))


def prepare_silva_record(
    header: str,
    seq: str,
    ref_id: str,
    keep_organelles: bool,
    organelle_keywords: tuple[str, ...],
) -> PreparedRecord:
    if not keep_organelles and has_silva_organelle_lineage(header, organelle_keywords):
        return PreparedRecord(target="skip", record=None, sequence="", skip_reason="organelle")

    tax = parse_silva_taxonomy(header)
    dkey = domain_key(tax.domain)

    if dkey in {"bacteria", "archaea"}:
        record = RefTaxonomy(ref_id=ref_id, marker="16S", taxonomy=tax)
        return PreparedRecord(target="arc_bac", record=record, sequence=normalize_sequence(seq))

    if dkey == "eukaryota":
        record = RefTaxonomy(ref_id=ref_id, marker="18S", taxonomy=tax)
        return PreparedRecord(target="euk", record=record, sequence=normalize_sequence(seq))

    return PreparedRecord(target="skip", record=None, sequence="", skip_reason="domain")


def write_prepared_record(
    prepared: PreparedRecord,
    fasta_handles: dict[str, TextIO],
    tax_handles: dict[str, TextIO],
    all_records: list[RefTaxonomy],
    stats: ProcessStats,
) -> None:
    if prepared.target == "skip":
        if prepared.skip_reason == "organelle":
            stats.organelle_skipped += 1
        else:
            stats.skipped_domain += 1
        return

    if prepared.record is None:
        return

    handle_key = prepared.target
    write_fasta_record(fasta_handles[handle_key], prepared.record.ref_id, prepared.sequence)
    write_tax_row(tax_handles[handle_key], prepared.record)
    all_records.append(prepared.record)
    stats.written += 1

    if prepared.target == "arc_bac":
        stats.arc_bac += 1
    elif prepared.target == "euk":
        stats.euk += 1
    elif prepared.target == "unite":
        stats.unite += 1


def maybe_count_records(path: Path, count_records: bool, logger: logging.Logger) -> int | None:
    if not count_records:
        return None
    logger.info("Counting FASTA records for exact progress: %s", path)
    start = time.time()
    total = count_fasta_records(path)
    logger.info("Detected %s FASTA records in %.1f sec", f"{total:,}", time.time() - start)
    return total


def process_unite(
    fasta: Path,
    outdir: Path,
    id_registry: IdRegistry,
    default_domain: str,
    add_sh_suffix: bool,
    sh_suffix_ranks: tuple[str, ...],
    force: bool,
    threads: int,
    batch_size: int,
    progress,  # noqa: ANN001
    count_records: bool,
    logger: logging.Logger,
) -> tuple[list[RefTaxonomy], list[Path], ProcessStats]:
    start = time.time()
    base = strip_fasta_suffix(fasta)
    fasta_out = outdir / f"{base}.shortid.fasta"
    tax_out = outdir / f"{base}.taxonomy.tsv"
    for path in (fasta_out, tax_out):
        ensure_can_write(path, force)

    logger.info("Processing UNITE: %s", fasta)
    total = maybe_count_records(fasta, count_records, logger)
    task_id = progress.add_task(f"UNITE ITS: {fasta.name}", total=total)

    records: list[RefTaxonomy] = []
    stats = ProcessStats(source="UNITE", output_files=[fasta_out, tax_out])

    with fasta_out.open("wt", encoding="utf-8") as fasta_handle, tax_out.open(
        "wt", encoding="utf-8"
    ) as tax_handle:
        write_tax_header(tax_handle)
        fasta_handles = {"unite": fasta_handle}
        tax_handles = {"unite": tax_handle}

        with ThreadPoolExecutor(max_workers=threads) as executor:
            for batch in batched(read_fasta(fasta), batch_size):
                futures = []
                for header, seq in batch:
                    ref_id = id_registry.unique(short_id_from_header(header))
                    futures.append(
                        executor.submit(
                            prepare_unite_record,
                            header,
                            seq,
                            ref_id,
                            default_domain,
                            add_sh_suffix,
                            sh_suffix_ranks,
                        )
                    )

                for future in futures:
                    prepared = future.result()
                    write_prepared_record(prepared, fasta_handles, tax_handles, records, stats)
                    stats.total_seen += 1
                    progress.update(task_id, advance=1)

    stats.duplicate_ids = id_registry.duplicate_count
    stats.elapsed_sec = time.time() - start
    logger.info("UNITE records written: %s", f"{stats.unite:,}")
    logger.info("UNITE FASTA: %s", fasta_out)
    logger.info("UNITE taxonomy: %s", tax_out)
    return records, [fasta_out], stats


def process_silva(
    fasta: Path,
    outdir: Path,
    id_registry: IdRegistry,
    force: bool,
    keep_organelles: bool,
    organelle_keywords: tuple[str, ...],
    threads: int,
    batch_size: int,
    progress,  # noqa: ANN001
    count_records: bool,
    logger: logging.Logger,
) -> tuple[list[RefTaxonomy], list[Path], ProcessStats]:
    start = time.time()
    base = silva_output_prefix(fasta)
    arc_bac_fasta = outdir / f"{base}.arc_bac.shortid.fasta"
    arc_bac_tax = outdir / f"{base}.arc_bac.taxonomy.tsv"
    euk_fasta = outdir / f"{base}.euk.shortid.fasta"
    euk_tax = outdir / f"{base}.euk.taxonomy.tsv"
    for path in (arc_bac_fasta, arc_bac_tax, euk_fasta, euk_tax):
        ensure_can_write(path, force)

    logger.info("Processing SILVA: %s", fasta)
    if keep_organelles:
        logger.warning("SILVA organelle records are kept because --keep-organelles was used")
    else:
        logger.info("SILVA organelle filtering keywords: %s", ", ".join(organelle_keywords))

    total = maybe_count_records(fasta, count_records, logger)
    task_id = progress.add_task(f"SILVA SSU: {fasta.name}", total=total)

    records: list[RefTaxonomy] = []
    stats = ProcessStats(
        source="SILVA",
        output_files=[arc_bac_fasta, arc_bac_tax, euk_fasta, euk_tax],
    )

    with arc_bac_fasta.open("wt", encoding="utf-8") as ab_fa, arc_bac_tax.open(
        "wt", encoding="utf-8"
    ) as ab_tax, euk_fasta.open("wt", encoding="utf-8") as euk_fa, euk_tax.open(
        "wt", encoding="utf-8"
    ) as eu_tax:
        write_tax_header(ab_tax)
        write_tax_header(eu_tax)
        fasta_handles = {"arc_bac": ab_fa, "euk": euk_fa}
        tax_handles = {"arc_bac": ab_tax, "euk": eu_tax}

        with ThreadPoolExecutor(max_workers=threads) as executor:
            for batch in batched(read_fasta(fasta), batch_size):
                futures = []
                for header, seq in batch:
                    ref_id = id_registry.unique(short_id_from_header(header))
                    futures.append(
                        executor.submit(
                            prepare_silva_record,
                            header,
                            seq,
                            ref_id,
                            keep_organelles,
                            organelle_keywords,
                        )
                    )

                for future in futures:
                    prepared = future.result()
                    write_prepared_record(prepared, fasta_handles, tax_handles, records, stats)
                    stats.total_seen += 1
                    progress.update(task_id, advance=1)

    stats.duplicate_ids = id_registry.duplicate_count
    stats.elapsed_sec = time.time() - start
    logger.info("SILVA arc/bac 16S records written: %s", f"{stats.arc_bac:,}")
    logger.info("SILVA euk 18S records written: %s", f"{stats.euk:,}")
    logger.info("SILVA organelle records skipped: %s", f"{stats.organelle_skipped:,}")
    logger.info("SILVA records skipped because domain was not Bacteria/Archaea/Eukaryota: %s", f"{stats.skipped_domain:,}")
    logger.info("SILVA arc/bac FASTA: %s", arc_bac_fasta)
    logger.info("SILVA arc/bac taxonomy: %s", arc_bac_tax)
    logger.info("SILVA euk FASTA: %s", euk_fasta)
    logger.info("SILVA euk taxonomy: %s", euk_tax)
    return records, [arc_bac_fasta, euk_fasta], stats


def write_combined_taxonomy(path: Path, records: Iterable[RefTaxonomy], force: bool, logger: logging.Logger) -> int:
    ensure_can_write(path, force)
    count = 0
    with path.open("wt", encoding="utf-8") as handle:
        write_tax_header(handle)
        for record in records:
            write_tax_row(handle, record)
            count += 1
    logger.info("Combined taxonomy records written: %s", f"{count:,}")
    logger.info("Combined taxonomy: %s", path)
    return count


def build_one_minimap2_index(
    fasta: Path,
    force: bool,
    minimap2_path: str,
    logger: logging.Logger,
) -> Path | None:
    if not fasta.exists() or fasta.stat().st_size == 0:
        logger.warning("Skipping minimap2 index for empty FASTA: %s", fasta)
        return None

    index = Path(f"{fasta}.mmi")
    ensure_can_write(index, force)
    cmd = [minimap2_path, "-d", str(index), str(fasta)]
    logger.info("Building minimap2 index: %s", index)
    logger.debug("CMD: %s", " ".join(map(str, cmd)))
    subprocess.run(cmd, check=True)
    return index


def build_minimap2_indexes(
    fastas: list[Path],
    force: bool,
    index_jobs: int,
    minimap2_path: str,
    progress,  # noqa: ANN001
    logger: logging.Logger,
) -> list[Path]:
    if not fastas:
        return []

    mm2 = shutil.which(minimap2_path) if os.sep not in minimap2_path else minimap2_path
    if mm2 is None:
        die("minimap2 not found in PATH. Install it or rerun with --skip-index.")

    task_id = progress.add_task("minimap2 indexing", total=len(fastas))
    built: list[Path] = []
    workers = max(1, min(index_jobs, len(fastas)))

    with ThreadPoolExecutor(max_workers=workers) as executor:
        future_to_fasta = {
            executor.submit(build_one_minimap2_index, fasta, force, str(mm2), logger): fasta
            for fasta in fastas
        }
        for future in as_completed(future_to_fasta):
            fasta = future_to_fasta[future]
            index = future.result()
            if index is not None:
                built.append(index)
            progress.update(task_id, advance=1)
    return built


def parse_organelle_keywords(raw: str) -> tuple[str, ...]:
    keywords = tuple(k.strip().lower() for k in re.split(r"[,;]", raw) if k.strip())
    if not keywords:
        raise argparse.ArgumentTypeError("At least one organelle keyword is required.")
    return keywords


def positive_int(value: str) -> int:
    if not re.fullmatch(r"[1-9][0-9]*", str(value)):
        raise argparse.ArgumentTypeError(f"Expected positive integer, got: {value}")
    return int(value)

def positive_float(value: str) -> float:
    text = str(value)
    if not re.fullmatch(r"(0?\.[0-9]+|[1-9][0-9]*(\.[0-9]+)?)", text):
        raise argparse.ArgumentTypeError(f"Expected positive float, got: {value}")
    x = float(text)
    if x <= 0:
        raise argparse.ArgumentTypeError(f"Expected positive float, got: {value}")
    return x


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=PROGRAM,
        description="Generate short-id FASTA and taxonomy TSV files for marker-count workflows.",
        formatter_class=RichHelpFormatter,
    )
    input_group = parser.add_argument_group("Input databases")
    input_group.add_argument("--silva", type=Path, help="SILVA FASTA or FASTA.GZ file.")
    input_group.add_argument("--unite", type=Path, help="UNITE FASTA or FASTA.GZ file.")

    output_group = parser.add_argument_group("Output files")
    output_group.add_argument(
        "-o",
        "--outdir",
        type=Path,
        default=DEFAULT_REF_DIR,
        help="Output directory.",
    )
    output_group.add_argument(
        "--combined-taxonomy",
        type=Path,
        help="Combined taxonomy TSV path. Default: OUTDIR/ref_taxonomy.tsv",
    )
    output_group.add_argument(
        "--skip-combined",
        action="store_true",
        help="Only write per-database FASTA/taxonomy files.",
    )
    output_group.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing output files.",
    )

    taxonomy_group = parser.add_argument_group("Taxonomy parsing and filtering")
    taxonomy_group.add_argument(
        "--unite-default-domain",
        default="Fungi",
        help="Domain used for UNITE records whose header has no parseable domain.",
    )
    taxonomy_group.add_argument(
        "--no-unite-sh-suffix",
        action="store_true",
        help="Do not append SH IDs to unresolved UNITE rank labels.",
    )
    taxonomy_group.add_argument(
        "--unite-sh-suffix-ranks",
        type=parse_rank_names,
        default=("phylum", "class", "order", "family", "genus", "species"),
        help="Comma-separated UNITE ranks that receive SH suffixes when unresolved.",
    )
    taxonomy_group.add_argument(
        "--keep-organelles",
        action="store_true",
        help="Keep SILVA records containing organelle keywords. Default skips them.",
    )
    taxonomy_group.add_argument(
        "--organelle-keywords",
        type=parse_organelle_keywords,
        default=("Mitochondria", "Chloroplast", "Plastid"),
        help="Comma-separated keywords used to detect SILVA organelle records.",
    )

    index_group = parser.add_argument_group("Index building")
    index_group.add_argument(
        "--skip-index",
        action="store_true",
        help="Do not build minimap2 .mmi indexes.",
    )
    index_group.add_argument(
        "--minimap2",
        default="minimap2",
        help="minimap2 executable name or path.",
    )
    index_group.add_argument(
        "--index-jobs",
        type=positive_int,
        default=2,
        help="Number of FASTA index-building jobs to run in parallel.",
    )

    performance_group = parser.add_argument_group("Performance")
    performance_group.add_argument(
        "-t",
        "--threads",
        type=positive_int,
        default=max(1, min(4, os.cpu_count() or 1)),
        help="Worker threads for FASTA record preparation.",
    )
    performance_group.add_argument(
        "--batch-size",
        type=positive_int,
        default=5000,
        help="Number of FASTA records submitted to worker threads per batch.",
    )
    performance_group.add_argument(
        "--count-records",
        action="store_true",
        help="Pre-scan FASTA files to show exact progress totals. This costs an extra pass.",
    )

    log_group = parser.add_argument_group("Logging and display")
    log_group.add_argument(
        "--log-file",
        type=Path,
        help="Log file path. Default: OUTDIR/logs/build_ref_files.YYYYmmdd_HHMMSS.log",
    )
    log_group.add_argument(
        "--quiet",
        action="store_true",
        help="Write logs to file only.",
    )
    log_group.add_argument(
        "--no-rich",
        action="store_true",
        help="Disable Rich colored output and progress bars.",
    )
    log_group.add_argument(
        "--log-every",
        type=positive_int,
        default=100000,
        help="Progress logging interval for file logs.",
    )
    log_group.add_argument(
        "--display-width",
        type=positive_int,
        default=None,
        help=(
            "Fixed Rich display width. Default: auto, current terminal width "
            "multiplied by --display-width-ratio."
        ),
    )
    log_group.add_argument(
        "--display-width-ratio",
        type=positive_float,
        default=0.5,
        help="Rich display width ratio relative to terminal width. Default: 0.5.",
    )
    log_group.add_argument(
        "--progress-bar-width",
        type=positive_int,
        default=None,
        help="Fixed Rich progress-bar width. Default: auto from Rich display width.",
    )

    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.silva and not args.unite:
        die("At least one of --silva or --unite is required.")

    for input_path in (args.silva, args.unite):
        if input_path and not input_path.is_file():
            die(f"Input FASTA not found: {input_path}")

    if not (0 < args.display_width_ratio <= 1):
        die("--display-width-ratio must be > 0 and <= 1.")

    if args.batch_size < args.threads:
        args.batch_size = args.threads


def print_run_plan(args: argparse.Namespace, console, logger: logging.Logger) -> None:  # noqa: ANN001
    plan = {
        "SILVA": str(args.silva) if args.silva else "not provided",
        "UNITE": str(args.unite) if args.unite else "not provided",
        "Output directory": str(args.outdir),
        "Combined taxonomy": str(args.combined_taxonomy or (args.outdir / "ref_taxonomy.tsv")),
        "Threads": str(args.threads),
        "Batch size": str(args.batch_size),
        "Build indexes": "no" if args.skip_index else "yes",
        "Index jobs": str(args.index_jobs),
        "Keep organelles": "yes" if args.keep_organelles else "no",
        "Organelle keywords": ", ".join(args.organelle_keywords),
        "UNITE SH suffix": "no" if args.no_unite_sh_suffix else "yes",
        "UNITE SH ranks": ", ".join(args.unite_sh_suffix_ranks),
        "Display width": str(args.display_width) if args.display_width else f"terminal × {args.display_width_ratio:g}",
    }

    if console is not None:
        table = Table.grid(padding=(0, 2))
        table.add_column(style="bold cyan")
        table.add_column()
        for key, value in plan.items():
            table.add_row(key, value)
        console.print(Panel(table, title=PROGRAM, border_style="cyan"))
    else:
        logger.info("Run plan")
        for key, value in plan.items():
            logger.info("  %-18s: %s", key, value)


def print_summary(stats_list: list[ProcessStats], index_paths: list[Path], console, logger: logging.Logger) -> None:  # noqa: ANN001
    if console is not None:
        table = Table(title="Reference build summary", show_lines=False)
        table.add_column("Source", style="bold cyan")
        table.add_column("Seen", justify="right")
        table.add_column("Written", justify="right")
        table.add_column("16S arc/bac", justify="right")
        table.add_column("18S euk", justify="right")
        table.add_column("ITS", justify="right")
        table.add_column("Organelle skipped", justify="right")
        table.add_column("Domain skipped", justify="right")
        table.add_column("Time", justify="right")
        for st in stats_list:
            table.add_row(
                st.source,
                f"{st.total_seen:,}",
                f"{st.written:,}",
                f"{st.arc_bac:,}",
                f"{st.euk:,}",
                f"{st.unite:,}",
                f"{st.organelle_skipped:,}",
                f"{st.skipped_domain:,}",
                f"{st.elapsed_sec:.1f}s",
            )
        console.print(table)
        if index_paths:
            index_table = Table(title="minimap2 indexes")
            index_table.add_column("Index path", style="green")
            for path in index_paths:
                index_table.add_row(str(path))
            console.print(index_table)
    else:
        for st in stats_list:
            logger.info(
                "%s summary: seen=%s written=%s arc_bac=%s euk=%s unite=%s organelle_skipped=%s domain_skipped=%s time=%.1fs",
                st.source,
                f"{st.total_seen:,}",
                f"{st.written:,}",
                f"{st.arc_bac:,}",
                f"{st.euk:,}",
                f"{st.unite:,}",
                f"{st.organelle_skipped:,}",
                f"{st.skipped_domain:,}",
                st.elapsed_sec,
            )
        for path in index_paths:
            logger.info("Index: %s", path)


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    validate_args(args)
    preflight_output_conflicts(args)

    args.outdir.mkdir(parents=True, exist_ok=True)
    logger = setup_logging(args.outdir, args.log_file, args.quiet, args.no_rich)
    console = make_console(
        args.no_rich or args.quiet,
        display_width=args.display_width,
        display_width_ratio=args.display_width_ratio,
    )
    print_run_plan(args, console, logger)

    combined_taxonomy = args.combined_taxonomy or (args.outdir / "ref_taxonomy.tsv")
    id_registry = IdRegistry()
    all_records: list[RefTaxonomy] = []
    index_fastas: list[Path] = []
    stats_list: list[ProcessStats] = []

    with create_progress(console, logger, args.log_every, args.progress_bar_width) as progress:
        if args.silva:
            records, fastas, stats = process_silva(
                args.silva,
                args.outdir,
                id_registry,
                args.force,
                keep_organelles=args.keep_organelles,
                organelle_keywords=args.organelle_keywords,
                threads=args.threads,
                batch_size=args.batch_size,
                progress=progress,
                count_records=args.count_records,
                logger=logger,
            )
            all_records.extend(records)
            index_fastas.extend(fastas)
            stats_list.append(stats)

        if args.unite:
            records, fastas, stats = process_unite(
                args.unite,
                args.outdir,
                id_registry,
                default_domain=args.unite_default_domain,
                add_sh_suffix=not args.no_unite_sh_suffix,
                sh_suffix_ranks=args.unite_sh_suffix_ranks,
                force=args.force,
                threads=args.threads,
                batch_size=args.batch_size,
                progress=progress,
                count_records=args.count_records,
                logger=logger,
            )
            all_records.extend(records)
            index_fastas.extend(fastas)
            stats_list.append(stats)

        if not args.skip_combined:
            write_combined_taxonomy(combined_taxonomy, all_records, args.force, logger)

        built_indexes: list[Path] = []
        if not args.skip_index:
            built_indexes = build_minimap2_indexes(
                index_fastas,
                args.force,
                args.index_jobs,
                args.minimap2,
                progress,
                logger,
            )

    if id_registry.duplicate_count:
        logger.warning(
            "Duplicate short IDs renamed with numeric suffixes: %s",
            f"{id_registry.duplicate_count:,}",
        )

    print_summary(stats_list, built_indexes, console, logger)
    logger.info("Done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
