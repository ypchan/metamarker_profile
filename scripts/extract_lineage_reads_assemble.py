#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Extract marker reads assigned to a target lineage and run an internal OLC assembly.

The selection logic mirrors metamarker_profile abundance assignment:
  1. parse mate-labeled PAF query IDs ending in /1 or /2;
  2. keep the best valid hit per pair_id + mate;
  3. arbitrate R1/R2 conflicts at pair level;
  4. extract reads assigned to the requested rank/lineage.

Outputs include mate-specific FASTQ, merged_segments.fq.gz, selected-read tables,
and Python-only OLC contigs assembled from merged_segments.fq.gz.
"""

from __future__ import annotations

import argparse
import ast
import csv
import gzip
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Iterator, List, Optional, Sequence, Tuple


RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
MARKER_ORDER = ["16S", "18S", "ITS"]
DEFAULT_MIN_IDENTITY = {"16S": 0.97, "18S": 0.97, "ITS": 0.95}
DEFAULT_MIN_ALN_LEN = {"16S": 80, "18S": 80, "ITS": 80}
DEFAULT_MIN_QCOV = {"16S": 0.60, "18S": 0.60, "ITS": 0.60}
DEFAULT_MIN_MAPQ = 0


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
    mate_lineage: str
    assigned_lineage: str
    assigned_target_id: str
    read_weight: int


@dataclass
class OlcContig:
    name: str
    seq: str
    qual: str
    support: int


def open_text_auto(path: Path, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode, encoding=None if "b" in mode else "utf-8")
    return open(path, mode, encoding=None if "b" in mode else "utf-8")


def open_gzip_text(path: Path, mode: str = "wt"):
    return gzip.open(path, mode, encoding="utf-8")


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
        required = {"ref_id", "marker", *RANKS}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Taxonomy table missing columns: {','.join(sorted(missing))}")
        for row in reader:
            ref_id = (row.get("ref_id") or "").strip()
            if not ref_id:
                continue
            rec = {k: (row.get(k) or "").strip() for k in ["marker", *RANKS]}
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
            if not lineage_match(mate_lineage):
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
                    mate_lineage=mate_lineage,
                    assigned_lineage=mate_lineage,
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
        if not lineage_match(assigned_lineage):
            continue
        read_weight = len(pair_hits)
        for hit in pair_hits:
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
                    mate_lineage=hit.taxonomy.get(rank, "unidentified"),
                    assigned_lineage=assigned_lineage,
                    assigned_target_id=winner.target_id,
                    read_weight=read_weight,
                )
            )
    return selected


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
                "mate_lineage",
                "assigned_lineage",
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
                    r.mate_lineage,
                    r.assigned_lineage,
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


def extract_fastq_records(
    outdir: Path,
    pipeline_outdir: Path,
    selected: Sequence[SelectedRead],
    orient: bool,
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
        print(f"WARNING: {missing} selected reads were not found in marker FASTQ files.", file=sys.stderr)
    return summary


def qscore(char: str) -> int:
    return max(0, ord(char) - 33)


def qchar(score: int) -> str:
    return chr(max(0, min(93, int(score))) + 33)


def mean_quality(qual: str) -> float:
    if not qual:
        return 0.0
    return sum(qscore(c) for c in qual) / float(len(qual))


def read_olc_segments(path: Path, min_len: int, max_reads: int) -> Tuple[List[OlcContig], int]:
    if not path.exists():
        raise FileNotFoundError(f"merged segments FASTQ not found: {path}")

    by_seq: Dict[str, OlcContig] = {}
    raw_reads = 0
    for header, seq, _plus, qual in iter_fastq(path):
        raw_reads += 1
        if max_reads > 0 and raw_reads > max_reads:
            print(f"WARNING: OLC read cap reached at {max_reads:,} reads; remaining records skipped.", file=sys.stderr)
            break
        seq = seq.upper()
        if len(seq) < min_len:
            continue
        if len(seq) != len(qual):
            raise ValueError(f"FASTQ sequence/quality length mismatch in {path}: {header}")
        existing = by_seq.get(seq)
        if existing is None:
            by_seq[seq] = OlcContig(name=normalize_read_id(header), seq=seq, qual=qual, support=1)
        else:
            existing.support += 1
            if mean_quality(qual) > mean_quality(existing.qual):
                existing.name = normalize_read_id(header)
                existing.qual = qual

    contigs = list(by_seq.values())
    contigs.sort(key=lambda c: (len(c.seq), c.support, c.name), reverse=True)
    return contigs, min(raw_reads, max_reads) if max_reads > 0 else raw_reads


def overlap_identity(left: str, right: str, overlap: int) -> float:
    matches = 0
    left_start = len(left) - overlap
    for i in range(overlap):
        if left[left_start + i] == right[i]:
            matches += 1
    return matches / float(overlap)


def best_pair_overlap(a: OlcContig, b: OlcContig, min_overlap: int, min_identity: float) -> Tuple[int, float]:
    max_overlap = min(len(a.seq), len(b.seq))
    for overlap in range(max_overlap, min_overlap - 1, -1):
        identity = overlap_identity(a.seq, b.seq, overlap)
        if identity >= min_identity:
            return overlap, identity
    return 0, 0.0


def find_best_overlap_full(
    contigs: Sequence[OlcContig],
    min_overlap: int,
    min_identity: float,
) -> Optional[Tuple[int, int, int, float]]:
    best: Optional[Tuple[Tuple[int, float, int], int, int, int, float]] = None
    for i, a in enumerate(contigs):
        if len(a.seq) < min_overlap:
            continue
        for j, b in enumerate(contigs):
            if i == j or len(b.seq) < min_overlap:
                continue
            overlap, identity = best_pair_overlap(a, b, min_overlap, min_identity)
            if overlap <= 0:
                continue
            score = (overlap, identity, a.support + b.support)
            if best is None or score > best[0]:
                best = (score, i, j, overlap, identity)
    if best is None:
        return None
    return best[1], best[2], best[3], best[4]


def find_best_overlap_seeded(
    contigs: Sequence[OlcContig],
    min_overlap: int,
    min_identity: float,
    seed_len: int,
) -> Optional[Tuple[int, int, int, float]]:
    seed_len = max(1, min(seed_len, min_overlap))
    prefix_index: Dict[str, List[int]] = {}
    for j, contig in enumerate(contigs):
        if len(contig.seq) >= min_overlap:
            prefix_index.setdefault(contig.seq[:seed_len], []).append(j)

    best: Optional[Tuple[Tuple[int, float, int], int, int, int, float]] = None
    max_target_len = max((len(c.seq) for c in contigs), default=0)
    for i, a in enumerate(contigs):
        max_overlap = min(len(a.seq), max_target_len)
        for overlap in range(max_overlap, min_overlap - 1, -1):
            seed_start = len(a.seq) - overlap
            seed = a.seq[seed_start : seed_start + seed_len]
            for j in prefix_index.get(seed, []):
                if i == j or len(contigs[j].seq) < overlap:
                    continue
                identity = overlap_identity(a.seq, contigs[j].seq, overlap)
                if identity < min_identity:
                    continue
                score = (overlap, identity, a.support + contigs[j].support)
                if best is None or score > best[0]:
                    best = (score, i, j, overlap, identity)
            if best is not None and best[0][0] == overlap:
                break
    if best is None:
        return None
    return best[1], best[2], best[3], best[4]


def consensus_base(left_base: str, left_q: str, right_base: str, right_q: str) -> Tuple[str, str]:
    if left_base == right_base:
        return left_base, qchar(max(qscore(left_q), qscore(right_q)))
    if left_base == "N":
        return right_base, right_q
    if right_base == "N":
        return left_base, left_q
    if qscore(left_q) > qscore(right_q):
        return left_base, left_q
    if qscore(right_q) > qscore(left_q):
        return right_base, right_q
    return "N", qchar(min(qscore(left_q), qscore(right_q)))


def merge_contigs(left: OlcContig, right: OlcContig, overlap: int, merge_no: int) -> OlcContig:
    left_prefix_len = len(left.seq) - overlap
    seq_parts = [left.seq[:left_prefix_len]]
    qual_parts = [left.qual[:left_prefix_len]]

    for i in range(overlap):
        base, qual = consensus_base(
            left.seq[left_prefix_len + i],
            left.qual[left_prefix_len + i],
            right.seq[i],
            right.qual[i],
        )
        seq_parts.append(base)
        qual_parts.append(qual)

    seq_parts.append(right.seq[overlap:])
    qual_parts.append(right.qual[overlap:])
    seq = "".join(seq_parts)
    qual = "".join(qual_parts)
    return OlcContig(name=f"olc_merge_{merge_no}", seq=seq, qual=qual, support=left.support + right.support)


def run_python_olc(
    segments_fastq: Path,
    outdir: Path,
    min_overlap: int,
    min_identity: float,
    min_read_len: int,
    min_contig_len: int,
    max_reads: int,
    full_scan_limit: int,
    seed_len: int,
) -> Dict[str, int | float | str]:
    if min_overlap < 1:
        raise ValueError("--olc-min-overlap must be positive.")
    if not (0 < min_identity <= 1):
        raise ValueError("--olc-min-identity must be in (0, 1].")

    contigs, input_reads = read_olc_segments(segments_fastq, min_read_len, max_reads)
    unique_segments = len(contigs)
    merges = 0

    while len(contigs) > 1:
        if len(contigs) <= full_scan_limit:
            best = find_best_overlap_full(contigs, min_overlap, min_identity)
        else:
            best = find_best_overlap_seeded(contigs, min_overlap, min_identity, seed_len)
        if best is None:
            break
        i, j, overlap, identity = best
        if i == j:
            break
        left = contigs[i]
        right = contigs[j]
        merges += 1
        merged = merge_contigs(left, right, overlap, merges)
        for idx in sorted([i, j], reverse=True):
            del contigs[idx]
        contigs.append(merged)
        if merges % 100 == 0:
            print(f"OLC merges={merges:,}; contigs={len(contigs):,}; last_overlap={overlap}; identity={identity:.4f}", file=sys.stderr)

    contigs.sort(key=lambda c: (len(c.seq), c.support, c.name), reverse=True)
    kept = [c for c in contigs if len(c.seq) >= min_contig_len]

    fasta_out = outdir / "olc_contigs.fa"
    fastq_out = outdir / "olc_contigs.fq.gz"
    summary_out = outdir / "olc_summary.tsv"
    with open(fasta_out, "w", encoding="utf-8") as fa, open_gzip_text(fastq_out, "wt") as fq:
        for idx, contig in enumerate(kept, start=1):
            name = f"contig_{idx} len={len(contig.seq)} support={contig.support}"
            fa.write(f">{name}\n{wrap_fasta(contig.seq)}\n")
            fq.write(f"@contig_{idx} support={contig.support}\n{contig.seq}\n+\n{contig.qual}\n")

    summary: Dict[str, int | float | str] = {
        "segments_fastq": str(segments_fastq),
        "input_reads": input_reads,
        "unique_segments": unique_segments,
        "merges": merges,
        "raw_contigs": len(contigs),
        "kept_contigs": len(kept),
        "min_overlap": min_overlap,
        "min_identity": min_identity,
        "min_read_len": min_read_len,
        "min_contig_len": min_contig_len,
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


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Extract reads assigned to a target lineage from metamarker_profile outputs "
            "and assemble merged_segments.fq.gz with a Python-only OLC assembler."
        )
    )
    p.add_argument(
        "--segments-fastq",
        help="Run OLC directly on an existing merged_segments.fq.gz and skip lineage extraction.",
    )
    p.add_argument("--outdir", default="metamarker_profile_out", help="metamarker_profile output directory.")
    p.add_argument("--taxonomy", help="Override taxonomy TSV. Default: path from OUTDIR/run_config.tsv.")
    p.add_argument("--rank", default="genus", choices=RANKS, help="Taxonomic rank to match. Default: genus.")
    p.add_argument("--lineage", help="Lineage value to extract at --rank.")
    p.add_argument("--markers", help="Comma-separated markers. Default: markers from run_config.tsv, otherwise all.")
    p.add_argument("--sample-id", help="Comma-separated sample IDs. Default: all samples in sample_manifest.tsv.")
    p.add_argument(
        "--mode",
        choices=["assigned", "mate-hit"],
        default="assigned",
        help=(
            "assigned matches the abundance table after pair arbitration; "
            "mate-hit keeps only mates whose own best hit matches the lineage. Default: assigned."
        ),
    )
    p.add_argument("--match", choices=["exact", "contains", "regex"], default="exact", help="Lineage match mode.")
    p.add_argument("--ignore-case", action="store_true", help="Case-insensitive lineage matching.")
    p.add_argument("--output-dir", help="Output directory. Default: OUTDIR/lineage_extract/<rank>_<lineage>.")
    p.add_argument(
        "--no-orient",
        action="store_true",
        help="Do not reverse-complement reads whose best PAF hit is on the minus strand.",
    )
    p.add_argument(
        "--assembler",
        choices=["olc", "none"],
        default="olc",
        help="Run the built-in Python OLC assembler on merged_segments.fq.gz. Default: olc.",
    )
    p.add_argument("--olc-min-overlap", type=int, default=40, help="Minimum suffix-prefix overlap. Default: 40.")
    p.add_argument("--olc-min-identity", type=float, default=0.98, help="Minimum overlap identity. Default: 0.98.")
    p.add_argument("--olc-min-read-len", type=int, default=40, help="Minimum input segment length. Default: 40.")
    p.add_argument("--olc-min-contig-len", type=int, default=0, help="Minimum output contig length. Default: 0.")
    p.add_argument("--olc-max-reads", type=int, default=50000, help="Maximum segments to assemble; 0 disables cap. Default: 50000.")
    p.add_argument(
        "--olc-full-scan-limit",
        type=int,
        default=500,
        help="Use exact all-vs-all overlap search up to this many unique segments. Default: 500.",
    )
    p.add_argument("--olc-seed-len", type=int, default=12, help="Seed length for large seeded overlap search. Default: 12.")
    return p


def main() -> int:
    args = build_parser().parse_args()
    if args.segments_fastq:
        segments_fastq = Path(args.segments_fastq)
        output_dir = Path(args.output_dir) if args.output_dir else segments_fastq.parent
        output_dir.mkdir(parents=True, exist_ok=True)
        if args.assembler != "none":
            summary = run_python_olc(
                segments_fastq,
                output_dir,
                args.olc_min_overlap,
                args.olc_min_identity,
                args.olc_min_read_len,
                args.olc_min_contig_len,
                args.olc_max_reads,
                args.olc_full_scan_limit,
                args.olc_seed_len,
            )
            print(
                f"OLC assembled {summary['kept_contigs']} contigs from "
                f"{summary['input_reads']} reads ({summary['unique_segments']} unique segments)."
            )
            print(f"Output directory: {output_dir}")
        else:
            print("No assembly requested because --assembler none was set.")
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

    selected: List[SelectedRead] = []
    for sample in samples:
        sample_id = sample["sample_id"]
        for marker in markers:
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
            selected.extend(select_reads(hits, args.rank, lineage_match, args.mode))

    write_selected_table(output_dir / "selected_reads.tsv", selected)
    if not selected:
        write_summary(output_dir / "summary.tsv", [])
        print(f"No reads matched lineage={args.lineage!r} rank={args.rank}; wrote empty tables to {output_dir}")
        return 0

    summary = extract_fastq_records(output_dir, pipeline_outdir, selected, orient=not args.no_orient)
    write_summary(output_dir / "summary.tsv", summary)

    segments_fastq = output_dir / "merged_segments.fq.gz"
    if args.assembler != "none":
        run_python_olc(
            segments_fastq,
            output_dir,
            args.olc_min_overlap,
            args.olc_min_identity,
            args.olc_min_read_len,
            args.olc_min_contig_len,
            args.olc_max_reads,
            args.olc_full_scan_limit,
            args.olc_seed_len,
        )

    total_reads = sum(int(row["selected_reads"]) for row in summary)
    total_pairs = sum(int(row["selected_pairs"]) for row in summary)
    print(f"Selected {total_reads} reads from {total_pairs} read IDs.")
    print(f"Output directory: {output_dir}")
    print("Key files: reads.R1.fq.gz, reads.R2.fq.gz, merged_segments.fq.gz, olc_contigs.fa, selected_reads.tsv")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except BrokenPipeError:
        raise SystemExit(0)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
