#!/usr/bin/env python3
"""Build short-id marker references and taxonomy tables from SILVA/UNITE FASTA.

The output taxonomy schema is the one consumed by meta_marker_count.sh:

ref_id  marker  domain  phylum  class  order  family  genus  species
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, TextIO


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

MISSING = "Unclassified"


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


class IdRegistry:
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


def log(message: str) -> None:
    print(f"[build_ref_files] {message}", file=sys.stderr)


def open_text(path: Path) -> TextIO:
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


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


def has_silva_organelle_lineage(header: str) -> bool:
    return any(part.strip() in {"Mitochondria", "Chloroplast"} for part in taxonomy_parts(header))


def domain_key(domain: str) -> str:
    return domain.lower().replace("_", "").replace(" ", "")


def write_tax_header(handle: TextIO) -> None:
    handle.write("\t".join(TAX_HEADER) + "\n")


def write_tax_row(handle: TextIO, record: RefTaxonomy) -> None:
    handle.write("\t".join(record.row()) + "\n")


def process_unite(
    fasta: Path,
    outdir: Path,
    id_registry: IdRegistry,
    default_domain: str,
    force: bool,
) -> list[RefTaxonomy]:
    base = strip_fasta_suffix(fasta)
    fasta_out = outdir / f"{base}.shortid.fasta"
    tax_out = outdir / f"{base}.taxonomy.tsv"
    ensure_can_write(fasta_out, force)
    ensure_can_write(tax_out, force)

    records: list[RefTaxonomy] = []
    count = 0
    with fasta_out.open("wt", encoding="utf-8") as fasta_handle, tax_out.open(
        "wt", encoding="utf-8"
    ) as tax_handle:
        write_tax_header(tax_handle)
        for header, seq in read_fasta(fasta):
            ref_id = id_registry.unique(short_id_from_header(header))
            tax = parse_taxonomy(header, default_domain=default_domain)
            record = RefTaxonomy(ref_id=ref_id, marker="ITS", taxonomy=tax)
            write_fasta_record(fasta_handle, ref_id, seq)
            write_tax_row(tax_handle, record)
            records.append(record)
            count += 1

    log(f"UNITE records written: {count}")
    log(f"  FASTA: {fasta_out}")
    log(f"  taxonomy: {tax_out}")
    return records


def process_silva(
    fasta: Path,
    outdir: Path,
    id_registry: IdRegistry,
    force: bool,
    keep_organelles: bool,
) -> list[RefTaxonomy]:
    base = strip_fasta_suffix(fasta)
    arc_bac_fasta = outdir / f"{base}.dna.arc_bac.shortid.fasta"
    arc_bac_tax = outdir / f"{base}.dna.arc_bac.taxonomy.tsv"
    euk_fasta = outdir / f"{base}.dna.euk.shortid.fasta"
    euk_tax = outdir / f"{base}.dna.euk.taxonomy.tsv"
    for path in (arc_bac_fasta, arc_bac_tax, euk_fasta, euk_tax):
        ensure_can_write(path, force)

    records: list[RefTaxonomy] = []
    counts = {"arc_bac": 0, "euk": 0, "organelle": 0, "skipped": 0}
    with arc_bac_fasta.open("wt", encoding="utf-8") as ab_fa, arc_bac_tax.open(
        "wt", encoding="utf-8"
    ) as ab_tax, euk_fasta.open("wt", encoding="utf-8") as euk_fa, euk_tax.open(
        "wt", encoding="utf-8"
    ) as eu_tax:
        write_tax_header(ab_tax)
        write_tax_header(eu_tax)
        for header, seq in read_fasta(fasta):
            if not keep_organelles and has_silva_organelle_lineage(header):
                counts["organelle"] += 1
                continue
            tax = parse_silva_taxonomy(header)
            dkey = domain_key(tax.domain)
            if dkey in {"bacteria", "archaea"}:
                ref_id = id_registry.unique(short_id_from_header(header))
                record = RefTaxonomy(ref_id=ref_id, marker="16S", taxonomy=tax)
                write_fasta_record(ab_fa, ref_id, seq)
                write_tax_row(ab_tax, record)
                records.append(record)
                counts["arc_bac"] += 1
            elif dkey == "eukaryota":
                ref_id = id_registry.unique(short_id_from_header(header))
                record = RefTaxonomy(ref_id=ref_id, marker="18S", taxonomy=tax)
                write_fasta_record(euk_fa, ref_id, seq)
                write_tax_row(eu_tax, record)
                records.append(record)
                counts["euk"] += 1
            else:
                counts["skipped"] += 1

    log(f"SILVA arc/bac records written: {counts['arc_bac']}")
    log(f"  FASTA: {arc_bac_fasta}")
    log(f"  taxonomy: {arc_bac_tax}")
    log(f"SILVA euk records written: {counts['euk']}")
    log(f"  FASTA: {euk_fasta}")
    log(f"  taxonomy: {euk_tax}")
    if counts["skipped"]:
        log(f"SILVA records skipped because domain was not Bacteria/Archaea/Eukaryota: {counts['skipped']}")
    if counts["organelle"]:
        log(f"SILVA organelle records skipped by exact lineage item Mitochondria/Chloroplast: {counts['organelle']}")
    return records


def ensure_can_write(path: Path, force: bool) -> None:
    if path.exists() and not force:
        raise SystemExit(f"Refusing to overwrite existing file without --force: {path}")


def write_combined_taxonomy(path: Path, records: Iterable[RefTaxonomy], force: bool) -> int:
    ensure_can_write(path, force)
    count = 0
    with path.open("wt", encoding="utf-8") as handle:
        write_tax_header(handle)
        for record in records:
            write_tax_row(handle, record)
            count += 1
    log(f"Combined taxonomy records written: {count}")
    log(f"  taxonomy: {path}")
    return count


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate short-id FASTA and taxonomy TSV files for meta_marker_count.sh."
    )
    parser.add_argument("--silva", type=Path, help="SILVA FASTA or FASTA.GZ file.")
    parser.add_argument("--unite", type=Path, help="UNITE FASTA or FASTA.GZ file.")
    parser.add_argument(
        "-o",
        "--outdir",
        type=Path,
        default=Path("ref"),
        help="Output directory. Default: ref",
    )
    parser.add_argument(
        "--combined-taxonomy",
        type=Path,
        help="Combined taxonomy TSV path. Default: OUTDIR/ref_taxonomy.tsv",
    )
    parser.add_argument(
        "--unite-default-domain",
        default="Fungi",
        help="Domain used for UNITE records whose header has no parseable domain. Default: Fungi",
    )
    parser.add_argument(
        "--skip-combined",
        action="store_true",
        help="Only write per-database FASTA/taxonomy files.",
    )
    parser.add_argument(
        "--keep-organelles",
        action="store_true",
        help="Keep SILVA records with exact lineage item Mitochondria or Chloroplast. Default skips them.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing output files.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if not args.silva and not args.unite:
        raise SystemExit("At least one of --silva or --unite is required.")

    for input_path in (args.silva, args.unite):
        if input_path and not input_path.is_file():
            raise SystemExit(f"Input FASTA not found: {input_path}")

    args.outdir.mkdir(parents=True, exist_ok=True)
    combined_taxonomy = args.combined_taxonomy or (args.outdir / "ref_taxonomy.tsv")
    id_registry = IdRegistry()
    all_records: list[RefTaxonomy] = []

    if args.silva:
        all_records.extend(
            process_silva(
                args.silva,
                args.outdir,
                id_registry,
                args.force,
                keep_organelles=args.keep_organelles,
            )
        )
    if args.unite:
        all_records.extend(
            process_unite(
                args.unite,
                args.outdir,
                id_registry,
                default_domain=args.unite_default_domain,
                force=args.force,
            )
        )

    if not args.skip_combined:
        write_combined_taxonomy(combined_taxonomy, all_records, args.force)

    if id_registry.duplicate_count:
        log(f"Duplicate short IDs renamed with numeric suffixes: {id_registry.duplicate_count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
