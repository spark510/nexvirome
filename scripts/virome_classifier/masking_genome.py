#!/usr/bin/env python3
"""
Compare a RefSeq FASTA against CHM13 and UniVec references and summarise
overlapping regions as BED files.

Usage example:
    python scripts/postprocessing/prepare_vdb.py \
        --refseq-file <repo>/resources/db/ncbi/refseq/viral.1.1.genomic.fna \
        --chm13-file /home/share/bowtie2_db/chm13v2.0/chm13v2.0.fa \
        --univec-file <repo>/resources/db/ncbi/univec/UniVec
"""

from __future__ import annotations

import argparse
import collections
import re
import subprocess
import sys
from pathlib import Path
from typing import Dict, Iterable, Tuple

import pandas as pd
from Bio import SeqIO
from Bio.Blast import NCBIXML

def build_blast_db(fasta: Path, dbtype: str = "nucl", force: bool = False) -> Path:
    """
    FASTA를 BLAST DB로 생성/확인. `force=False`면 기존 .nin/.nhr/.nsq 존재 시 건너뜀.
    반환값은 BLAST -db 인자로 쓸 프리픽스(Path).
    """
    if not fasta.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta}")

    prefix = fasta  # -out 프리픽스
    expected = [Path(f"{prefix}.{ext}") for ext in ("nhr", "nin", "nsq")]
    if not force and all(p.exists() for p in expected):
        print(f"Using existing BLAST database: {prefix}", file=sys.stderr)
        return prefix

    print(f"Building BLAST database: {prefix}", file=sys.stderr)
    cmd = [
        "makeblastdb",
        "-in", str(fasta),
        "-dbtype", dbtype,
        "-parse_seqids",
        "-out", str(prefix),
    ]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"makeblastdb failed for {fasta}\n"
            f"stdout: {e.stdout}\n"
            f"stderr: {e.stderr}"
        ) from e
    return prefix


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run blastn between refseq/chm13/univec and export BED summaries.",
    )
    parser.add_argument("--refseq-file", required=True, help="RefSeq FASTA to screen.")
    parser.add_argument("--chm13-file", required=True, help="CHM13 FASTA or blast DB.")
    parser.add_argument("--univec-file", required=True, help="UniVec FASTA or blast DB.")
    parser.add_argument(
        "--output-dir",
        default="results_vdb",
        help="Directory for XML and BED outputs (default: results_vdb).",
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=150,
        help="Minimum HSP length to keep (default: 150).",
    )
    parser.add_argument(
        "--min-identity",
        type=float,
        default=90.0,
        help="Minimum percent identity to keep (default: 90).",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=16,
        help="Threads passed to blastn (default: 16).",
    )
    parser.add_argument(
        "--task",
        default="megablast",
        help="blastn task (default: megablast).",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-run blast even if XML already exists.",
    )
    parser.add_argument(
        "--force-db",
        action="store_true",
        help="Rebuild BLAST databases even if index files already exist.",
    )
    parser.add_argument(
        "--skip-partial-cds",
        action="store_true",
        help="Do not emit BED for sequences annotated as partial CDS in RefSeq.",
    )
    return parser.parse_args()


def run_blast(query: Path, db: Path, out_xml: Path, *, threads: int, task: str, force: bool) -> Path:
    if out_xml.exists() and not force:
        print(f"Using existing BLAST results: {out_xml}", file=sys.stderr)
        return out_xml

    print(f"Running BLAST: {query.name} vs {db.name} -> {out_xml.name}", file=sys.stderr)
    out_xml.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "blastn",
        "-query",
        str(query),
        "-db",
        str(db),
        "-out",
        str(out_xml),
        "-outfmt",
        "5",
        "-task",
        task,
        "-num_threads",
        str(threads),
    ]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"blastn failed: {query} vs {db}\n"
            f"stdout: {e.stdout}\n"
            f"stderr: {e.stderr}"
        ) from e
    return out_xml


def parse_xml(xml_file: Path, *, min_len: int, min_iden: float, filter_self_hits: bool = False):
    query_data = collections.defaultdict(lambda: {"length": 0, "hsps": []})
    hit_data = collections.defaultdict(lambda: {"length": 0, "hsps": []})

    # 빈 XML이면 바로 반환 (BLAST hit 없음 등)
    if xml_file.stat().st_size == 0:
        return query_data, hit_data

    with xml_file.open() as handle:
        blast_records = NCBIXML.parse(handle)
        for record in blast_records:
            query_name = record.query
            query_data[query_name]["length"] = record.query_length

            for alignment in record.alignments:
                hit_name = alignment.title
                hit_data[hit_name]["length"] = alignment.length

                for hsp in alignment.hsps:
                    # Filter self-hits if requested
                    if filter_self_hits and query_name == hit_name:
                        continue

                    identity_perc = (hsp.identities / hsp.align_length) * 100
                    if hsp.align_length < min_len or identity_perc < min_iden:
                        continue

                    q_start, q_end = sorted((hsp.query_start, hsp.query_end))
                    query_data[query_name]["hsps"].append((q_start, q_end, hit_name))

                    h_start, h_end = sorted((hsp.sbjct_start, hsp.sbjct_end))
                    hit_data[hit_name]["hsps"].append((h_start, h_end, query_name))

    return query_data, hit_data


def parse_name(name: str) -> Tuple[str, str]:
    # Try to match common NCBI accession formats
    # Examples: NC_123456.1, NM_123456.1, AC_123456.1, gnl|BL_ORD_ID|123
    patterns = [
        r"([A-Z]{1,3}_[0-9]+\.[0-9]+)",  # Standard RefSeq: NC_123456.1
        r"([A-Z]{1,4}[0-9]+\.[0-9]+)",   # GenBank: AB123456.1
        r"(gnl\|[^|]+\|[^\s]+)",         # gnl format: gnl|BL_ORD_ID|123
    ]

    for pattern in patterns:
        matches = re.findall(pattern, name)
        if matches:
            acc = matches[0] if len(matches) == 1 else matches[-1]
            idx = name.find(acc)
            desc = name[idx + len(acc):].strip(" |,")
            return acc, desc

    # Fallback: split on first whitespace
    parts = name.split(maxsplit=1)
    if len(parts) == 2:
        return parts[0].strip(), parts[1].strip()

    print(f"[WARN] parse_name failed: {name}", file=sys.stderr)
    return name.strip(), ""


def write_intervals(
    fh,
    accver: str,
    desc: str,
    seq_length: int,
    intervals: Iterable[Tuple[int, int]],
    other_accver: str,
    other_desc: str,
):
    """
    Write merged intervals with depth information.
    Uses event-based algorithm instead of array for memory efficiency.
    """
    if not intervals:
        return

    # Create events: (position, +1 for start, -1 for end)
    events = []
    for s, e in intervals:
        events.append((s, 1))
        events.append((e, -1))

    # Sort events by position
    events.sort()

    # Sweep through events to find merged regions and depths
    current_depth = 0
    region_start = None
    depth_positions = []  # (position, depth)

    for pos, delta in events:
        if current_depth == 0 and delta == 1:
            # Start of new region
            region_start = pos
            depth_positions = [(pos, 1)]
        elif current_depth > 0:
            depth_positions.append((pos, current_depth + delta))

        current_depth += delta

        if current_depth == 0:
            # End of region - calculate average depth
            region_end = pos
            total_depth = 0
            prev_pos = region_start

            for i, (p, d) in enumerate(depth_positions[1:], 1):
                segment_len = p - prev_pos
                segment_depth = depth_positions[i - 1][1]
                total_depth += segment_len * segment_depth
                prev_pos = p

            region_len = region_end - region_start
            avg_depth = total_depth / region_len if region_len > 0 else 0

            fh.write(
                f"{accver}\t{region_start}\t{region_end}\t{avg_depth:.2f}\t{desc}\t{other_accver}\t{other_desc}\n"
            )
            depth_positions = []


def write_bed(data: Dict[str, dict], out_file: Path, *, mode: str):
    with out_file.open("w") as fh:
        if mode == "query":
            _write_bed_query_mode(data, fh)
        elif mode == "hit":
            _write_bed_hit_mode(data, fh)
        else:
            raise ValueError("mode must be 'query' or 'hit'")


def _write_bed_query_mode(data: Dict[str, dict], fh):
    for query_name, info in data.items():
        length = info["length"]
        hsps = info["hsps"]
        query_accver, query_desc = parse_name(query_name)

        hit_intervals = {}
        for start, end, hit_name in hsps:
            # Convert from 1-based inclusive (BLAST) to 0-based half-open (BED)
            s, e = start - 1, end
            hit_intervals.setdefault(hit_name, []).append((s, e))

        for hit_name, intervals in hit_intervals.items():
            hit_accver, hit_desc = parse_name(hit_name)
            write_intervals(fh, query_accver, query_desc, length, intervals, hit_accver, hit_desc)


def _write_bed_hit_mode(data: Dict[str, dict], fh):
    for hit_name, info in data.items():
        length = info["length"]
        hsps = info["hsps"]
        hit_accver, hit_desc = parse_name(hit_name)

        query_intervals = {}
        for start, end, query_name in hsps:
            # Convert from 1-based inclusive (BLAST) to 0-based half-open (BED)
            s, e = start - 1, end
            query_intervals.setdefault(query_name, []).append((s, e))

        for query_name, intervals in query_intervals.items():
            query_accver, query_desc = parse_name(query_name)
            write_intervals(fh, hit_accver, hit_desc, length, intervals, query_accver, query_desc)


def merge_bed(in_file: Path) -> Path:
    out_file = in_file.parent / f"{in_file.stem}_merged{in_file.suffix}"

    # Check if input file is empty or doesn't exist
    if not in_file.exists() or in_file.stat().st_size == 0:
        # Create empty output file
        out_file.touch()
        return out_file

    df = pd.read_csv(
        in_file,
        sep="\t",
        header=None,
        names=["acc", "start", "end", "depth", "desc", "other_acc", "other_desc"],
    )

    # Handle empty dataframe
    if df.empty:
        out_file.touch()
        return out_file

    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)

    merged_rows = []
    for (acc, desc), group in df.groupby(["acc", "desc"]):
        g = group.sort_values("start")

        merged = []
        for _, row in g.iterrows():
            s, e = row["start"], row["end"]
            if not merged or s > merged[-1][1]:
                # No overlap, start new interval
                merged.append([s, e, [str(row["depth"])], [row["other_acc"]], [row["other_desc"]]])
            else:
                # Overlapping or adjacent, merge
                merged[-1][1] = max(merged[-1][1], e)
                merged[-1][2].append(str(row["depth"]))
                merged[-1][3].append(row["other_acc"])
                merged[-1][4].append(row["other_desc"])

        for s, e, depths, other_accs, other_descs in merged:
            merged_rows.append(
                [
                    acc,
                    s,
                    e,
                    ",".join(depths),
                    desc,
                    ",".join(map(str, other_accs)),
                    ",".join(str(x) for x in other_descs if pd.notna(x)),
                ]
            )

    out_df = pd.DataFrame(
        merged_rows,
        columns=["acc", "start", "end", "depth", "desc", "other_accs", "other_descs"],
    ).sort_values(by=["acc", "start"])  # Sort by accession and start position (BED standard)

    out_df.to_csv(out_file, sep="\t", index=False, header=False)
    return out_file


def fasta_partial_cds_to_bed(fasta_file: Path, bed_file: Path):
    with bed_file.open("w") as out:
        for record in SeqIO.parse(fasta_file, "fasta"):
            desc = record.description
            if "partial cds" in desc.lower():
                accver, desc = _parse_fasta_name(desc)
                start, end = 0, len(record.seq)
                avg_depth = 0.0
                other_accver, other_desc = ".", "."
                out.write(
                    f"{accver}\t{start}\t{end}\t{avg_depth:.2f}\t{desc}\t{other_accver}\t{other_desc}\n"
                )


def _parse_fasta_name(name: str) -> Tuple[str, str]:
    parts = name.split(maxsplit=1)
    accver = parts[0].lstrip(">")
    desc = parts[1] if len(parts) > 1 else ""
    return accver, desc


def generate_all(targetseq: Path, chm13: Path, univec: Path, args: argparse.Namespace):
    # Validate input files exist
    targetseq_name = targetseq.name
    for fasta_file, name in [(targetseq, targetseq_name), (chm13, "CHM13"), (univec, "UniVec")]:
        if not fasta_file.exists():
            raise FileNotFoundError(f"{name} file not found: {fasta_file}")
        if fasta_file.stat().st_size == 0:
            raise ValueError(f"{name} file is empty: {fasta_file}")

    out_dir = Path(args.output_dir).resolve()  # 절대경로로 변환
    out_dir.mkdir(parents=True, exist_ok=True)

    targetseq_db = build_blast_db(targetseq, dbtype="nucl", force=args.force_db)
    chm13_db = build_blast_db(chm13, dbtype="nucl", force=args.force_db)
    univec_db = build_blast_db(univec, dbtype="nucl", force=args.force_db)



    # chm13_to_refseq_xml = out_dir / "chm13_to_refseq.xml"
    targetseq_to_chm13_xml = out_dir / f"{targetseq_name}_to_chm13.xml"
    targetseq_to_univec_xml = out_dir / f"{targetseq_name}_to_univec.xml"
    
    univec_to_targetseq_xml = out_dir / f"univec_to_{targetseq_name}.xml"
    chm13_to_targetseq_xml = out_dir / f"chm13_to_{targetseq_name}.xml"

    # run_blast(chm13_db, refseq_db, chm13_to_refseq_xml, threads=args.threads, task=args.task, force=args.force)
    run_blast(targetseq_db, chm13_db, targetseq_to_chm13_xml, threads=args.threads, task=args.task, force=args.force)
    run_blast(targetseq_db, univec_db, targetseq_to_univec_xml, threads=args.threads, task=args.task, force=args.force)
    run_blast(univec_db, targetseq_db, univec_to_targetseq_xml, threads=args.threads, task=args.task, force=args.force)    
    run_blast(chm13_db, targetseq_db, chm13_to_targetseq_xml, threads=args.threads, task=args.task, force=args.force)

    datasets = {
        # "refseq_to_chm13": (refseq_to_chm13_xml, "query", False),
        "targetseq_to_chm13": (targetseq_to_chm13_xml, "query", False),
        "targetseq_to_univec": (targetseq_to_univec_xml, "query", False),
        "univec_to_targetseq": (univec_to_targetseq_xml, "hit", False),
        "chm13_to_targetseq": (chm13_to_targetseq_xml, "hit", False),  # Filter self-hits
    }

    bed_paths = []
    for name, (xml_path, mode, filter_self) in datasets.items():
        query_data, hit_data = parse_xml(
            xml_path,
            min_len=args.min_length,
            min_iden=args.min_identity,
            filter_self_hits=filter_self
        )
        data = query_data if mode == "query" else hit_data

        bed_file = out_dir / f"{name}.{mode}.bed"
        write_bed(data, bed_file, mode=mode)
        bed_paths.append(bed_file)

        merged_file = merge_bed(bed_file)
        bed_paths.append(merged_file)

    # if not args.skip_partial_cds:
    #     partial_bed = out_dir / "refseq_partial_cds.bed"
    #     fasta_partial_cds_to_bed(refseq, partial_bed)
    #     bed_paths.append(partial_bed)

    return bed_paths


def main():
    args = parse_args()
    refseq = Path(args.refseq_file)
    chm13 = Path(args.chm13_file)
    univec = Path(args.univec_file)

    bed_paths = generate_all(refseq, chm13, univec, args)
    for p in bed_paths:
        print(p)


if __name__ == "__main__":
    main()
