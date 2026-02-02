#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import gzip
import os
import sys
from typing import Dict, List, Optional, Tuple
import pandas as pd

FASTA_EXTS = [".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz"]

def open_maybe_gzip(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

def reverse_complement(seq: str) -> str:
    comp = str.maketrans("ACGTacgtnN", "TGCAtgcanN")
    return seq.translate(comp)[::-1]

def read_fasta_to_dict(path: str) -> Dict[str, str]:
    seqs: Dict[str, List[str]] = {}
    cur: Optional[str] = None
    with open_maybe_gzip(path) as fh:
        for ln in fh:
            if not ln:
                continue
            if ln.startswith(">"):
                token = ln[1:].strip().split()[0]
                cur = token
                if cur not in seqs:
                    seqs[cur] = []
            else:
                if cur is not None:
                    seqs[cur].append(ln.strip())
    return {k: "".join(v) for k, v in seqs.items()}

def sanitize_contig_name(name: str) -> Tuple[str, str]:
    n = str(name).strip()
    if n.endswith(".0") and n.replace(".", "", 1).isdigit():
        n = n.split(".", 1)[0]
    n2 = n
    if n2.lower().startswith("chr"):
        n2 = n2[3:]
    n2 = n2.lstrip("0")
    return n, n2

def match_contig(seqs: Dict[str, str], chrom: str) -> Optional[str]:
    if chrom in seqs:
        return chrom
    orig, stripped = sanitize_contig_name(chrom)
    for k in seqs.keys():
        if k.lower() == orig.lower():
            return k
    for k in seqs.keys():
        if k.lower() in (f"chr{stripped}".lower(), stripped.lower()):
            return k
    base = orig.split()[0].split(".")[0]
    for k in seqs.keys():
        if k.split()[0].split(".")[0].lower() == base.lower():
            return k
    return None

def find_by_basename(root_dir: str, basename: str) -> List[str]:
    matches: List[str] = []
    for dirpath, _, filenames in os.walk(root_dir):
        for fn in filenames:
            if fn == basename:
                matches.append(os.path.join(dirpath, fn))
    return matches

def twin_name(basename: str) -> Optional[str]:
    if basename.endswith(".gz"):
        return basename[:-3]
    else:
        return basename + ".gz"

def choose_best_path(paths: List[str]) -> Optional[str]:
    if not paths:
        return None
    def score(p: str) -> int:
        name = os.path.basename(p).lower()
        s = 0
        if "genome" in name or "assembly" in name:
            s += 2
        if name.endswith((".fa", ".fa.gz", ".fna", ".fna.gz")):
            s += 1
        s -= p.count(os.sep) // 10
        return s
    return sorted(paths, key=lambda p: (-score(p), len(p)))[0]

def parse_args():
    ap = argparse.ArgumentParser(description="Extract genomic sequences using genome names from the coords table.")
    ap.add_argument("--coords", required=True, help="CSV produced by extract_gene_coords.py")
    ap.add_argument("--genome-dir", required=True, help="Root directory containing genome FASTA files (recursive).")
    ap.add_argument("--out-fasta", required=True, help="Output FASTA filepath for extracted sequences.")
    ap.add_argument("--out-table", required=True, help="Output CSV filepath for result table (+status).")
    ap.add_argument("--id-col", default="ID", help="Column for the unique ID to use in FASTA headers (default: ID).")
    ap.add_argument("--genome-col", default="Genome File Name", help="Column with genome filename (default: 'Genome File Name').")
    ap.add_argument("--chrom-col", default="chrom", help="Column for chromosome/seqid (default: chrom).")
    ap.add_argument("--start-col", default="start", help="Column for start coordinate (1-based, inclusive).")
    ap.add_argument("--end-col", default="end", help="Column for end coordinate (1-based, inclusive).")
    ap.add_argument("--strand-col", default="strand", help="Column for strand (+/-).")
    ap.add_argument("--seqname-col", default="seqname", help="Column for seqname label to report back (default: seqname).")
    ap.add_argument("--fasta-header", choices=["id", "full"], default="id",
                    help="FASTA header style: 'id' -> >ID ; 'full' -> >ID|chrom:used_start-used_end(strand)|seqname=...")
    ap.add_argument("--upstream", type=int, default=0, help="Bases upstream of the gene to include (strand-aware).")
    ap.add_argument("--downstream", type=int, default=0, help="Bases downstream of the gene to include (strand-aware).")
    return ap.parse_args()

def main():
    args = parse_args()

    df = pd.read_csv(args.coords)
    required = [args.id_col, args.genome_col, args.chrom_col, args.start_col, args.end_col, args.strand_col]
    for c in required:
        if c not in df.columns:
            sys.stderr.write(f"ERROR: required column '{c}' not found in {args.coords}\n")
            sys.exit(2)

    upstream = max(0, int(args.upstream))
    downstream = max(0, int(args.downstream))

    def to_int_safe(x):
        try:
            return int(float(x))
        except Exception:
            return None

    out_rows = []
    fasta_out = []

    df["_genome_name"] = df[args.genome_col].astype(str).apply(lambda s: os.path.basename(s))

    for gname, sub in df.groupby("_genome_name", dropna=True):
        paths = find_by_basename(args.genome_dir, gname)
        if not paths:
            twin = twin_name(gname)
            if twin:
                paths = find_by_basename(args.genome_dir, twin)
        genome_path = choose_best_path(paths)
        if genome_path is None:
            sys.stderr.write(f"WARN: No genome FASTA found matching basename '{gname}' under {args.genome_dir}\n")
            genome_seqs = {}
        else:
            genome_seqs = read_fasta_to_dict(genome_path)

        for _, row in sub.iterrows():
            rid = str(row[args.id_col])
            chrom_raw = row[args.chrom_col]
            start_raw = row[args.start_col]
            end_raw = row[args.end_col]
            strand = str(row[args.strand_col]) if pd.notna(row[args.strand_col]) else ""
            seqname = str(row[args.seqname_col]) if args.seqname_col in row else ""

            status = "ok"
            error = ""
            seq = ""
            chrom = str(chrom_raw) if pd.notna(chrom_raw) else ""
            start = to_int_safe(start_raw)
            end = to_int_safe(end_raw)

            used_start = None
            used_end = None
            contig_key = ""

            if genome_path is None:
                status = "no_genome_found"
            elif not chrom or start is None or end is None or not strand or strand not in ["+", "-"]:
                status = "missing_fields"
                error = f"chrom/start/end/strand invalid: chrom={chrom_raw}, start={start_raw}, end={end_raw}, strand={strand}"
            else:
                contig_key = match_contig(genome_seqs, chrom)
                if contig_key is None:
                    status = "contig_not_found"
                    error = f"chrom '{chrom}' not found in FASTA headers of {os.path.basename(genome_path) if genome_path else 'N/A'}"
                else:
                    gseq = genome_seqs[contig_key]
                    if strand == "+":
                        win_start = start - upstream
                        win_end = end + downstream
                    else:
                        win_start = start - downstream
                        win_end = end + upstream
                    used_start = max(1, win_start)
                    used_end = min(len(gseq), win_end)
                    if used_start > used_end:
                        status = "invalid_range"
                        error = f"window invalid after clamp: win_start={win_start}, win_end={win_end}, len={len(gseq)}"
                    else:
                        subseq = gseq[used_start-1:used_end]
                        if strand == "-":
                            subseq = reverse_complement(subseq)
                        seq = subseq

            if status == "ok" and seq:
                if args.fasta_header == "id":
                    header = f">{rid}"
                else:
                    header = f">{rid}|{chrom}:{used_start}-{used_end}({strand})|seqname={seqname}"
                fasta_out.append((header, seq))

            out_row = dict(row)
            out_row.update({
                "genome_path": genome_path if genome_path else "",
                "genome_file": os.path.basename(genome_path) if genome_path else "",
                "contig_used": contig_key,
                "used_start": used_start if used_start is not None else "",
                "used_end": used_end if used_end is not None else "",
                "flank_upstream": upstream,
                "flank_downstream": downstream,
                "extracted_len": len(seq) if seq else 0,
                "status": status,
                "error": error,
            })
            out_rows.append(out_row)

    with open(args.out_fasta, "w") as fh:
        for header, seq in fasta_out:
            fh.write(header + "\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i+60] + "\n")

    out_df = pd.DataFrame(out_rows)
    out_df.to_csv(args.out_table, index=False)

    print(f"Wrote {len(fasta_out)} sequences to {args.out_fasta}", file=sys.stderr)
    print(f"Wrote {len(out_df)} rows to {args.out_table}", file=sys.stderr)

if __name__ == "__main__":
    main()
