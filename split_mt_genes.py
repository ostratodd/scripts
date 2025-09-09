#!/usr/bin/env python3
"""
Download cephalopod mitochondrial genomes by accession and split into per-gene FASTA files.

Default: extracts CDS and rRNA (skips tRNAs).
- Normalize many gene-name spellings (COX/CO/ND/NAD, CYTB/COB, 12S/16S variants).
- Translate CDS to protein with --translate (mito invertebrate code 5).
- Report any unmapped raw labels with --report-unknowns.

Usage:
  python split_mt_genes.py accessions.txt --email you@uni.edu -o outdir [--translate] [--include-trna] [--report-unknowns]
"""

import argparse
import re
import sys
import time
from collections import defaultdict, Counter
from pathlib import Path
from io import StringIO

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# -------------------------- Synonyms (curated) -------------------------------

GENE_SYNONYMS = {
    # Cytochrome c oxidase
    "cox1": {
        "cox1","coxi","co1","coi","cytochrome c oxidase subunit i",
        "cytochrome oxidase subunit i","cytochrome c oxidase subunit 1",
    },
    "cox2": {
        "cox2","co2","coii","cytochrome c oxidase subunit ii",
        "cytochrome oxidase subunit ii","cytochrome c oxidase subunit 2",
    },
    "cox3": {
        "cox3","co3","coiii","cytochrome c oxidase subunit iii",
        "cytochrome oxidase subunit iii","cytochrome c oxidase subunit 3",
    },

    # Cytochrome b
    "cob": {"cob","cytb","cytochrome b","cytochrome-b"},

    # ATP synthase
    "atp6": {"atp6","atpase6","atp synthase f0 subunit 6","atp synthase subunit 6"},
    "atp8": {"atp8","atpase8","atp synthase f0 subunit 8","atp synthase subunit 8"},

    # NADH dehydrogenase (Complex I)
    "nad1": {"nad1","nd1","nadh dehydrogenase subunit 1"},
    "nad2": {"nad2","nd2","nadh dehydrogenase subunit 2"},
    "nad3": {"nad3","nd3","nadh dehydrogenase subunit 3"},
    "nad4": {"nad4","nd4","nadh dehydrogenase subunit 4"},
    "nad4l": {"nad4l","nd4l","nadh dehydrogenase subunit 4l"},
    "nad5": {"nad5","nd5","nadh dehydrogenase subunit 5"},
    "nad6": {"nad6","nd6","nadh dehydrogenase subunit 6"},

    # rRNAs (12S / 16S) — include many variants & your observed oddballs
    "rrnS": {
        "rrns","12s","12s rrna","ssu rrna","small subunit rrna",
        "12s ribosomal rna","rrn s","s-rrna","12 s rrna",
        "small subunit ribosomal rna",
        # oddball labels seen in your outputs (mapped to rrnS)
        "aaw62 gr01","al384 gr01","apq60 gr01","ask61 gr01","at075 gr01",
        "h876 mgr01","keg13 r01","yc02 gr01",
    },
    "rrnL": {
        "rrnl","16s","16s rrna","lsu rrna","large subunit rrna",
        "16s ribosomal rna","rrn l","l-rrna","16 s rrna",
        "large subunit ribosomal rna",
        # oddball labels seen in your outputs (mapped to rrnL)
        "aaw62 gr02","al384 gr02","apq60 gr02","ask61 gr02","at075 gr02",
        "h876 mgr02","keg13 r02","yc02 gr02",
    },
}

# Build lowercase lookup
LOWER_TO_CANON = {}
for canon, variants in GENE_SYNONYMS.items():
    for v in variants:
        LOWER_TO_CANON[v.lower()] = canon

# ----------------------------- CLI / Helpers ---------------------------------

def parse_args():
    ap = argparse.ArgumentParser(description="Split cephalopod mitogenomes into per-gene FASTAs.")
    ap.add_argument("accessions_file", help="Text file with one NCBI nucleotide accession per line.")
    ap.add_argument("-o", "--outdir", default="mtgenes_out", help="Output directory (default: mtgenes_out)")
    ap.add_argument("--email", required=True, help="Email for NCBI Entrez (required).")
    ap.add_argument("--api-key", default=None, help="NCBI API key (optional).")
    ap.add_argument("--translate", action="store_true", help="Also write amino-acid FASTAs for CDS genes.")
    ap.add_argument("--delay", type=float, default=0.35, help="Delay between NCBI requests (s).")
    ap.add_argument("--include-trna", action="store_true", help="Include tRNA features (default: OFF).")
    ap.add_argument("--report-unknowns", action="store_true",
                    help="Write unmapped raw labels and counts to unknown_labels.tsv")
    return ap.parse_args()


def read_accessions(path):
    accs = []
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            accs.append(s)
    return accs


def fetch_genbank(acc, rettype="gbwithparts", retmode="text"):
    try:
        with Entrez.efetch(db="nucleotide", id=acc, rettype=rettype, retmode=retmode) as handle:
            return handle.read()
    except Exception as e:
        sys.stderr.write(f"[WARN] Failed to fetch {acc}: {e}\n")
        return None


def record_species_name(rec):
    org = rec.annotations.get("organism")
    if org:
        return org
    for feat in rec.features:
        if feat.type == "source":
            vals = feat.qualifiers.get("organism", [])
            if vals:
                return vals[0]
    return "Unknown_organism"


def best_feature_name(feat):
    # Prefer explicit gene names; fall back to product, then notes
    for key in ("gene", "locus_tag", "product", "note"):
        vals = feat.qualifiers.get(key)
        if vals:
            return vals[0]
    return ""


def normalize_gene_name(raw: str, feature_type: str, unknown_counter: Counter = None) -> str:
    """
    Normalize to stable keys:
      - lowercase, underscores->spaces, collapse spaces
      - map via LOWER_TO_CANON
      - fallback: sanitized raw label (for visibility)
    Unknown raw labels are counted if unknown_counter is provided.
    """
    raw = raw or ""
    raw_clean = raw.strip().lower().replace("_", " ")
    raw_clean = re.sub(r"\s+", " ", raw_clean)

    if raw_clean in LOWER_TO_CANON:
        return LOWER_TO_CANON[raw_clean]

    # Count unknowns for reporting
    if unknown_counter is not None:
        unknown_counter[raw_clean] += 1

    # Fallback sanitized
    fallback = re.sub(r"[^a-z0-9]+", "_", raw_clean).strip("_")
    return fallback or feature_type.lower()


def extract_features_by_gene(rec, include_trna=False, unknown_counter=None):
    """
    Extract sequences for CDS/rRNA (and optional tRNA), grouped by normalized gene name.
    Returns dict: gene_name -> list of (SeqRecord, feature_type, is_cds_bool)
    """
    allowed = {"CDS", "rRNA"} | ({"tRNA"} if include_trna else set())
    gene_bins = defaultdict(list)

    full_seq = rec.seq
    organism = record_species_name(rec)
    acc = rec.id

    for feat in rec.features:
        if feat.type not in allowed:
            continue

        try:
            subseq = feat.location.extract(full_seq)
        except Exception as e:
            sys.stderr.write(f"[WARN] {acc}: failed to extract {feat.type} ({e})\n")
            continue

        raw_name = best_feature_name(feat)
        gene_name = normalize_gene_name(raw_name, feat.type, unknown_counter=unknown_counter)

        species_slim = re.sub(r"\s+", "_", organism).replace("/", "_")
        header_id = f"{species_slim}|{acc}|{gene_name}"
        nrec = SeqRecord(Seq(str(subseq)), id=header_id, description="")
        gene_bins[gene_name].append((nrec, feat.type, feat.type == "CDS"))

    return gene_bins


def write_fastas(per_gene, outdir, write_proteins=False):
    outdir = Path(outdir)
    nuc_dir = outdir / "nucleotide"
    nuc_dir.mkdir(parents=True, exist_ok=True)

    prot_dir = outdir / "protein" if write_proteins else None
    if prot_dir:
        prot_dir.mkdir(parents=True, exist_ok=True)

    for gene, entries in sorted(per_gene.items()):
        # nucleotide
        nfa = nuc_dir / f"{gene}.fasta"
        with open(nfa, "a") as fh:
            for nrec, _, _ in entries:
                SeqIO.write(nrec, fh, "fasta")

        # proteins
        if prot_dir and any(is_cds for (_, _, is_cds) in entries):
            pfa = prot_dir / f"{gene}.faa"
            with open(pfa, "a") as fh:
                for nrec, _, is_cds in entries:
                    if not is_cds:
                        continue
                    try:
                        prot = Seq(nrec.seq).translate(table=5, to_stop=True)
                        SeqIO.write(SeqRecord(prot, id=nrec.id, description=""), fh, "fasta")
                    except Exception as e:
                        sys.stderr.write(f"[WARN] translate fail {nrec.id}: {e}\n")


def main():
    args = parse_args()
    Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key

    accs = read_accessions(args.accessions_file)
    if not accs:
        sys.stderr.write("No accessions found. Exiting.\n")
        sys.exit(1)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    all_per_gene = defaultdict(list)
    unknown_counter = Counter() if args.report_unknowns else None

    for i, acc in enumerate(accs, 1):
        sys.stderr.write(f"[INFO] ({i}/{len(accs)}) Fetching {acc}...\n")
        gbtxt = fetch_genbank(acc, rettype="gbwithparts", retmode="text")
        if gbtxt is None:
            continue

        try:
            rec = next(SeqIO.parse(StringIO(gbtxt), "genbank"))
        except Exception as e:
            sys.stderr.write(f"[WARN] Failed to parse GenBank for {acc}: {e}\n")
            continue

        bins = extract_features_by_gene(
            rec,
            include_trna=args.include_trna,
            unknown_counter=unknown_counter
        )
        for gene, items in bins.items():
            all_per_gene[gene].extend(items)

        time.sleep(args.delay)

    # Write outputs
    write_fastas(all_per_gene, outdir, write_proteins=args.translate)

    # Manifest
    manifest = outdir / "manifest.tsv"
    with open(manifest, "w") as fh:
        fh.write("gene\tfeature_types\tcount\n")
        for gene, items in sorted(all_per_gene.items()):
            types = sorted(set(t for _, t, _ in items))
            fh.write(f"{gene}\t{','.join(types)}\t{len(items)}\n")

    # Unknowns (optional)
    if unknown_counter:
        unk_path = outdir / "unknown_labels.tsv"
        with open(unk_path, "w") as fh:
            fh.write("raw_label_normalized\tcount\n")
            for label, cnt in unknown_counter.most_common():
                fh.write(f"{label}\t{cnt}\n")
        sys.stderr.write(f"[INFO] Wrote unknown label report: {unk_path}\n")

    sys.stderr.write(f"[DONE] Wrote per-gene FASTAs to: {outdir}\n")
    sys.stderr.write("Subfolders: nucleotide/ (always), protein/ (if --translate)\n")
    sys.stderr.write(f"Manifest: {manifest}\n")


if __name__ == "__main__":
    main()

