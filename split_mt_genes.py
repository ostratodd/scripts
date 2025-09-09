#!/usr/bin/env python3
"""
Download cephalopod mitochondrial genomes by accession and split into per-gene FASTA files.

Usage:
  python split_mt_genes.py accessions.txt --email you@univ.edu [--api-key YOUR_NCBI_KEY] -o outdir [--translate]

Notes:
- 'accessions.txt' should contain one nucleotide accession per line (e.g., NC_012345, AY123456, OQ123456).
- By default, outputs nucleotide sequences. Use --translate to also write amino-acid FASTAs for CDS genes.
- Designed to be tolerant to COX1/COI naming, cytochrome b/cytb/cob, NADH/NAD names, and common tRNA/rRNA labels.
"""

import argparse
import os
import re
import sys
import time
from collections import defaultdict
from pathlib import Path

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# --- Gene name normalization -------------------------------------------------

# Map many-to-one common mitochondrial gene name variants
GENE_SYNONYMS = {
    # Cytochrome oxidase subunits
    "cox1": {"cox1", "coxi", "co1", "coi", "coxi subunit i", "cytochrome c oxidase subunit i"},
    "cox2": {"cox2", "co2", "coii", "cytochrome c oxidase subunit ii"},
    "cox3": {"cox3", "co3", "coiii", "cytochrome c oxidase subunit iii"},
    # Cytochrome b
    "cob": {"cob", "cob", "cytb", "cytochrome b"},
    # ATP synthase
    "atp6": {"atp6", "atpase6", "atp synthase f0 subunit 6"},
    "atp8": {"atp8", "atpase8", "atp synthase f0 subunit 8"},
    # NADH dehydrogenase (Complex I)
    "nad1": {"nad1", "nd1", "nadh dehydrogenase subunit 1"},
    "nad2": {"nad2", "nd2", "nadh dehydrogenase subunit 2"},
    "nad3": {"nad3", "nd3", "nadh dehydrogenase subunit 3"},
    "nad4": {"nad4", "nd4", "nadh dehydrogenase subunit 4"},
    "nad4l": {"nad4l", "nd4l", "nadh dehydrogenase subunit 4l"},
    "nad5": {"nad5", "nd5", "nadh dehydrogenase subunit 5"},
    "nad6": {"nad6", "nd6", "nadh dehydrogenase subunit 6"},
    # rRNA (12S/16S)
    "rrnS": {"rrns", "12s", "12s rrna", "ssu rrna", "small subunit rrna", "rrn s"},
    "rrnL": {"rrnl", "16s", "16s rrna", "lsu rrna", "large subunit rrna", "rrn l"},
}

# Invert synonyms for quick lookup
LOWER_TO_CANON = {}
for canon, variants in GENE_SYNONYMS.items():
    for v in variants:
        LOWER_TO_CANON[v.lower()] = canon

# Helper to normalize gene names
def normalize_gene_name(raw: str, feature_type: str) -> str:
    """
    Normalize gene names to stable keys:
    - protein-coding/rRNAs via synonyms map
    - tRNAs as 'trnX' optionally with anticodon, e.g. 'trnL-uur'
    """
    if not raw:
        raw = ""

    raw_clean = re.sub(r"\s+", " ", raw.strip()).lower()

    # Try direct map first (covers most protein-coding + rRNA)
    if raw_clean in LOWER_TO_CANON:
        return LOWER_TO_CANON[raw_clean]

    # Handle common spellings like "cytochrome c oxidase subunit i" → cox1
    # (fall back already done via LOWER_TO_CANON)

    # tRNA handling
    if feature_type == "tRNA" or "trna" in raw_clean or "tRNA" in raw:
        # Try to extract amino acid letter or name + anticodon
        # Examples of product strings:
        #   "tRNA-Leu (UUR)", "tRNA-Lys", "transfer RNA-Val (UAC)"
        aa = None
        anticodon = None

        # Extract AA long/short
        m = re.search(r"trna[-\s]*(\w+)", raw_clean)
        if m:
            aa = m.group(1)  # could be 'leu', 'lys', 'val', or single-letter
        # Extract anticodon in parentheses
        m2 = re.search(r"\(([^)]+)\)", raw_clean)
        if m2:
            anticodon = m2.group(1).replace(" ", "").lower()

        # Map 3-letter/word AA to one-letter, default to first letter if unknown
        AA_TO_ONE = {
            "ala":"a","arg":"r","asn":"n","asp":"d","cys":"c","gln":"q","glu":"e","gly":"g",
            "his":"h","ile":"i","leu":"l","lys":"k","met":"m","phe":"f","pro":"p","ser":"s",
            "thr":"t","trp":"w","tyr":"y","val":"v",
        }
        # Accept already one-letter
        aa_one = None
        if aa:
            aa3 = aa[:3]
            if aa in AA_TO_ONE.values() and len(aa) == 1:
                aa_one = aa
            elif aa3 in AA_TO_ONE:
                aa_one = AA_TO_ONE[aa3]
            else:
                aa_one = aa[0]  # best-effort

        if aa_one:
            if anticodon:
                # Normalize common ones like UUR/AGN to lowercase
                anticodon_norm = anticodon.replace("uur","uur").replace("agn","agn")
                return f"trn{aa_one.upper()}-{anticodon_norm}"
            return f"trn{aa_one.upper()}"

        # If we can't parse, just return 'trn' + sanitized raw
        return "trn_" + re.sub(r"[^a-z0-9]+", "_", raw_clean).strip("_")

    # Final fallback: sanitize the raw name
    fallback = re.sub(r"[^a-z0-9]+", "_", raw_clean).strip("_")
    return fallback or (feature_type.lower())
# -----------------------------------------------------------------------------


def parse_args():
    ap = argparse.ArgumentParser(description="Split cephalopod mitochondrial genomes into per-gene FASTAs.")
    ap.add_argument("accessions_file", help="Text file with one NCBI nucleotide accession per line.")
    ap.add_argument("-o", "--outdir", default="mtgenes_out", help="Output directory (default: mtgenes_out)")
    ap.add_argument("--email", required=True, help="Email for NCBI Entrez (required).")
    ap.add_argument("--api-key", default=None, help="NCBI API key (optional but recommended).")
    ap.add_argument("--translate", action="store_true",
                    help="Also write amino-acid FASTAs for CDS genes.")
    ap.add_argument("--delay", type=float, default=0.35,
                    help="Polite delay between NCBI requests in seconds (default 0.35).")
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


def fetch_genbank(acc, rettype="gb", retmode="text"):
    from Bio import Entrez
    try:
        with Entrez.efetch(db="nucleotide", id=acc, rettype=rettype, retmode=retmode) as handle:
            return handle.read()
    except Exception as e:
        sys.stderr.write(f"[WARN] Failed to fetch {acc}: {e}\n")
        return None


def record_species_name(rec):
    # Prefer full organism name; fall back to source feature if needed
    org = rec.annotations.get("organism")
    if org:
        return org
    for feat in rec.features:
        if feat.type == "source":
            orgq = feat.qualifiers.get("organism", [])
            if orgq:
                return orgq[0]
    return "Unknown_organism"


def best_feature_name(feat):
    """
    Return a raw gene name string from a feature using common qualifiers.
    """
    q = feat.qualifiers
    for key in ("gene", "locus_tag", "product", "note"):
        vals = q.get(key)
        if vals:
            return vals[0]
    return ""


def extract_features_by_gene(rec):
    """
    Extract sequences for CDS/tRNA/rRNA and group by normalized gene name.
    Returns dict: gene_name -> list of (SeqRecord, feature_type, is_cds_bool)
    """
    gene_bins = defaultdict(list)
    full_seq = rec.seq
    organism = record_species_name(rec)
    acc = rec.id

    for feat in rec.features:
        if feat.type not in {"CDS", "tRNA", "rRNA"}:
            continue

        # Extract nucleotide sequence respecting strand and joins
        try:
            subseq = feat.location.extract(full_seq)
        except Exception as e:
            sys.stderr.write(f"[WARN] {acc}: failed to extract {feat.type} ({e})\n")
            continue

        raw_name = best_feature_name(feat)
        gene_name = normalize_gene_name(raw_name, feat.type)

        # Build a clean header: Species | accession | gene
        species_slim = re.sub(r"\s+", "_", organism).replace("/", "_")
        header_id = f"{species_slim}|{acc}|{gene_name}"

        # Nucleotide record
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
        # Nucleotide
        nfa = nuc_dir / f"{gene}.fasta"
        with open(nfa, "a") as fh:
            for nrec, ftype, is_cds in entries:
                SeqIO.write(nrec, fh, "fasta")

        # Protein (optional; CDS only, translate with standard mito invertebrate table = 5)
        if prot_dir and any(is_cds for (_, _, is_cds) in entries):
            pfa = prot_dir / f"{gene}.faa"
            with open(pfa, "a") as fh:
                for nrec, ftype, is_cds in entries:
                    if not is_cds:
                        continue
                    try:
                        # Translate using invertebrate mitochondrial code (NCBI table 5)
                        # Stop codons trimmed at first '*' for cleanliness
                        prot = Seq(nrec.seq).translate(table=5, to_stop=True)
                        prec = SeqRecord(prot, id=nrec.id, description="")
                        SeqIO.write(prec, fh, "fasta")
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

    # Aggregate across all accessions
    all_per_gene = defaultdict(list)

    for i, acc in enumerate(accs, 1):
        sys.stderr.write(f"[INFO] ({i}/{len(accs)}) Fetching {acc}...\n")
        gbtxt = fetch_genbank(acc, rettype="gbwithparts", retmode="text")
        if gbtxt is None:
            continue

        # Parse as GenBank record
        from io import StringIO
        try:
            rec = next(SeqIO.parse(StringIO(gbtxt), "genbank"))
        except Exception as e:
            sys.stderr.write(f"[WARN] Failed to parse GenBank for {acc}: {e}\n")
            continue

        # Skip if not mitochondrial? (Optional check)
        # We assume these are mitochondrial; but you can enforce via source/organelle if desired.

        bins = extract_features_by_gene(rec)
        # Merge into global bins
        for gene, items in bins.items():
            all_per_gene[gene].extend(items)

        time.sleep(args.delay)

    # Write outputs
    write_fastas(all_per_gene, outdir, write_proteins=args.translate)

    # Write a manifest of genes and counts
    manifest = outdir / "manifest.tsv"
    with open(manifest, "w") as fh:
        fh.write("gene\tfeature_types\tcount\n")
        for gene, items in sorted(all_per_gene.items()):
            types = sorted(set(t for _, t, _ in items))
            fh.write(f"{gene}\t{','.join(types)}\t{len(items)}\n")

    sys.stderr.write(f"[DONE] Wrote per-gene FASTAs to: {outdir}\n")
    sys.stderr.write("Subfolders: nucleotide/ (always), protein/ (if --translate)\n")
    sys.stderr.write(f"Manifest: {manifest}\n")


if __name__ == "__main__":
    main()

