#!/usr/bin/env python3
"""
Join the CoverM abundance table with the HMM annotation of every catalogue gene.

Gene IDs:
  abundance / rep_seq / cluster.tsv : SAMPLE_contig_12_3   (sample prefix added when pooling)
  hmm_hits.tsv                      : sample=SAMPLE, orf_id=contig_12_3
  -> key = sample + "_" + orf_id

Output: one row per catalogue gene (MMseqs2 representative):
  gene_id, origin_sample, contig_id, profile, evalue, bitscore,
  cluster_size, member_samples, member_profiles, protein_length,
  <all CoverM columns>, protein_sequence
"""
import argparse
import csv
import sys
from collections import defaultdict


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--abundance", required=True, help="CoverM table (<batch>_abundance.tsv)")
    ap.add_argument("--hits",      required=True, help="all_samples_hmm_hits.tsv (or several per-sample hmm_hits.tsv)", nargs="+")
    ap.add_argument("--clusters",  help="MMseqs2 <batch>_cluster.tsv (optional: adds cluster size and member samples/profiles)")
    ap.add_argument("--proteins",  help="MMseqs2 <batch>_rep_seq.fasta (adds length and sequence)")
    ap.add_argument("--output",    required=True)
    a = ap.parse_args()

    # 1. HMM annotation per pooled ID
    hits = {}
    for path in a.hits:
        with open(path) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row.get("sample") == "sample":          # repeated header if files were concatenated
                    continue
                gene = f"{row['sample']}_{row['orf_id']}"
                hits[gene] = row

    # 2. cluster members of every representative
    members = defaultdict(list)
    if a.clusters:
        with open(a.clusters) as fh:
            for line in fh:
                rep, mem = line.rstrip("\n").split("\t")[:2]
                members[rep].append(mem)

    proteins = read_fasta(a.proteins) if a.proteins else {}

    # 3. walk the abundance table, add annotation in front
    missing = 0
    with open(a.abundance) as fh, open(a.output, "w", newline="") as out:
        reader = csv.reader(fh, delimiter="\t")
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        header = next(reader)
        writer.writerow(
            ["gene_id", "origin_sample", "contig_id", "profile", "evalue", "bitscore",
             "cluster_size", "member_samples", "member_profiles", "protein_length"]
            + header[1:] + (["protein_sequence"] if proteins else [])
        )
        for row in reader:
            gene = row[0]
            h = hits.get(gene)
            if h is None:
                missing += 1
                h = {}
            mem = members.get(gene, [gene])
            mem_samples  = sorted({hits[m]["sample"]  for m in mem if m in hits})
            mem_profiles = sorted({hits[m]["profile"] for m in mem if m in hits})
            cluster_cols = [len(mem), ",".join(mem_samples), ",".join(mem_profiles)] if a.clusters else ["", "", ""]
            seq = proteins.get(gene, "")
            writer.writerow(
                [gene,
                 h.get("sample", ""),
                 f"{h['sample']}_{h['contig_id']}" if h else "",
                 h.get("profile", ""),
                 h.get("evalue", ""),
                 h.get("bitscore", ""),
                 *cluster_cols,
                 len(seq.rstrip("*")) if seq else ""]
                + row[1:] + ([seq] if proteins else [])
            )

    if missing:
        print(f"WARNING: {missing} catalogue genes had no HMM annotation (check IDs)", file=sys.stderr)


if __name__ == "__main__":
    main()
