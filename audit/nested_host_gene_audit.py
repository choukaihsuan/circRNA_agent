"""Audit of host-gene assignment: is the assigned gene nested inside a larger gene?

assign_isoforms.py assigns each circRNA the SMALLEST gene that contains it. This audit flags
an assignment as NESTED when another GTF gene (a) contains the circRNA and (b) fully contains
the assigned gene, i.e. the assigned gene sits inside a larger host. Such a case means the
smallest-gene rule chose the embedded element over the gene that actually hosts it.

Inputs (server): <GSE>_results/circRNA/isoform_groups.tsv for each dataset, and the GENCODE
v19 GTF (genes.gtf). Gene features only are read (type == gene).

usage: python3 nested_host_gene_audit.py <gtf> <out_tsv> <GSE> [<GSE> ...]
(results directories are taken as ~/<GSE>_results)
"""
import collections
import os
import re
import sys

import numpy as np
import pandas as pd

GTF, OUT = sys.argv[1], sys.argv[2]
GSES = sys.argv[3:]

ATTR_GID = re.compile(r'gene_id "([^"]+)"')
ATTR_GNAME = re.compile(r'gene_name "([^"]+)"')
ATTR_GTYPE = re.compile(r'gene_type "([^"]+)"')

genes = []  # (gene_id, gene_name, gene_type, chrom, start, end)
with open(GTF) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t", 8)
        if len(f) < 9 or f[2] != "gene":
            continue
        a = f[8]
        gid = ATTR_GID.search(a); gn = ATTR_GNAME.search(a); gt = ATTR_GTYPE.search(a)
        genes.append((gid.group(1) if gid else "", gn.group(1) if gn else "",
                      gt.group(1) if gt else "", f[0], int(f[3]), int(f[4])))
G = pd.DataFrame(genes, columns=["gene_id", "gene_name", "gene_type", "chrom", "start", "end"])
print(f"[audit] GTF gene features: {len(G)}", file=sys.stderr)

by_chr = {}
for c, sub in G.groupby("chrom"):
    sub = sub.reset_index(drop=True)
    by_chr[c] = dict(start=sub.start.values, end=sub.end.values, idx=sub.index.values, df=sub)
gid_to_row = {}
for i, r in G.iterrows():
    gid_to_row.setdefault(r.gene_id, i)  # first occurrence; gene_id is unique in GENCODE

rows = []
for g in GSES:
    path = os.path.expanduser(f"~/{g}_results/circRNA/isoform_groups.tsv")
    iso = pd.read_csv(path, sep="\t", dtype=str)
    n_total = len(iso)
    n_intergenic = int((iso.gene_name == "intergenic").sum())
    for r in iso.itertuples(index=False):
        if r.gene_name == "intergenic" or not isinstance(r.gene_name, str):
            rows.append(dict(dataset=g, circ_id=r.circ_id, assigned_gene=r.gene_name, assigned_type="",
                             nested=False, containing_larger=""))
            continue
        cs, ce = int(r.start), int(r.end)
        chrom = r.chrom
        assigned_row = gid_to_row.get(r.gene_id)
        if assigned_row is None:
            # gene_id not found in GTF; fall back to gene_name, containing this circRNA
            cand = G[(G.gene_name == r.gene_name) & (G.chrom == chrom)]
            cand = cand[(cand.start <= cs) & (cand.end >= ce)]
            if cand.empty:
                rows.append(dict(dataset=g, circ_id=r.circ_id, assigned_gene=r.gene_name, assigned_type="",
                                 nested=False, containing_larger="UNRESOLVED(gene_id not in GTF)"))
                continue
            assigned_row = cand.index[0]
        a = G.loc[assigned_row]
        blk = by_chr.get(chrom)
        if blk is None:
            rows.append(dict(dataset=g, circ_id=r.circ_id, assigned_gene=a.gene_name, assigned_type=a.gene_type,
                             nested=False, containing_larger=""))
            continue
        # genes that contain the circRNA span
        m = (blk["start"] <= cs) & (blk["end"] >= ce)
        cand_idx = blk["idx"][m]
        bigger = []
        for k in cand_idx:
            row = blk["df"].loc[k]
            if row.gene_id == a.gene_id:
                continue
            # strictly larger AND fully contains the assigned gene
            if row.start <= a.start and a.end <= row.end and (row.end - row.start) > (a.end - a.start):
                bigger.append(f"{row.gene_name}({row.gene_type})")
        rows.append(dict(dataset=g, circ_id=r.circ_id, assigned_gene=a.gene_name, assigned_type=a.gene_type,
                         nested=bool(bigger), containing_larger=";".join(bigger)))
    print(f"[audit] {g}: {n_total} circRNAs, {n_intergenic} intergenic", file=sys.stderr)

out = pd.DataFrame(rows)
out.to_csv(OUT, sep="\t", index=False)

# summary
assigned = out[out.assigned_gene.notna() & (out.assigned_gene != "intergenic")]
summ = []
for g, sub in assigned.groupby("dataset"):
    nested = sub[sub.nested]
    summ.append(dict(dataset=g, genic_circ=len(sub), nested=len(nested),
                     nested_pct=round(100 * len(nested) / len(sub), 2) if len(sub) else 0.0,
                     nested_small_types=int(nested.assigned_type.isin(
                         ["snRNA", "snoRNA", "miRNA", "misc_RNA", "scaRNA", "rRNA", "scRNA"]).sum()),
                     nested_pseudogene=int(nested.assigned_type.str.contains("pseudogene", na=False).sum())))
S = pd.DataFrame(summ)
S.to_csv(OUT.replace(".tsv", "_summary.tsv"), sep="\t", index=False)
print(S.to_string(index=False))
ut = out.drop_duplicates(["circ_id"])
print("unique circ_id across datasets:", len(ut), "nested:", int(ut.nested.sum()))
