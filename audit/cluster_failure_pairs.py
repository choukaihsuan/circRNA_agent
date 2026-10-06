"""D-2: count pairs of output (consensus) loci within one sample whose distance max(|dstart|,|dend|) <= slop.
Correct clustering would leave none. Reads existing --dump-dist files only. usage: cluster_failure_pairs.py <dump_root> [slop=10] [dump_dirname=dump]"""
import glob, os, sys
import numpy as np, pandas as pd

root, slop = sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 10
dname = sys.argv[3] if len(sys.argv) > 3 else "dump"
rows, ex = [], []
for g in sorted(os.listdir(root)):
    files = [f for f in glob.glob(f"{root}/{g}/{dname}/*.tsv") if not f.endswith("_s.tsv")]
    if not files:
        continue
    n_loc = n_pair = n_samp = 0
    for f in sorted(files):
        d = pd.read_csv(f, sep="\t", na_values="NA")
        d = d[d.confidence_score.notna()]
        n_loc += len(d)
        for chrom, c in d.groupby("chrom"):
            c = c.sort_values("start")
            st, en = c.start.to_numpy(), c.end.to_numpy()
            for i in range(len(st)):
                j = i + 1
                while j < len(st) and st[j] - st[i] <= slop:
                    if abs(en[j] - en[i]) <= slop:
                        n_pair += 1
                        if len(ex) < 12:
                            ex.append((g, os.path.basename(f), chrom, st[i], en[i], st[j], en[j]))
                    j += 1
    rows.append((g, n_loc, n_pair, n_pair / n_loc))
    print(f"{g:12s} loci {n_loc:>8,}  pairs<= {slop}bp {n_pair:>6,}  share {n_pair / n_loc:.4%}", flush=True)
t = sum(r[1] for r in rows), sum(r[2] for r in rows)
print(f"TOTAL loci {t[0]:,} pairs {t[1]:,} share {t[1] / t[0]:.4%}")
for e in ex:
    print("example", *e)
