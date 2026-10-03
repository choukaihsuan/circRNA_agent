import csv,glob,json,os
H='/home3/choukaihsuan/'
B=H+'GSE113230_results/benchmark/'
print(open(B+'compute_cost.tsv').read())
print(open(B+'pr_curve.tsv').read()[:3000])
for p in ['GSE55872_results/circRNA/SRR444655/find_circ/splice_sites.bed']:
    L=[l for l in open(H+p)]; print(p,len(L), sum(1 for l in L if 'CIRCULAR' in l))
