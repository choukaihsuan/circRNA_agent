import csv,math,collections
H='/home3/choukaihsuan/'
def rd(p): return list(csv.DictReader(open(p),delimiter='\t'))
bm=rd(H+'GSE133998_results/de/biomarker_candidates.tsv')
f=float
S=[min(-math.log10(f(r['pvalue'])),10) for r in bm]
F=[min(abs(f(r['log2FC'])),5) for r in bm]
C=[f(r['confidence_score']) for r in bm]
M=[f(r['n_mirna']) for r in bm]; R=[f(r['n_rbp']) for r in bm]
mm=lambda v,x:(x-min(v))/(max(v)-min(v)) if max(v)>min(v) else 0
print('Table2 recompute: rank gene S F C K M R score(file)')
for i,r in enumerate(bm[:10]):
    K=f(r['in_circbase'])
    comp=[mm(S,S[i]),mm(F,F[i]),mm(C,C[i]),K,mm(M,M[i]),mm(R,R[i])]
    print(i+1,r['circbase_gene'],' '.join('%.2f'%c for c in comp),'mean=%.4f'%(sum(comp)/6),r['biomarker_score'])
print('n bm',len(bm))
