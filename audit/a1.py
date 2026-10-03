import csv,collections,os
H='/home3/choukaihsuan/'
def rd(p): return list(csv.DictReader(open(p),delimiter='\t'))
R=H+'GSE133998_results/'
cm=open(R+'circRNA/count_matrix.tsv').readline().split('\t'); n=sum(1 for _ in open(R+'circRNA/count_matrix.tsv'))-1
print('GSE133998 count_matrix rows',n,'cols',len(cm)-1)
iso=rd(R+'de/isoform_switching.tsv'); print('iso cols',list(iso[0].keys())); print('iso rows',len(iso))
genes=set(r['gene_name'] for r in iso); print('genes',len(genes))
def f(x):
    try: return float(x)
    except: return None
for thr in (0.05,0.1):
    s=[r for r in iso if f(r.get('padj_within_gene')) is not None and f(r['padj_within_gene'])<thr and abs(f(r['delta_iui']) or 0)>0.1]
    print('FDR<',thr,'rows',len(s),'genes',len(set(r['gene_name'] for r in s)))
print('is_switching col count', sum(1 for r in iso if r.get('is_switching') in ('TRUE','True','true','1')))
# multi-isoform genes
gc=collections.Counter(r['gene_name'] for r in iso)
print('genes with >=2 rows',sum(1 for g,c in gc.items() if c>=2))
print(iso[0])
# venn
def sigset(path,pc='pvalue'):
    S=set()
    for r in rd(path):
        p=f(r[pc]);l=f(r['log2FC'])
        if p is not None and l is not None and p<0.05 and abs(l)>1: S.add(r['circ_id'])
    return S
E_=sigset(R+'de/de_results_edgeR_ciriquant.tsv');D=sigset(R+'de/de_results_deseq2.tsv');L=sigset(R+'de/de_results_limma.tsv')
print('edgeR',len(E_),'deseq',len(D),'limma',len(L),'E&D',len(E_&D),'E&L',len(E_&L),'D&L',len(D&L),'all3',len(E_&D&L))
print(open(R+'de/de_results_deseq2.tsv').readline())
bm=rd(R+'de/biomarker_candidates.tsv'); print(list(bm[0].keys()))
for r in bm[:10]: print({k:r[k] for k in r if k in('circ_id','gene','gene_name','biomarker_score','log2FC','n_mirna','n_rbp','sig_norm','fc_norm','conf_norm','mirna_norm','rbp_norm','known_bonus','circbase_id','circbase_gene','Type')})
