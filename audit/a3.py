import csv,collections
H='/home3/choukaihsuan/'
DS=['GSE113230','GSE58135','GSE133998','SRP156355','GSE77509','GSE130078','GSE248612','GSE221107','PRJNA553289','GSE229705','GSE148036','GSE121842','GSE136569','GSE143797','GSE108735','GSE171011','GSE97239','GSE192410','GSE192849']
def rd(p): return list(csv.DictReader(open(p),delimiter='\t'))
sig={};dat={}
for d in DS:
    S={};
    for r in rd(H+d+'_results/de/de_results.tsv'):
        p=float(r['pvalue']);l=float(r['log2FC'])
        if p<0.05 and abs(l)>1: S[r['circ_id']]=(l,r['Type'])
    sig[d]=S
c=collections.defaultdict(list)
for d in DS:
    for k in sig[d]: c[k].append(d)
ann={}
for d in DS:
    try:
        for r in rd(H+d+'_results/circRNA/circbase_annotated.tsv'): ann.setdefault(r['circ_id'],(r.get('circbase_id'),r.get('circbase_gene')))
    except Exception as e: print('noann',d,e)
iso={}
for d in DS:
    try:
        for r in rd(H+d+'_results/circRNA/isoform_groups.tsv'): iso.setdefault(r['circ_id'],(r.get('gene_name'),r.get('region')))
    except Exception as e: print('noiso',d)
def show(k):
    print(k,ann.get(k),iso.get(k),len(c[k]),[(d,round(sig[d][k][0],2),sig[d][k][1]) for d in c[k]])
print('== >=5'); 
for k,v in c.items():
    if len(v)>=5: show(k)
n4=[k for k,v in c.items() if len(v)==4]; print('exactly4',len(n4))
print('== exactly 4 list')
for k in n4: print(k,ann.get(k),iso.get(k)[0] if k in iso else None,c[k])
print('== sensitivity: recurrence of the >=4 set excluding GSE192410')
for k,v in c.items():
    if len(v)>=4:
        w=[d for d in v if d!='GSE192410']
        if len(w)<3 or len(v)>=5: print(k,len(v),'->',len(w))
for k in ['chr13:33091994|33101669','chr12:46622936|46648719']: show(k)
for g in ['hsa_circ_0000118','hsa_circ_0006156','hsa_circ_0000211','hsa_circ_0008934','hsa_circ_0001772','hsa_circ_0001246','hsa_circ_0035796']:
    ks=[k for k,a in ann.items() if a[0]==g]
    print(g,ks,[ (k,len(c[k]),c[k]) for k in ks])
# MAN1A2 
print('MAN1A2 circs sig:',[(k,c[k]) for k in c if iso.get(k,('',))[0]=='MAN1A2'])
# top-ranked
for d in ['GSE108735','GSE136569']:
    b=rd(H+d+'_results/de/biomarker_candidates.tsv')[:3]
    for r in b: print(d,r['rank'],r['circ_id'],r['circbase_id'],r['circbase_gene'],r['log2FC'],r['Type'],r['biomarker_score'])
