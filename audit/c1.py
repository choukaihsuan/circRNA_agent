import csv,json,glob,os,time
H='/home3/choukaihsuan/'
DS=['GSE113230','GSE58135','GSE133998','SRP156355','GSE77509','GSE130078','GSE248612','GSE221107','PRJNA553289','GSE229705','GSE148036','GSE121842','GSE136569','GSE143797','GSE108735','GSE171011','GSE97239','GSE192410','GSE192849']
tot=0
for d in DS:
    n=0;m=0
    for r in csv.DictReader(open(H+d+'_results/de/de_results.tsv'),delimiter='\t'):
        try: pa=float(r['padj']); l=float(r['log2FC'])
        except: continue
        if pa<0.05: m+=1
        if pa<0.05 and abs(l)>1: n+=1
    tot+=n; print(d,'padj<0.05:',m,'& |lfc|>1:',n)
print('total padj-sig',tot)
for r in csv.DictReader(open(H+'GSE113230_results/circRNA/isoform_groups.tsv'),delimiter='\t'):
    if r['circ_id'] in('chr1:117944808|117963271','chr1:117944808|117948267') or r['circ_id']=='chr12:46622936|46648719': print(r['circ_id'],r['gene_name'],r['strand'])
for l in open(H+'circRNA_agent/reference/circbase_hg19.txt'):
    if 'hsa_circ_0000211' in l or 'hsa_circ_0000118\t' in l: print(l.rstrip()[:160])
x=json.load(open(H+'GSE113230_results/de/interactions.json'))
k='chr12:46622936|46648719'
if k in x: 
    v=x[k]; print('SLC38A1 keys',list(v.keys())[:6])
    for kk in v:
        if 'mirna' in kk.lower() or 'rbp' in kk.lower(): print(kk,len(v[kk]) if hasattr(v[kk],'__len__') else v[kk])
print('GSE323364 interactions mtime',time.strftime('%F',time.localtime(os.path.getmtime(H+'GSE323364_results/de/interactions.json'))))
