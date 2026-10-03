import csv
H='/home3/choukaihsuan/GSE113230_results/'
rows=list(csv.DictReader(open(H+'de/de_results.tsv'),delimiter='\t'))
pa=[float(r['padj']) for r in rows if r['padj'] not in('NA','')]
print('GSE113230 de rows',len(rows),'min padj',min(pa),'cols',list(rows[0].keys()))
iso=list(csv.DictReader(open(H+'de/isoform_switching.tsv'),delimiter='\t'))
print('iso rows',len(iso),'global min',min(float(r['padj_global']) for r in iso if r['padj_global'] not in('NA','')))
for thr in (0.05,0.1):
    s=[r for r in iso if r['padj_within_gene'] not in('NA','') and float(r['padj_within_gene'])<thr and abs(float(r['delta_iui']))>0.1]
    print(thr,len(s),len(set(r['gene_name'] for r in s)))
print('is_switching true',sum(r['is_switching']=='TRUE' for r in iso))
sg=[r for r in iso if r['is_switching']=='TRUE']; print('is_sw genes',len(set(r['gene_name'] for r in sg)))
# global BH significant count
print('global<0.05',sum(1 for r in iso if r['padj_global'] not in('NA','') and float(r['padj_global'])<0.05))
