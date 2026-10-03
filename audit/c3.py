import csv,collections
H='/home3/choukaihsuan/'
for d in ('GSE113230','GSE133998'):
    rows=list(csv.DictReader(open(H+d+'_results/de/isoform_switching.tsv'),delimiter='\t'))
    ok=lambda r,c: r[c] not in('NA','')
    big=[r for r in rows if ok(r,'p_value') and abs(float(r['delta_iui']))>0.1]
    raw=[r for r in big if float(r['p_value'])<0.05]
    for thr in (0.05,0.1):
        w=[r for r in big if ok(r,'padj_within_gene') and float(r['padj_within_gene'])<thr]
        two=[r for r in w if int(float(r['n_isoforms']))==2]
        print(d,'within<',thr,len(w),'of which from 2-isoform genes',len(two),'| genes with n_iso>=3:',len(set(r['gene_name'] for r in w if int(float(r['n_isoforms']))>=3)))
    print(d,'raw p<0.05 & |dIUI|>0.1 (no correction):',len(raw),'genes',len(set(r['gene_name'] for r in raw)))
    nis=collections.Counter(int(float(r['n_isoforms'])) for r in rows if ok(r,'p_value'))
    print(d,'rows with p by n_isoforms (rows):',dict(sorted(nis.items())[:6]))
    mp=min(float(r['p_value']) for r in rows if ok(r,'p_value')); print(d,'min raw p',mp)
