import glob,re,collections,statistics
H='/home3/choukaihsuan/'
for d in ('GSE133998','GSE113230','GSE108735'):
    by=collections.defaultdict(list)
    for g in glob.glob(H+d+'_results/circRNA/SRR*/SRR*.gtf'):
        for l in open(g):
            if l.startswith('#'): continue
            f=l.split('\t')
            if len(f)<9 or f[2]!='circRNA': continue
            t=re.search(r'circ_type "([^"]+)"',f[8]); t=t.group(1) if t else 'NA'
            b=float(re.search(r'bsj ([\d.]+)',f[8]).group(1)); s=float(re.search(r'fsj ([\d.]+)',f[8]).group(1))
            by[t].append((b,s))
    print(d)
    for t,v in sorted(by.items()):
        n=len(v); z=sum(1 for b,s in v if s==0); fs=[s for b,s in v]
        print('  circ_type=%-11s n=%7d  FSJ==0: %5.1f%%  median FSJ=%g  mean FSJ=%.1f  | BSJ>=5 subset: n=%d FSJ==0 %5.1f%%'%(t,n,100*z/n,statistics.median(fs),sum(fs)/n,sum(1 for b,s in v if b>=5),100*sum(1 for b,s in v if b>=5 and s==0)/max(1,sum(1 for b,s in v if b>=5))))
