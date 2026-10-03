import csv,glob,collections
H='/home3/choukaihsuan/'
DS=['GSE113230','GSE58135','GSE133998','SRP156355','GSE77509','GSE130078','GSE248612','GSE221107','PRJNA553289','GSE229705','GSE148036','GSE121842','GSE136569','GSE143797','GSE108735','GSE171011','GSE97239','GSE192410','GSE192849']
for d in DS:
    c=collections.Counter(); per=[]
    for f in sorted(glob.glob(H+d+'_results/circRNA/SRR*/consensus_summary.tsv')):
        n=collections.Counter(r['n_tools'] for r in csv.DictReader(open(f),delimiter='\t'))
        c.update(n); per.append((n.get('2',0)))
    tot=sum(c.values())
    print(d,'samples',len(per),'rows',tot,'n_tools dist',dict(c),'samples_with_2tool',sum(1 for x in per if x>0))
