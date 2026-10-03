import csv,math,re
H='/home3/choukaihsuan/'
DS=['GSE113230','GSE58135','GSE133998','SRP156355','GSE77509','GSE130078','GSE248612','GSE221107','PRJNA553289','GSE229705','GSE148036','GSE121842','GSE136569','GSE143797','GSE108735','GSE171011','GSE97239','GSE192410','GSE192849']
tot=[0,0,0,0]
print('dataset DEC | DEC with >=1 "FSJ-collapse" sample (BSJ>=5 & FSJ<=2) | sign(log2FC) != sign(BSJ-only CPM log-ratio) | both')
for d in DS:
    def mat(fn):
        r=list(csv.reader(open(H+d+'_results/circRNA/'+fn),delimiter='\t')); return r[0][1:],{l[0]:[float(x) for x in l[1:]] for l in r[1:]}
    cs,B=mat('count_matrix.tsv'); _,F=mat('fsj_count_matrix.tsv')
    cfg=open(H+'circRNA_agent/config/projects/%s.yaml'%d).read(); grp=H+'circRNA_agent/'+re.search(r'^groups:\s*(\S+)',cfg,re.M).group(1)
    cond={r['srr_id']:r['condition'] for r in csv.DictReader(open(grp))}
    ti=[i for i,s in enumerate(cs) if cond.get(s)=='tumor']; ni=[i for i,s in enumerate(cs) if cond.get(s)=='normal']
    lib=[sum(v[i] for v in B.values()) for i in range(len(cs))]
    n=c=s_=cs_=0
    for r in csv.DictReader(open(H+d+'_results/de/de_results.tsv'),delimiter='\t'):
        p=float(r['pvalue']);l=float(r['log2FC'])
        if not(p<0.05 and abs(l)>1): continue
        n+=1; b=B[r['circ_id']]; f=F[r['circ_id']]
        col=any(b[i]>=5 and f[i]<=2 for i in range(len(cs)))
        mt=sum(b[i]/lib[i] for i in ti)/len(ti); mn=sum(b[i]/lib[i] for i in ni)/len(ni)
        bs=math.log2((mt+1e-9)/(mn+1e-9)); dis=(bs>0)!=(l>0)
        c+=col; s_+=dis; cs_+=(col and dis)
    print(d,n,c,s_,cs_); tot[0]+=n;tot[1]+=c;tot[2]+=s_;tot[3]+=cs_
print('TOTAL',tot)
