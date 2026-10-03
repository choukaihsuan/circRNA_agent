import csv,collections,sys
DS=['GSE113230','GSE58135','GSE323364','GSE133998','SRP156355','GSE77509','GSE130078','GSE248612','GSE221107','PRJNA553289','GSE229705','GSE148036','GSE121842','GSE136569','GSE143797','GSE108735','GSE171011','GSE97239','GSE192410','GSE192849']
def sig(ds):
    S=set();n=0;up=dn=0;ty=collections.Counter()
    for r in csv.DictReader(open('/home3/choukaihsuan/%s_results/de/de_results.tsv'%ds),delimiter='\t'):
        n+=1
        try: p=float(r['pvalue']); l=float(r['log2FC'])
        except Exception as e: continue
        if p<0.05 and abs(l)>1:
            S.add(r['circ_id']); up+=l>0; dn+=l<0; ty[r.get('Type','')]+=1
    return S,n,up,dn,ty
res={d:sig(d) for d in DS}
tot=0
for d in DS:
    S,n,u,dn,ty=res[d]; tot+=len(S)
    print(d,'tested',n,'sig',len(S),'up',u,'dn',dn,dict(ty))
print('sum20',tot)
def rec(ds):
    c=collections.Counter()
    for d in ds:
        for x in res[d][0]: c[x]+=1
    return {k:sum(1 for v in c.values() if v>=k) for k in (2,3,4,5,6)}, c
r19=[d for d in DS if d!='GSE323364']
print('19',rec(r19)[0],'sum',sum(len(res[d][0]) for d in r19))
print('20',rec(DS)[0])
print('19 excl 192410',rec([d for d in r19 if d!='GSE192410'])[0])
