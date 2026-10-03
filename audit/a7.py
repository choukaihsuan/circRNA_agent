import yaml,collections
H='/home3/choukaihsuan/circRNA_agent/config/projects/'
DS=['GSE113230','GSE58135','GSE323364','GSE133998','SRP156355','GSE77509','GSE130078','GSE248612','GSE221107','PRJNA553289','GSE229705','GSE148036','GSE121842','GSE136569','GSE143797','GSE108735','GSE171011','GSE97239','GSE192410','GSE192849']
def flat(d,p=''):
    o={}
    for k,v in d.items():
        if isinstance(v,dict): o.update(flat(v,p+k+'.'))
        else: o[p+k]=repr(v)
    return o
F={d:flat(yaml.safe_load(open(H+d+'.yaml'))) for d in DS}
keys=sorted(set(k for f in F.values() for k in f))
for k in keys:
    vals={d:F[d].get(k,'<absent>') for d in DS if d!='GSE323364'}
    c=collections.Counter(vals.values())
    if len(c)>1:
        if len(c)>6: print('##',k,'(%d distinct values; path-like)'%len(c), list(c)[:2]); continue
        print('##',k)
        for v,n in c.most_common():
            print('   ',n,v if len(v)<120 else v[:120],[d for d in vals if vals[d]==v] if n<=8 else '')
print('keys total',len(keys))
print('GSE323364 vs others differing keys:')
for k in keys:
    if F['GSE323364'].get(k)!=F['GSE113230'].get(k) and not any(x in k for x in('dir','metadata','groups','project_id','study_title','_file','fasta')): print(k,F['GSE323364'].get(k),'|',F['GSE113230'].get(k))
