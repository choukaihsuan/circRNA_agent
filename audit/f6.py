import glob,re,collections
H='/home3/choukaihsuan/'
for d in ('GSE133998','GSE113230','GSE108735'):
    z_sib=z=p_sib=p=0; z_sib_pos=0
    for g in glob.glob(H+d+'_results/circRNA/SRR*/SRR*.gtf'):
        rows=[]
        for l in open(g):
            if l.startswith('#'): continue
            f=l.split('\t')
            if len(f)<9 or f[2]!='circRNA': continue
            b=float(re.search(r'bsj ([\d.]+)',f[8]).group(1)); s=float(re.search(r'fsj ([\d.]+)',f[8]).group(1))
            rows.append((f[0],int(f[3]),int(f[4]),b,s))
        st=collections.defaultdict(list); en=collections.defaultdict(list)
        for i,r in enumerate(rows): st[(r[0],r[1])].append(i); en[(r[0],r[2])].append(i)
        for i,(c,a,e,b,s) in enumerate(rows):
            if b<5: continue
            sib=[j for j in st[(c,a)]+en[(c,e)] if j!=i]
            has=len(sib)>0
            if s==0:
                z+=1; z_sib+=has; z_sib_pos+= any(rows[j][4]>0 for j in sib)
            else:
                p+=1; p_sib+=has
    print(d,'cells BSJ>=5: FSJ==0 n=%d, with boundary-sharing sibling %.1f%% (sibling with FSJ>0: %.1f%%) | FSJ>0 n=%d, with sibling %.1f%%'%(z,100*z_sib/max(z,1),100*z_sib_pos/max(z,1),p,100*p_sib/max(p,1)))
