import glob,collections,os
H='/home3/choukaihsuan/'
DS=['GSE113230','GSE58135','GSE323364','GSE133998','SRP156355','GSE77509','GSE130078','GSE248612','GSE221107','PRJNA553289','GSE229705','GSE148036','GSE121842','GSE136569','GSE143797','GSE108735','GSE171011','GSE97239','GSE192410','GSE192849']
for d in DS:
    mins=[];hdr=0;n=0
    for f in sorted(glob.glob(H+d+'_results/circRNA/SRR*/DCC/CircRNACount')):
        vals=[]
        for l in open(f):
            p=l.rstrip('\n').split('\t')
            if len(p)>=4 and p[0]!='Chr':
                try: vals.append(int(float(p[3])))
                except: pass
        n+=1
        if not vals: hdr+=1
        else: mins.append(min(vals))
    cc=[os.path.getsize(f) for f in glob.glob(H+d+'_results/circRNA/SRR*/DCC/CircCoordinates')]
    print(d,'samples',n,'empty',hdr,'min-count-per-sample',sorted(set(mins)),'CircCoordinates sizes<100B:',sum(1 for s in cc if s<100))
