import csv
H='/home3/choukaihsuan/'
DS=['GSE113230','GSE58135','GSE133998','SRP156355','GSE77509','GSE130078','GSE248612','GSE221107','PRJNA553289','GSE229705','GSE148036','GSE121842','GSE136569','GSE143797','GSE108735','GSE171011','GSE97239','GSE192410','GSE192849']
tot=[0,0,0,0]
for d in DS:
    a=b=c=e=0
    for r in csv.DictReader(open(H+d+'_results/de/de_results.tsv'),delimiter='\t'):
        p=float(r['pvalue']);l=float(r['log2FC'])
        if p<0.05 and abs(l)>1:
            try: pf=float(r['pvalue_fsj']);lf=float(r['logFC_fsj'])
            except: pf=lf=None
            code_II = pf is not None and pf<0.05
            paper_II = code_II and abs(lf)>1
            a+=code_II; b+=paper_II; c+=1
            e+= (r['Type']=='Type_II')
    print(d,'sig',c,'TypeII(file)',e,'TypeII(code rule recomputed)',a,'TypeII(paper rule P<.05&|lfc|>1)',b,'TypeI share file %.1f%% paper-rule %.1f%%'%(100*(c-e)/c,100*(c-b)/c))
    tot[0]+=c;tot[1]+=e;tot[2]+=b
print(tot)
