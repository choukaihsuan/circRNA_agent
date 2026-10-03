import csv,collections
H='/home3/choukaihsuan/'
B=H+'GSE113230_results/benchmark/'
rows=list(csv.DictReader(open(B+'pr_curve.tsv'),delimiter='\t'))
for m in ['Our_adaptive','CirComPara2_4tools','nfcore_3tools']:
    pts=[(float(r['Recall']),float(r['Precision'])) for r in rows if r['method']==m]
    pts.sort()
    a=sum((pts[i+1][0]-pts[i][0])*(pts[i+1][1]+pts[i][1])/2 for i in range(len(pts)-1))
    # variant incl (0,first precision)
    a2=a+pts[0][0]*pts[0][1]
    print(m,'trapz %.4f'%a,'with-origin-rect %.4f'%a2,'n',len(pts))
fc=H+'GSE55872_results/circRNA/SRR444655/find_circ/splice_sites.bed'
c=collections.Counter()
for l in open(fc):
    p=l.rstrip('\n').split('\t')
    c[p[17] if len(p)>17 else 'short']+=1
print(c.most_common(8))
gt=list(csv.DictReader(open(B+'rnaser/ground_truth.tsv'),delimiter='\t'))
b=[float(r['bsj_total']) for r in gt]
print('GT rows',len(b),'frac bsj_total<=2: %.3f'%(sum(1 for x in b if x<=2)/len(b)), 'frac<=2 among pos %.3f'%(sum(1 for r in gt if r['is_true']=='1' and float(r['bsj_total'])<=2)/sum(1 for r in gt if r['is_true']=='1')))
print('frac total bsj in {1,2}: %.3f'%(sum(1 for x in b if x in (1.0,2.0))/len(b)))
er=[float(r['enrichment_ratio']) for r in gt]
print('ER==1.5 count',sum(1 for x in er if x==1.5),'ER==0.5',sum(1 for x in er if x==0.5), 'is_true=1 min ER',min(float(r['enrichment_ratio']) for r in gt if r['is_true']=='1'),'is_true=0 max',max(float(r['enrichment_ratio']) for r in gt if r['is_true']=='0'))
