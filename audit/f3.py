import csv,glob,statistics
H='/home3/choukaihsuan/'
DS=['GSE113230','GSE58135','GSE133998','SRP156355','GSE77509','GSE130078','GSE248612','GSE221107','PRJNA553289','GSE229705','GSE148036','GSE121842','GSE136569','GSE143797','GSE108735','GSE171011','GSE97239','GSE192410','GSE192849']
rows=[]
for d in DS:
    def mat(fn):
        r=list(csv.reader(open(H+d+'_results/circRNA/'+fn),delimiter='\t')); return r[0][1:],[ [float(x) for x in l[1:]] for l in r[1:]]
    cs,B=mat('count_matrix.tsv'); _,F=mat('fsj_count_matrix.tsv')
    tested=sum(1 for _ in open(H+d+'_results/de/de_results.tsv'))-1
    bs=[sum(r[i] for r in B) for i in range(len(cs))]; fs=[sum(r[i] for r in F) for i in range(len(cs))]
    rows.append((d,len(B),tested,100*tested/len(B),statistics.median(bs),statistics.median(fs),statistics.median(fs)/statistics.median(bs)))
    print(d,len(B),tested,'%.2f%%'%(100*tested/len(B)),'BSJsum_med=%d'%statistics.median(bs),'FSJsum_med=%d'%statistics.median(fs),'FSJ/BSJ=%.1f'%(statistics.median(fs)/statistics.median(bs)))
def rank(v): 
    s=sorted(range(len(v)),key=lambda i:v[i]); r=[0]*len(v)
    for k,i in enumerate(s): r[i]=k
    return r
def sp(a,b):
    ra,rb=rank(a),rank(b); n=len(a); ma=sum(ra)/n; mb=sum(rb)/n
    num=sum((x-ma)*(y-mb) for x,y in zip(ra,rb)); den=(sum((x-ma)**2 for x in ra)*sum((y-mb)**2 for y in rb))**.5
    return num/den
pr=[r[3] for r in rows]
print('Spearman(pass%, median BSJ library sum)=%.2f'%sp(pr,[r[4] for r in rows]))
print('Spearman(pass%, median FSJ library sum)=%.2f'%sp(pr,[r[5] for r in rows]))
print('Spearman(pass%, FSJ/BSJ ratio)=%.2f'%sp(pr,[r[6] for r in rows]))
