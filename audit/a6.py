import csv,glob,json,os,statistics,collections
H='/home3/choukaihsuan/'
DS=['GSE113230','GSE58135','GSE323364','GSE133998','SRP156355','GSE77509','GSE130078','GSE248612','GSE221107','PRJNA553289','GSE229705','GSE148036','GSE121842','GSE136569','GSE143797','GSE108735','GSE171011','GSE97239','GSE192410','GSE192849']
for d in DS:
    R=H+d+'_results/'
    cfg=open(H+'circRNA_agent/config/projects/%s.yaml'%d).read()
    import re
    grp=H+'circRNA_agent/'+re.search(r'^groups:\s*(\S+)',cfg,re.M).group(1)
    g=list(csv.DictReader(open(grp)))
    cond=collections.Counter(r['condition'] for r in g); pid='patient_id' in g[0]
    reads=[];rl=[];q30=[]
    for j in glob.glob(R+'qc/fastp/*.json'):
        x=json.load(open(j))
        s=x['summary']
        reads.append(s['before_filtering']['total_reads']/2); rl.append(s['before_filtering']['read1_mean_length']); q30.append(s['before_filtering']['q30_rate'])
    sr=set(r['srr_id'] for r in g)
    # only fastp for samples in groups
    reads2=[];rl2=[]
    for j in glob.glob(R+'qc/fastp/*.json'):
        n=os.path.basename(j).replace('.json','')
        if n in sr:
            x=json.load(open(j))['summary']['before_filtering']; reads2.append(x['total_reads']/2); rl2.append(x['read1_mean_length'])
    cm=sum(1 for _ in open(R+'circRNA/count_matrix.tsv'))-1
    print(d,dict(cond),'pid',pid,'ncm',cm,'pairs_min_reads(M)=%.1f'%(min(reads2)/1e6 if reads2 else -1),'max=%.1f'%(max(reads2)/1e6 if reads2 else -1),'readlen=%s'%(sorted(set(round(x) for x in rl2))[:3]),'nfastp',len(reads2))
