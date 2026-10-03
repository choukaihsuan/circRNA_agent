import csv,glob,json,os,collections,statistics
H='/home3/choukaihsuan/'
DS=['GSE113230','GSE58135','GSE133998','SRP156355','GSE77509','GSE130078','GSE248612','GSE221107','PRJNA553289','GSE229705','GSE148036','GSE121842','GSE136569','GSE143797','GSE108735','GSE171011','GSE97239','GSE192410','GSE192849']
SINGLE=['GSE58135','GSE229705','GSE192410','GSE192849']
def rd(p): return list(csv.DictReader(open(p),delimiter='\t'))
sig={};ty={};info={}
for d in DS:
    S={};T=collections.Counter()
    for r in rd(H+d+'_results/de/de_results.tsv'):
        p=float(r['pvalue']);l=float(r['log2FC'])
        if p<0.05 and abs(l)>1: S[r['circ_id']]=r['Type']; T[r['Type']]+=1
    sig[d]=S; ty[d]=T
    # consensus rows n_tools
    rows1=rows2=0; dec_support=set()
    sg=rd(H+d+'_results/../'+d+'_results/circRNA/count_matrix.tsv') if False else None
    two_ids=set(); allrows=0
    for f in glob.glob(H+d+'_results/circRNA/SRR*/consensus_summary.tsv'):
        for r in csv.DictReader(open(f),delimiter='\t'):
            allrows+=1
            if r['n_tools']=='1': rows1+=1
            else:
                rows2+=1
                if r['circ_id'] in S: two_ids.add(r['circ_id'])
    # read length / depth
    cfg=open(H+'circRNA_agent/config/projects/%s.yaml'%d).read()
    import re
    grp=H+'circRNA_agent/'+re.search(r'^groups:\s*(\S+)',cfg,re.M).group(1)
    g=rd(grp.replace('.tsv','.tsv')) if grp.endswith('.tsv') else list(csv.DictReader(open(grp)))
    cond=collections.Counter(x['condition'] for x in g); pid='patient_id' in g[0]
    ids=set(x['srr_id'] for x in g); rl=[]; rr=[]
    for j in glob.glob(H+d+'_results/qc/fastp/*.json'):
        n=os.path.basename(j)[:-5]
        if n in ids:
            s=json.load(open(j))['summary']['before_filtering']; rl.append(s['read1_mean_length']); rr.append(s['total_reads'])
    info[d]=dict(n_t=cond.get('tumor',0),n_n=cond.get('normal',0),pid=pid,readlen='%d-%d'%(min(rl),max(rl)) if rl else 'NA',reads_M='%.0f-%.0f'%(min(rr)/1e6,max(rr)/1e6) if rr else 'NA',rows=allrows,single_pct=100*rows1/allrows if allrows else float('nan'),dec_two=len(two_ids),
        cfg_tools=re.search(r'tools:\s*\n((?:\s*-\s*\w+\n)+)',cfg).group(1).replace('\n','').replace(' ','').replace('-','+').strip('+'),cfg_min=re.search(r'min_tools:\s*(\d)',cfg).group(1))
def tiers(ds):
    c=collections.Counter()
    for d in ds:
        for k in sig[d]: c[k]+=1
    return c
def summ(c): return {k:sum(1 for v in c.values() if v>=k) for k in (2,3,4,5)}
c19=tiers(DS); c15=tiers([d for d in DS if d not in SINGLE]); c18=tiers([d for d in DS if d!='GSE192410'])
print('19',summ(c19)); print('excl192410',summ(c18)); print('15',summ(c15)); print('15 max',max(c15.values()))
# 7 candidates
c7=[k for k,v in c19.items() if v>=5]
print('7 candidates:')
for k in c7:
    ds=[d for d in DS if k in sig[d]]
    print(k,'19set',len(ds),'15set',len([d for d in ds if d not in SINGLE]),'single-tool datasets among its 5:',[d for d in ds if d in SINGLE])
# per-dataset table
out=open('/home/choukaihsuan/audit_tmp/per_dataset.tsv','w')
w=csv.writer(out,delimiter='\t')
w.writerow(['accession','n_tumor','n_normal','paired_model(patient_id)','read_len_mean_range','reads_per_sample_M(range,both mates)','config_tools','config_min_tools','consensus_rows','single_tool_rows_pct','effective_min_tools','DEC_total','DEC_with_any_2tool_support','TypeI_pct','DEC_in_n>=2','DEC_in_n>=3','DEC_in_n>=4','DEC_in_n>=5','DEC_in_n>=2_excl4single','DEC_in_n>=3_excl4single'])
for d in DS:
    i=info[d]; n=len(sig[d]); 
    eff='2' if i['single_pct']<1 else ('1' if i['single_pct']>99 else 'mixed(%.1f%% single)'%i['single_pct'])
    if 0<i['single_pct']<1: eff='2 (%.2f%% single-tool rows)'%i['single_pct']
    row=[d,i['n_t'],i['n_n'],'yes' if i['pid'] else 'no',i['readlen'],i['reads_M'],i['cfg_tools'],i['cfg_min'],i['rows'],'%.1f'%i['single_pct'],eff,n,i['dec_two'],'%.1f'%(100*ty[d]['Type_I']/n)]
    for k in (2,3,4,5): row.append(sum(1 for x in sig[d] if c19[x]>=k))
    row.append(sum(1 for x in sig[d] if c15[x]>=2) if d not in SINGLE else 'NA'); row.append(sum(1 for x in sig[d] if c15[x]>=3) if d not in SINGLE else 'NA')
    w.writerow(row)
out.close()
print(open('/home/choukaihsuan/audit_tmp/per_dataset.tsv').read())
# isoform tables
def iso(d):
    rows=rd(H+d+'_results/de/isoform_switching.tsv')
    res=[]
    for kind,col in (('global','padj_global'),('within-gene','padj_within_gene')):
        for thr in (0.05,0.1):
            s=[r for r in rows if r[col] not in('NA','') and float(r[col])<thr and abs(float(r['delta_iui']))>0.1]
            res.append((d,kind,thr,len(set(r['gene_name'] for r in s)),len(s)))
    return res,len(rows)
for d in ('GSE113230','GSE133998'):
    r,n=iso(d); print('iso rows',d,n); [print(x) for x in r]
