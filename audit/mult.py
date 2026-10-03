import csv,collections
B='/home/choukaihsuan/GSE113230_results/benchmark/'
gt=[]
for r in csv.DictReader(open(B+'rnaser/ground_truth.tsv'),delimiter='\t'):
    gt.append((r['chr'],int(r['start']),int(r['end']),int(r['is_true']),float(r['enrichment_ratio'])))
print('gt rows',len(gt),collections.Counter(g[3] for g in gt))
det=[l.split('\t') for l in open(B+'detection/our_method.bed') if l.strip()]
det=[(d[0],int(d[1]),int(d[2])) for d in det]
print('detections',len(det))
idx=collections.defaultdict(list)
for i,g in enumerate(gt): idx[g[0]].append(i)
def m(d,g,s=10): return max(abs(d[1]-g[1]),abs(d[2]-g[2]))<=s
none=0; cnt=collections.Counter(); inter_only=0; multi=0; mult_gt_total=0
matched_gt=set()
for d in det:
    ms=[i for i in idx[d[0]] if m(d,gt[i])]
    matched_gt.update(ms)
    labs=[gt[i][3] for i in ms]
    if not ms: none+=1
    cnt[(len(ms))]+=1
    if ms and all(l==-1 for l in labs): inter_only+=1
    if len(ms)>1: multi+=1
print('detections matching 0 GT:',none,' by #GT matched:',dict(cnt))
print('detections matching only intermediate GT:',inter_only,' detections matching >1 GT:',multi)
ml=collections.Counter(gt[i][3] for i in matched_gt)
print('distinct GT entries matched by label (1=pos,0=neg,-1=inter):',dict(ml),'total',len(matched_gt))
# detections whose best match includes pos/neg
dp=sum(1 for d in det if any(gt[i][3]==1 for i in idx[d[0]] if m(d,gt[i])))
dn=sum(1 for d in det if any(gt[i][3]==0 for i in idx[d[0]] if m(d,gt[i])))
print('detections touching >=1 pos GT',dp,' >=1 neg GT',dn)
# exact-match detections
ex=sum(1 for d in det if any(gt[i][1]==d[1] and gt[i][2]==d[2] for i in idx[d[0]]))
print('exact',ex)
