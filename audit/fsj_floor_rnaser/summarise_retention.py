"""How much of the baseline DEC set survives an FSJ floor? (inputs: fsj_floor_sweep*_<GSE>.tsv from ~/audit_tmp/g3.R / g3p.R on the server)
GSE113230 and GSE133998 sweeps (unpaired = their production design) come from g3.R; GSE97239 and GSE136569 (RNase R in protocol; paired = production design) from g3p.R.
Baseline row (floor 0) must equal the published n_sig: 240, 79, 140, 199 -- checked below."""
import csv, os
H = os.path.dirname(os.path.abspath(__file__)); PUB = {'GSE113230': 240, 'GSE133998': 79, 'GSE97239': 140, 'GSE136569': 199}
FILES = {'GSE113230': 'fsj_floor_sweep_GSE113230.tsv', 'GSE133998': 'fsj_floor_sweep_GSE133998.tsv', 'GSE97239': 'fsj_floor_sweep_paired_GSE97239.tsv', 'GSE136569': 'fsj_floor_sweep_paired_GSE136569.tsv'}
CFG = [('1', 'all'), ('2', 'all'), ('5', 'all'), ('1', 'half_per_group'), ('5', 'half_per_group'), ('20', 'half_per_group')]
out = []
for d, fn in FILES.items():
    rows = list(csv.DictReader(open(os.path.join(H, fn)), delimiter='\t')); b = [r for r in rows if r['floor'] == '0'][0]
    assert int(b['n_sig']) == PUB[d], (d, b['n_sig'], PUB[d])
    base = int(b['n_sig']); base_up = int(b['up'])
    for fl, ru in CFG:
        r = [x for x in rows if x['floor'] == fl and x['rule'] == ru][0]
        out.append([d, 'yes' if d in ('GSE97239', 'GSE136569') else 'no', fl, ru, r['loci_passing'], b['loci_passing'], base, r['n_sig'], r['baseline_sig_retained'], '%.1f' % (100*int(r['baseline_sig_retained'])/base), '%s/%s' % (b['up'], b['down']), '%s/%s' % (r['up'], r['down']), r['pct_typeI'], b['pct_typeI']])
with open(os.path.join(H, 'retention_summary.tsv'), 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n')
    w.writerow(['dataset', 'rnaseR_in_protocol', 'fsj_floor', 'rule', 'loci_passing', 'loci_tested_baseline', 'baseline_n_sig', 'n_sig_after', 'baseline_sig_retained', 'pct_baseline_retained', 'up/down_baseline', 'up/down_after', 'pct_typeI_after', 'pct_typeI_baseline']); w.writerows(out)
for r in out: print('\t'.join(map(str, r)))
