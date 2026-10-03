import yaml,csv,os
H='/home3/choukaihsuan/circRNA_agent/config/projects/'
DS=['GSE113230','GSE58135','GSE133998','SRP156355','GSE77509','GSE130078','GSE248612','GSE221107','PRJNA553289','GSE229705','GSE148036','GSE121842','GSE136569','GSE143797','GSE108735','GSE171011','GSE97239','GSE192410','GSE192849']
def g(c,*k):
    for x in k:
        if not isinstance(c,dict) or x not in c: return '<absent>'
        c=c[x]
    return c
w=csv.writer(open('/home/choukaihsuan/audit_tmp/config_diff.tsv','w'),delimiter='\t')
w.writerow(['accession','consensus.tools','consensus.min_tools','de.method','circbase_file','genome.hisat2_index(basename)','threads','download.max_concurrent','log_dir','notify.email_to','config_mtime'])
import time
for d in DS:
    c=yaml.safe_load(open(H+d+'.yaml'))
    cb=g(c,'circbase_file'); cb='local_cache' if cb else 'auto_download'
    hi=os.path.basename(str(g(c,'genome','hisat2_index')))
    w.writerow([d,'+'.join(g(c,'consensus','tools')),g(c,'consensus','min_tools'),g(c,'de','method'),cb,hi,g(c,'threads'),g(c,'download','max_concurrent'),g(c,'log_dir'),g(c,'notify','email_to'),time.strftime('%F',time.localtime(os.path.getmtime(H+d+'.yaml')))])
print(open('/home/choukaihsuan/audit_tmp/config_diff.tsv').read())
