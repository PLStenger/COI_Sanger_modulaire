from common import *
import math, itertools
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def distance(a,b):
    pairs=[(x,y) for x,y in zip(a,b) if x in 'ACGT' and y in 'ACGT']
    n=len(pairs)
    if not n:return dict(comparable_bp=0,differences=None,p_distance=None,K2P=None)
    transitions=sum(x!=y and {x,y} in [{'A','G'},{'C','T'}] for x,y in pairs)
    diff=sum(x!=y for x,y in pairs);p=transitions/n;q=(diff-transitions)/n
    u,v=1-2*p-q,1-2*q
    k=-.5*math.log(u)-.25*math.log(v) if u>0 and v>0 else None
    return dict(comparable_bp=n,differences=diff,p_distance=diff/n,K2P=k)


def main():
    c,out=load();qfile=validated_query(c,out)
    provenance=json.loads((out/'04_references'/'provenance.json').read_text())
    source=out/'04_references'/'oriented_COI.fasta';metafile=out/'04_references'/'metadata.json'
    for key,current in [('query_sha256',sha(qfile)),('panel_sha256',sha(ROOT/c['references_file'])),
                        ('fasta_sha256',sha(source)),('metadata_sha256',sha(metafile)),('config_sha256',config_sha(c))]:
        require(provenance[key]==current,'Références/configuration changées : relancer 04 puis 05.')
    metadata=json.loads(metafile.read_text());qid='QUERY_'+c['sample']
    raw=out/'05_alignment'/'mafft_full.fasta'
    command([exe('mafft'),'--localpair','--maxiterate','1000','--thread',str(c['threads']),source],out/'logs'/'05_mafft.log',stdout=raw)
    aln=AlignIO.read(raw,'fasta');q=str(next(r for r in aln if r.id==qid).seq).upper()
    # Conserver toute la région query, y compris N ; les insertions par rapport à la query sont auditées à part.
    columns=[i for i,b in enumerate(q) if b!='-'];qpos={col:i+1 for i,col in enumerate(columns)}
    qregion=[];coverage=[]
    for r in aln:
        seq=''.join(str(r.seq)[i] for i in columns).upper()
        n=sum(x in 'ACGT' and y in 'ACGT' for x,y in zip(seq,q.replace('-','')))
        keep=(r.id==qid or (n>=c['min_reference_comparable_bp'] and n/sum(b in 'ACGT' for b in q)>=c['min_reference_coverage']))
        coverage.append(dict(id=r.id,comparable_to_query=n,retained=keep))
        if keep:qregion.append(SeqRecord(Seq(seq),id=r.id,description=''))
    tsv(out/'05_alignment'/'coverage.tsv',coverage,list(coverage[0]))
    require(len(qregion)>=5,'Trop peu de séquences après MAFFT.')
    require(any(metadata[r.id]['role']=='outgroup' for r in qregion if r.id!=qid),'Outgroup exclu après alignement : revoir le panel.')
    SeqIO.write(qregion,out/'05_alignment'/'query_region.fasta','fasta')
    common=[i for i in range(len(columns)) if all(str(r.seq)[i] in 'ACGT' for r in qregion)]
    require(len(common)>=c['min_alignment_columns'],f'Seulement {len(common)} sites ACGT communs : augmenter recouvrement ou revoir les références, pas abaisser automatiquement le seuil.')
    analysis=[SeqRecord(Seq(''.join(str(r.seq)[i] for i in common)),id=r.id,description='') for r in qregion]
    final=out/'05_alignment'/'analysis_common_sites.fasta';SeqIO.write(analysis,final,'fasta')
    tsv(out/'05_alignment'/'column_map.tsv',[
        dict(query_position=i+1,mafft_column=col+1,retained_common=i in set(common)) for i,col in enumerate(columns)],
        ['query_position','mafft_column','retained_common'])
    checks=[]
    for r in aln:
        if r.id not in {x.id for x in qregion}:continue
        seq=str(r.seq).upper();non=[i for i in range(len(seq)) if seq[i]!='-']
        qnon=[i for i,b in enumerate(q) if b!='-']
        lo=max(min(non),min(qnon));hi=min(max(non),max(qnon))+1
        internal=sum(seq[i]=='-' for i in range(lo,hi)) if non else 0
        checks.append(dict(id=r.id,internal_gap_columns=internal,insertions_vs_query=sum(seq[i]!='-' and q[i]=='-' for i in range(lo,hi))))
    tsv(out/'05_alignment'/'indel_audit.tsv',checks,list(checks[0]))
    # Contrôle de traduction sur la région query AVANT suppression des sites communs.
    # La délétion de colonnes pour les distances/ML n'est PAS un alignement codon.
    approved=json.loads((out/'consensus_approved.json').read_text())
    off=approved['coding_offset_0based']; translation_audit=[]
    for r in qregion:
        seq=str(r.seq).replace('-','N')
        frag=seq[off:];frag=frag[:len(frag)//3*3]
        protein=str(Seq(frag).translate(table=c['genetic_code']))
        translation_audit.append(dict(id=r.id,offset_0based=off,stop_count=protein.count('*'),unknown_aa=protein.count('X'),protein=protein))
    tsv(out/'05_alignment'/'translation_audit.tsv',translation_audit,list(translation_audit[0]))
    queries=[]
    for dataset,records in [('pairwise_query_region',qregion),('common_sites',analysis)]:
        seqs={r.id:str(r.seq) for r in records}
        for name,s in seqs.items():
            if name==qid:continue
            d=distance(seqs[qid],s);m=metadata[name]
            queries.append(dict(dataset=dataset,id=name,accession=m['accession'],species=m['organism'],role=m['role'],**d))
        for metric in ['p_distance','K2P']:
            names=list(seqs)
            with (out/'05_alignment'/f'{dataset}_{metric}_matrix.csv').open('w',newline='') as h:
                w=csv.writer(h);w.writerow(['id']+names)
                for n in names:w.writerow([n]+[distance(seqs[n],seqs[m])[metric] for m in names])
    queries.sort(key=lambda x:(x['dataset'],x['p_distance'] if x['p_distance'] is not None else float('inf')))
    tsv(out/'05_alignment'/'query_distances.tsv',queries,list(queries[0]))
    values=sorted([x for x in queries if x['dataset']=='common_sites'],key=lambda x:x['p_distance'])
    fig,ax=plt.subplots(figsize=(10,max(4,len(values)*.4+1)))
    ax.barh(range(len(values)),[100*x['p_distance'] for x in values],color=['#888888' if x['role']=='outgroup' else '#2678B2' for x in values])
    ax.set_yticks(range(len(values)));ax.set_yticklabels([f"{x['species']} | {x['accession']}" for x in values]);ax.invert_yaxis()
    ax.set_xlabel(f"Distance p à la query (%) — mêmes {len(common)} sites pour toutes les paires")
    ax.set_title(c['sample']+' — distances COI, sans seuil automatique d’espèce');fig.tight_layout()
    fig.savefig(out/'07_plots'/'distances_common.pdf',bbox_inches='tight');fig.savefig(out/'07_plots'/'distances_common.png',dpi=200,bbox_inches='tight');plt.close(fig)
    for x in values[:10]:print(f"{x['accession']} | {x['species']} | {x['differences']}/{x['comparable_bp']} différences")
    fingerprint=hashlib.sha256((sha(final)+sha(metafile)+sha(qfile)+sha(ROOT/c['references_file'])).encode()).hexdigest()
    dump(out/'05_alignment'/'manifest.json',dict(fingerprint=fingerprint,alignment_sha256=sha(final),
        metadata_sha256=sha(metafile),query_sha256=sha(qfile),panel_sha256=sha(ROOT/c['references_file']),
        config_sha256=config_sha(c),n_sequences=len(analysis),query_columns=len(columns),common_columns=len(common)))
    print('Distances prêtes. Revoir audit GenBank, gaps, alignement complet et sites communs ; aucun arbre lancé.')

if __name__=='__main__':main()
