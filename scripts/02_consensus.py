from common import *
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import PairwiseAligner
import numpy as np


def aligner():
    a=PairwiseAligner();a.mode='global';a.match_score=2;a.mismatch_score=-3
    a.open_gap_score=-10;a.extend_gap_score=-1;a.end_gap_score=0
    return a


def columns(alignment):
    coords=alignment.coordinates
    for k in range(coords.shape[1]-1):
        f0,f1=map(int,coords[0,k:k+2]);r0,r1=map(int,coords[1,k:k+2])
        if f1>f0 and r1>r0:
            require(f1-f0==r1-r0,'Segment inattendu dans l’alignement.')
            yield from zip(range(f0,f1),range(r0,r1))
        elif f1>f0:yield from ((i,None) for i in range(f0,f1))
        else:yield from ((None,j) for j in range(r0,r1))


def build(f,r,c):
    fs,rs=str(f.seq).upper(),str(r.seq).upper()
    fq,rq=f.letter_annotations['phred_quality'],r.letter_annotations['phred_quality']
    a=aligner().align(fs,rs)[0];cols=list(columns(a))
    double=[k for k,(i,j) in enumerate(cols) if i is not None and j is not None]
    require(bool(double),'Pas de recouvrement des deux lectures.')
    lo,hi=min(double),max(double)+1
    rows=[];bases=[];position=0;indels=0;pairs=[]
    for k,(i,j) in enumerate(cols):
        fb,rb=fs[i] if i is not None else '-',rs[j] if j is not None else '-'
        qf,qr=fq[i] if i is not None else 0,rq[j] if j is not None else 0
        retained=lo<=k<hi
        if not retained:call='';status='outside_double_coverage'
        elif i is None or j is None:call='N';status='internal_indel';indels+=1
        else:
            vf,vr=fb in 'ACGT' and qf>=c['min_base_q'],rb in 'ACGT' and qr>=c['min_base_q']
            if vf and vr and fb==rb:call=fb;status='two_reads_agree'
            elif vf and vr:call='N';status='high_quality_conflict'
            elif vf:call=fb;status='forward_only_usable'
            elif vr:call=rb;status='reverse_only_usable'
            else:call='N';status='no_usable_base'
            if vf and vr:pairs.append((fb,rb))
        if retained:position+=1;bases.append(call)
        rows.append(dict(alignment_column=k+1,consensus_position=position if retained else '',
                         f_trimmed_position=i+1 if i is not None else '',r_oriented_position=j+1 if j is not None else '',
                         f_base=fb,r_base=rb,f_q=qf,r_q=qr,call=call,status=status))
    seq=''.join(bases)
    pid=sum(x==y for x,y in pairs)/len(pairs) if pairs else 0
    return seq,rows,dict(overlap_columns=hi-lo,high_quality_paired_sites=len(pairs),pair_identity=pid,
                          internal_indel_columns=indels,consensus_bp=len(seq),n_count=seq.count('N'),
                          n_fraction=seq.count('N')/len(seq) if seq else 1),a


def main():
    c,out=load();p=out/'01_qc'/'provenance.json'
    require(p.exists(),'Lancer 01_qc_ab1.py avant.')
    provenance=json.loads(p.read_text())
    require(provenance['config_sha256']==config_sha(c),'Configuration modifiée : relancer 01 et 02.')
    for role,key in [('forward','forward_ab1'),('reverse','reverse_ab1')]:
        require(sha(Path(c['raw_dir'])/c[key])==provenance['raw_sha256'][role],'AB1 changé : relancer 01.')
        require(sha(out/'01_qc'/f'{role}.fastq')==provenance['outputs'][role],'FASTQ changé : relancer 01.')
    f=SeqIO.read(out/'01_qc'/'forward.fastq','fastq');rawr=SeqIO.read(out/'01_qc'/'reverse.fastq','fastq')
    a=aligner();s0=a.score(str(f.seq).upper(),str(rawr.seq).upper());s1=a.score(str(f.seq).upper(),revcomp(str(rawr.seq).upper()))
    require(s1>s0,'HCO ne favorise pas le complément inverse : contrôler étiquettes/brins, pas de retournement forcé.')
    r=rawr.reverse_complement(id='reverse_oriented',description='')
    SeqIO.write([f,r],out/'02_consensus'/'reads_oriented.fastq','fastq')
    seq,rows,stats,aln=build(f,r,c)
    stats.update(dict(score_hco_as_is=s0,score_hco_reverse_complement=s1,config_sha256=config_sha(c)))
    tsv(out/'02_consensus'/'base_audit.tsv',rows,list(rows[0]));dump(out/'02_consensus'/'summary.json',stats)
    (out/'02_consensus'/'pair_alignment.txt').write_text(str(aln))
    record=SeqRecord(Seq(seq),id=f"QUERY_{c['sample']}",description='double_covered_conservative_consensus')
    candidate=out/'02_consensus'/'candidate.fasta';SeqIO.write(record,candidate,'fasta')
    translations=translation_rows(seq,c['genetic_code']);tsv(out/'02_consensus'/'translation_six_frames.tsv',translations,list(translations[0]))
    errors=[]
    if stats['high_quality_paired_sites']<c['min_overlap_bp']:errors.append('recouvrement de qualité insuffisant')
    if stats['pair_identity']<c['min_pair_identity']:errors.append('discordance globale LCO/HCO')
    if stats['internal_indel_columns']:errors.append('indel interne : examiner chromatogrammes')
    if stats['n_fraction']>c['max_consensus_n_fraction']:errors.append('trop de N')
    print(json.dumps(stats,indent=2,ensure_ascii=False))
    require(not errors,'; '.join(errors)+'. Pas de consensus validable ; inspecter AB1 et trimming.')
    dump(out/'02_consensus'/'automatic_pass.json',dict(sha256=sha(candidate),config_sha256=config_sha(c),passed=True))
    print('Consensus candidat prêt, NON encore approuvé. Voir base_audit.tsv, chromatogrammes et traduction six cadres.')

if __name__=='__main__':main()
