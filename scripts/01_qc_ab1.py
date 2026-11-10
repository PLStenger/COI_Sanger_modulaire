from common import *
import numpy as np
import Bio, platform, sys
from Bio import SeqIO
from Bio.Seq import Seq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def good_interval(q, window, threshold):
    require(len(q)>=window, 'Lecture plus courte que la fenêtre qualité.')
    means=np.convolve(q,np.ones(window)/window,mode='valid')
    good=means>=threshold
    best=(0,0); start=None
    for i,value in enumerate(list(good)+[False]):
        if value and start is None:start=i
        if not value and start is not None:
            end=i-1+window
            if end-start > best[1]-best[0]:best=(start,end)
            start=None
    require(best[1]>best[0], 'Aucune plage ne passe le filtre qualité.')
    return best


def primer_bounds(s, primers, search, mismatches):
    left,right=0,len(s);events=[]
    variants=set(primers)|{revcomp(p) for p in primers}
    for p in variants:
        n=len(p)
        for start in range(min(search,max(0,len(s)-n+1))):
            if sum(a!=b for a,b in zip(s[start:start+n],p))<=mismatches:
                if start+n>left:left=start+n;events.append('left:'+p)
        for start in range(max(0,len(s)-search-n),len(s)-n+1):
            if sum(a!=b for a,b in zip(s[start:start+n],p))<=mismatches:
                if start<right:right=start;events.append('right:'+p)
    require(right>left, 'Amorces : intervalle vide, vérifier les lectures.')
    return left,right,';'.join(events)


def main():
    c,out=load()
    require(c.get('mapping_confirmed') is True,
            'Vérifier DI/D1 et la correspondance des fichiers, puis mapping_confirmed=true dans la configuration.')
    summaries=[]; source_hash={}
    for role,key in [('forward','forward_ab1'),('reverse','reverse_ab1')]:
        path=Path(c['raw_dir']).expanduser()/c[key]
        require(path.exists(),f'Fichier absent : {path}')
        rec=SeqIO.read(path,'abi'); seq=str(rec.seq).upper()
        q=list(rec.letter_annotations.get('phred_quality',[]))
        require(len(q)==len(seq),f'Qualités PHRED absentes : {path}')
        source_hash[role]=sha(path)
        manual=c.get('raw_trim_1based',{}).get(role)
        a,b=(int(manual[0])-1,int(manual[1])) if manual else (0,len(seq))
        require(0<=a<b<=len(seq),'Bornes manuelles AB1 incorrectes.')
        x,y=good_interval(q[a:b],c['quality_window'],c['quality_window_mean'])
        start,end=a+x,a+y
        pstart,pend,events=primer_bounds(seq[start:end],c['primers'],c['primer_end_search_bp'],c['primer_max_mismatches'])
        old=start;start=old+pstart;end=old+pend
        trimmed=rec[start:end];trimmed.id=role;trimmed.description=''
        require(len(trimmed)>=c['min_read_bp'],f'{role} : lecture trop courte après trimming.')
        SeqIO.write(trimmed,out/'01_qc'/f'{role}.fastq','fastq')
        SeqIO.write(trimmed,out/'01_qc'/f'{role}.fasta','fasta')
        raw=rec.annotations.get('abif_raw',{})
        loc=raw.get('PLOC2',[]);order=raw.get('FWO_1',b'GATC')
        if isinstance(order,bytes):order=order.decode(errors='replace')
        traces={base:raw.get(f'DATA{9+i}',[]) for i,base in enumerate(order)}
        colors={'A':'#16853B','C':'#1262D6','G':'#222222','T':'#CE2424'}
        with PdfPages(out/'01_qc'/f'{role}_chromatogrammes.pdf') as pdf:
            for lo in range(0,len(seq),100):
                hi=min(len(seq),lo+100)
                fig,(ax,aq)=plt.subplots(2,1,figsize=(17,5),gridspec_kw={'height_ratios':[3,1]})
                if len(loc)>=hi and all(len(v)>0 for v in traces.values()):
                    t0=max(0,int(loc[lo])-10);t1=min(min(map(len,traces.values())),int(loc[hi-1])+12)
                    for base,values in traces.items():ax.plot(range(t0,t1),values[t0:t1],color=colors.get(base,'grey'),lw=0.8,label=base)
                    for k in range(lo,hi):
                        ax.text(loc[k],1.01,seq[k],transform=ax.get_xaxis_transform(),ha='center',fontsize=7,color=colors.get(seq[k],'grey'))
                    ax.set_xlim(t0,t1);ax.legend(loc='upper right',ncol=4)
                else:ax.text(.1,.5,'Traces DATA9–12 / PLOC2 indisponibles : ouvrir le .ab1 dans un lecteur dédié.',transform=ax.transAxes)
                aq.bar(range(lo+1,hi+1),q[lo:hi],width=1,color='grey')
                aq.axhline(c['min_base_q'],color='red',lw=.8)
                aq.axvspan(max(lo+1,start+1),min(hi,end),color='green',alpha=.12) if max(lo+1,start+1)<=min(hi,end) else None
                aq.set_xlim(lo+.5,hi+.5);aq.set_ylabel('PHRED');aq.set_xlabel('Position dans la lecture brute (1-based)')
                fig.suptitle(f"{c['sample']} | {path.name} | plage retenue : {start+1}–{end}")
                fig.tight_layout();pdf.savefig(fig);plt.close(fig)
        tsv(out/'01_qc'/f'{role}_positions.tsv',[
            dict(raw_position=i+1,base=base,phred=q[i],retained=start<=i<end)
            for i,base in enumerate(seq)],['raw_position','base','phred','retained'])
        summaries.append(dict(read=role,file=path.name,raw_bp=len(seq),start_1based=start+1,end_1based=end,trimmed_bp=end-start,
                              mean_phred=round(float(np.mean(q[start:end])),2),primers_removed=events))
    tsv(out/'01_qc'/'summary.tsv',summaries,list(summaries[0]))
    dump(out/'01_qc'/'provenance.json',dict(python_version=sys.version,biopython_version=Bio.__version__,platform=platform.platform(),config_sha256=config_sha(c),raw_sha256=source_hash,
        outputs={r:sha(out/'01_qc'/f'{r}.fastq') for r in ['forward','reverse']}))
    print('QC terminé. Examiner summary.tsv et les deux PDF de chromatogrammes AVANT validation du consensus.')

if __name__=='__main__':main()
