from common import *
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def main():
    c,out=load()
    p=argparse.ArgumentParser(description='Validation EXPLICITE après inspection des AB1 ; pas une correction automatique.')
    p.add_argument('config');p.add_argument('--strand',choices=['+','-'],required=True)
    p.add_argument('--offset',type=int,choices=[0,1,2],required=True)
    p.add_argument('--approve',action='store_true');p.add_argument('--comment',required=True)
    a=p.parse_args();require(a.approve,'Ajouter --approve après revue manuelle, pas avant.')
    path=out/'02_consensus'/'candidate.fasta';auto=out/'02_consensus'/'automatic_pass.json'
    require(auto.exists(),'02 doit réussir avant approbation.')
    audit=json.loads(auto.read_text())
    require(audit['sha256']==sha(path) and audit['config_sha256']==config_sha(c),'Consensus/configuration modifiés : relancer 01 et 02.')
    seq=str(SeqIO.read(path,'fasta').seq).upper()
    chosen=next(r for r in translation_rows(seq,c['genetic_code']) if r['strand']==a.strand and r['offset_0based']==a.offset)
    require(chosen['stop_count']==0,'Cadre avec codon STOP : arrêter et inspecter, ne pas choisir arbitrairement un autre cadre.')
    require(len(a.comment.strip())>=10,'Documenter la revue des chromatogrammes, conflits et cadre de lecture.')
    if a.strand=='-':seq=revcomp(seq)
    target=out/'03_blast'/'query.fasta'
    SeqIO.write(SeqRecord(Seq(seq),id='QUERY_'+c['sample'],description='manually_reviewed_consensus'),target,'fasta')
    dump(out/'consensus_approved.json',dict(approved=True,sha256=sha(target),candidate_sha256=sha(path),
        config_sha256=config_sha(c),coding_offset_0based=a.offset,strand_from_candidate=a.strand,
        comment=a.comment))
    (out/'03_blast'/'BLAST_instructions.txt').write_text(
        'https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch\n'
        'Coller query.fasta. Première recherche sans restriction taxonomique.\n'
        'Choisir megablast pour proches ; blastn si aucun hit convaincant.\n'
        'Exporter les résultats, noter date/base/paramètres/RID.\n'
        'Comparer identité, couverture, gaps, ex aequo et vouchers, pas seulement le premier hit.\n')
    print(f'Query validée pour BLAST : {target}\nCadre conservé pour contrôles COI, aucune identification imposée.')

if __name__=='__main__':main()
