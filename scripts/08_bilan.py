from common import *
from datetime import datetime, timezone
from Bio import SeqIO


def main():
    c,out=load();qpath=validated_query(c,out)
    qc=json.loads((out/'02_consensus'/'summary.json').read_text());seq=str(SeqIO.read(qpath,'fasta').seq)
    dfile=out/'05_alignment'/'query_distances.tsv';require(dfile.exists(),'Lancer 05 avant le bilan.')
    rows=table(dfile);common=sorted([r for r in rows if r['dataset']=='common_sites'],key=lambda r:float(r['p_distance']))
    manifest=json.loads((out/'05_alignment'/'manifest.json').read_text())
    require(manifest['query_sha256']==sha(qpath),'Distances provenant d’une ancienne query : reprendre 04/05.')
    require(manifest['panel_sha256']==sha(ROOT/c['references_file']),'Panel modifié depuis les distances : reprendre 04/05.')
    lines=[f"# {c['sample']} — bilan COI provisoire",'',
        f"Généré : {datetime.now(timezone.utc).isoformat()}",
        f"Identifiant de terrain : {c['field_id']}",
        f"Taxon présumé fourni : {c['presumed_taxon']} (non imposé à l’analyse)",'',
        '## Métadonnées fournies', '',
        *[f'- {k} : {v}' for k,v in c.get('metadata',{}).items()], '',
        '## Consensus', '',
        f"- Longueur : {len(seq)} bases ; N : {seq.count('N')}.",
        f"- Recouvrement HQ LCO/HCO : {qc['high_quality_paired_sites']} sites ; accord : {100*qc['pair_identity']:.3f} %.",
        f"- SHA256 query : {sha(qpath)}",'',
        '## Distances sur sites communs','',
        f"Toutes les références ci-dessous sont comparées sur les mêmes {manifest['common_columns']} sites.",'',
        '| Accession | Taxon du dépôt | Rôle | Différences | Distance p (%) |',
        '|---|---|---|---|---:|']
    for r in common:
        lines.append(f"| {r['accession']} | {r['species']} | {r['role']} | {r['differences']}/{r['comparable_bp']} | {100*float(r['p_distance']):.4f} |")
    lines+=['','## Analyses phylogénétiques','']
    for mode in ['explore','ufboot','standard']:
        p=out/'06_trees'/f'latest_{mode}.json'
        if p.exists():
            run=json.loads(p.read_text())
            status='alignement courant' if run['alignment_fingerprint']==manifest['fingerprint'] else 'ANCIEN ALIGNEMENT : ne pas interpréter'
            lines.append(f"- {mode} : {run['treefile']} ; {run['n_replicates']} réplicats ; {status}.")
    lines+=['','## Conclusion à rédiger après revue','',
        'Aucune identification spécifique automatique. Documenter BLAST, couverture, ex aequo interspécifiques, vouchers et support du clade pertinent.',
        'En cas de discordance ou de complexe non discriminé par COI : limiter le niveau taxonomique et proposer une vérification indépendante.',
        'Pour H1 présumé Bactrocera : ne pas assimiler haplotype COI et identification certaine du complexe dorsalis/carambolae.',
        '', 'Les fichiers bruts, audits, commandes et empreintes sont à conserver avec le compte rendu signé.']
    target=out/'BILAN_PROVISOIRE.md';target.write_text('\n'.join(lines)+'\n');print(target)

if __name__=='__main__':main()
