from common import *
import argparse,re


def current_manifest(c,out):
    validated_query(c,out)
    manifest=json.loads((out/'05_alignment'/'manifest.json').read_text())
    require(manifest['config_sha256']==config_sha(c),'Configuration modifiée : relancer les étapes amont concernées.')
    for key,path in [('alignment_sha256',out/'05_alignment'/'analysis_common_sites.fasta'),
        ('metadata_sha256',out/'04_references'/'metadata.json'),('query_sha256',out/'03_blast'/'query.fasta'),
        ('panel_sha256',ROOT/c['references_file'])]:
        require(manifest[key]==sha(path),f'Données modifiées : {path}. Relancer 04/05.')
    return manifest


def main():
    c,out=load();parser=argparse.ArgumentParser()
    parser.add_argument('config');parser.add_argument('--mode',choices=['explore','ufboot','standard'],default='explore')
    parser.add_argument('--approve',action='store_true');parser.add_argument('--comment',default='')
    args=parser.parse_args();m=current_manifest(c,out)
    if args.approve:
        require(len(args.comment.strip())>=10,'Documenter les références, vouchers, ambiguïtés et l’alignement.')
        dump(out/'alignment_approved.json',dict(sha256=m['fingerprint'],approved=True,comment=args.comment))
    approval(out,'alignment_approved',m['fingerprint'])
    translation=table(out/'05_alignment'/'translation_audit.tsv')
    require(not any(int(x['stop_count'])>0 for x in translation),
            'Codon STOP dans la région COI alignée : revoir le cadre, les séquences et les chromatogrammes avant l’arbre.')
    indels=table(out/'05_alignment'/'indel_audit.tsv')
    require(not any(int(x['internal_gap_columns']) or int(x['insertions_vs_query']) for x in indels),
            'Indel interne ou insertion détectée : revoir l’alignement et envisager un alignement guidé par protéines, sans forcer le bootstrap.')
    program=exe('iqtree3','iqtree2','iqtree')
    version=subprocess.run([program,'-version'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True).stdout
    require(not re.search(r'version\s+1\.',version,re.I),'Utiliser IQ-TREE >=2 pour ces paramètres.')
    alignment=out/'05_alignment'/'analysis_common_sites.fasta'
    opts=[program,'-s',str(alignment),'-st','DNA','-m','MFP','-T',str(c['threads']),'-seed',str(c['seed'])]
    if args.mode=='ufboot':opts+=['-B',str(c['bootstrap_replicates']),'-bnni','-alrt','1000']
    if args.mode=='standard':opts+=['-b',str(c['bootstrap_replicates'])]
    if args.mode=='ufboot':require(c['bootstrap_replicates']>=1000,'UFBoot : au moins 1000 réplicats.')
    signature=hashlib.sha256((m['fingerprint']+str(opts)+c['run_tag']+version).encode()).hexdigest()[:12]
    prefix=out/'06_trees'/f"{c['sample']}_{args.mode}_{c['run_tag']}_{signature}"
    opts+=['-pre',str(prefix)]
    tree=Path(str(prefix)+'.treefile')
    finished=Path(str(prefix)+'.iqtree').exists() and Path(str(prefix)+'.log').exists() and 'Total wall-clock time used:' in Path(str(prefix)+'.log').read_text()
    if not (tree.exists() and finished):command(opts,out/'logs'/f'06_{args.mode}.log')
    else:print('Arbre déjà présent : réutilisation sans relancer les bootstraps.')
    require(tree.exists(),'Arbre absent ; consulter le log IQ-TREE.')
    manifest=dict(treefile=str(tree),mode=args.mode,prefix=str(prefix),
         alignment_fingerprint=m['fingerprint'],n_replicates=c['bootstrap_replicates'] if args.mode!='explore' else 0,
         version=version,command=opts,alignment_sha256=sha(alignment),tree_sha256=sha(tree))
    dump(out/'06_trees'/f'latest_{args.mode}.json',manifest)
    # Le tracé est une étape indépendante : un problème graphique ne supprime pas le calcul.
    print('Arbre obtenu :',tree)
    print(f"Tracer : python scripts/07_tracer_arbre.py {sys.argv[1]} --mode {args.mode}")

if __name__=='__main__':main()
