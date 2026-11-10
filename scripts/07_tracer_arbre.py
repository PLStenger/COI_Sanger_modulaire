from common import *
import argparse
from Bio import Phylo
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def splitkey(names,allnames):
    names=frozenset(names);other=frozenset(allnames)-names
    return min(tuple(sorted(names)),tuple(sorted(other)))


def main():
    c,out=load();validated_query(c,out);p=argparse.ArgumentParser()
    p.add_argument('config');p.add_argument('--mode',choices=['explore','ufboot','standard'],default='explore')
    p.add_argument('--root',default='',help='Accession.version d’un outgroup retenu ; sinon arbre non enraciné.')
    a=p.parse_args();path=out/'06_trees'/f'latest_{a.mode}.json'
    require(path.exists(),'Lancer 06 pour ce mode avant le tracé.')
    run=json.loads(path.read_text());treefile=Path(run['treefile']);require(sha(treefile)==run['tree_sha256'],'Arbre modifié après calcul.')
    manifest=json.loads((out/'05_alignment'/'manifest.json').read_text())
    require(sha(out/'05_alignment'/'analysis_common_sites.fasta')==run['alignment_sha256'],'Alignement modifié après calcul : reprendre 05/06.')
    require(run['alignment_fingerprint']==manifest['fingerprint'],'L’arbre correspond à un ancien alignement ; relancer 06.')
    metadata=json.loads((out/'04_references'/'metadata.json').read_text())
    require(sha(out/'04_references'/'metadata.json')==manifest['metadata_sha256'],'Métadonnées modifiées : relancer 04/05/06.')
    tree=Phylo.read(treefile,'newick');qid='QUERY_'+c['sample']
    names={t.name for t in tree.get_terminals()};require(qid in names,'Query absente de l’arbre.')
    # Les supports annotent des bipartitions, pas une racine biologique.
    supports={}
    for cl in tree.get_nonterminals():
        label=cl.name if cl.name else str(cl.confidence) if cl.confidence is not None else None
        if label:supports[splitkey([t.name for t in cl.get_terminals()],names)]=label
    rooted=False
    if a.root:
        tipname='REF_'+a.root
        require(tipname in names and metadata[tipname]['role']=='outgroup','--root doit désigner une accession outgroup retenue.')
        tree.root_with_outgroup(next(t for t in tree.get_terminals() if t.name==tipname));rooted=True
    else:tree.rooted=False
    for cl in tree.get_nonterminals():
        cl.name=None;cl.confidence=None
    tree.ladderize()
    species=sorted({metadata[n]['organism'] for n in names if n!=qid and metadata[n]['role']=='reference'})
    palette=['#0072B2','#D55E00','#CC79A7','#E69F00','#56B4E9','#8B4513','#5555AA','#999933']
    colors={s:palette[i%len(palette)] for i,s in enumerate(species)}
    if len(species)>len(palette):print('Attention : palette réutilisée ; les noms et accessions restent prioritaires pour identifier les références.')
    labelcolors={}
    for tip in tree.get_terminals():
        if tip.name==qid:color='#008A34';label=f'>>> QUERY {c["sample"]} | {c["field_id"]} <<<'
        else:
            m=metadata[tip.name];color='#777777' if m['role']=='outgroup' else colors[m['organism']]
            label=f"{m['organism']} | {m['accession']}"+(' [outgroup]' if m['role']=='outgroup' else '')
        labelcolors[label]=color
        if tip.name==qid:tip.color=color;tip.width=2.5
    def lab(cl):
        if not cl.is_terminal():return None
        if cl.name==qid:return f'>>> QUERY {c["sample"]} | {c["field_id"]} <<<'
        m=metadata[cl.name];return f"{m['organism']} | {m['accession']}"+(' [outgroup]' if m['role']=='outgroup' else '')
    rootchildren={id(x) for x in tree.root.clades}
    already=set()
    def branch_label(cl):
        if cl.is_terminal() or cl is tree.root:return None
        key=splitkey([t.name for t in cl.get_terminals()],names)
        # Une bipartition scindée par la racine n’est affichée qu’une fois.
        if id(cl) in rootchildren:
            if key in already:return None
            already.add(key)
        return supports.get(key) if a.mode!='explore' else None
    labelmode={'explore':'ML exploratoire — sans bootstrap',
               'ufboot':f"ML — SH-aLRT / UFBoot ({run['n_replicates']} réplicats UFBoot ; -bnni)",
               'standard':f"ML — bootstrap classique ({run['n_replicates']} réplicats)"}[a.mode]
    fig,ax=plt.subplots(figsize=(18,max(7,len(names)*.38+3)))
    Phylo.draw(tree,axes=ax,do_show=False,label_func=lab,label_colors=labelcolors,show_confidence=False,branch_labels=branch_label)
    for t in ax.texts:
        if '>>> QUERY' in t.get_text():t.set_fontweight('bold');t.set_color('#008A34')
    ax.set_title(c['sample']+' | COI | '+labelmode+'\n'+('Enracinement d’affichage sur '+a.root if rooted else 'Non enraciné : la position graphique à gauche n’est pas un ancêtre'))
    ax.set_xlabel('Substitutions estimées par site');ax.set_ylabel('')
    handles=[Line2D([0],[0],color='#008A34',lw=3,label='Query'),Line2D([0],[0],color='#777777',lw=2,label='Outgroup')]
    ax.legend(handles=handles,loc='lower right');fig.tight_layout()
    stem=out/'07_plots'/f"{c['sample']}_{a.mode}"
    for ext in ['pdf','png','svg']:fig.savefig(str(stem)+'.'+ext,dpi=220,bbox_inches='tight')
    plt.close(fig)
    # Les supports originaux restent dans le .treefile IQ-TREE, non réécrit ici.
    dump(str(stem)+'.display.json',dict(source_tree_sha256=sha(treefile),root=a.root or None,
          supports=supports_to_json(supports),note='Supports conservés sur les bipartitions ; le Newick de référence reste le .treefile IQ-TREE.'))
    print('Figures :',str(stem)+'.pdf',str(stem)+'.png')


def supports_to_json(supports):return {'|'.join(k):v for k,v in supports.items()}

if __name__=='__main__':main()
