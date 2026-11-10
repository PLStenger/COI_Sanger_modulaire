from pathlib import Path
import csv, hashlib, json, shutil, subprocess, sys

ROOT = Path(__file__).resolve().parents[1]

def require(ok, message):
    if not ok:
        raise SystemExit('ERREUR : ' + message)

def load():
    require(len(sys.argv) >= 2, 'Indiquer une configuration : configs/H1.json')
    path = Path(sys.argv[1]).expanduser().resolve()
    c = json.loads(path.read_text())
    require(c['sample'] in {'H1','A2','B2','C2','E2'}, 'Code non prévu.')
    out = ROOT / 'analyses' / c['sample']
    for d in ['01_qc','02_consensus','03_blast','04_references','05_alignment','06_trees','07_plots','logs']:
        (out/d).mkdir(parents=True, exist_ok=True)
    return c, out

def dump(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, ensure_ascii=False)+'\n')

def tsv(path, rows, fields):
    with Path(path).open('w', newline='') as h:
        w=csv.DictWriter(h,fieldnames=fields,delimiter='\t',extrasaction='ignore')
        w.writeheader(); w.writerows(rows)

def table(path):
    with Path(path).open(newline='') as h:
        return list(csv.DictReader(h, delimiter='\t'))

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def config_sha(c):
    return hashlib.sha256(json.dumps(c,sort_keys=True).encode()).hexdigest()

def exe(*names):
    for n in names:
        p=shutil.which(n)
        if p:return p
    raise SystemExit('ERREUR : exécutable absent : ' + ' / '.join(names))

def command(cmd, log, stdout=None):
    print('+', ' '.join(map(str,cmd)), flush=True)
    with Path(log).open('a') as err:
        err.write('\nCOMMAND '+ ' '.join(map(str,cmd))+'\n');err.flush()
        if stdout:
            with Path(stdout).open('w') as out:
                subprocess.run(list(map(str,cmd)),stdout=out,stderr=err,check=True)
        else:
            subprocess.run(list(map(str,cmd)),stdout=err,stderr=err,check=True)

def approval(out, name, fingerprint):
    p=out/(name+'.json')
    require(p.exists(), f'Validation manuelle absente : {p}. Voir le guide.')
    a=json.loads(p.read_text())
    require(a.get('sha256')==fingerprint and a.get('approved') is True,
            f'Validation périmée ou refusée : {p}. Revérifier les résultats.')
    require(bool(a.get('comment','').strip()), 'Ajouter un commentaire de validation.')

def revcomp(s):
    from Bio.Seq import Seq
    return str(Seq(s).reverse_complement())

def translation_rows(seq, table_id):
    from Bio.Seq import Seq
    rows=[]
    for strand,s in [('+',seq),('-',revcomp(seq))]:
        for off in range(3):
            frag=s[off:];frag=frag[:len(frag)//3*3]
            aa=str(Seq(frag).translate(table=table_id))
            rows.append(dict(strand=strand,offset_0based=off,aa_length=len(aa),
                             stop_count=aa.count('*'),unknown_aa=aa.count('X'),protein=aa))
    return rows


def validated_query(c,out):
    p=out/'03_blast'/'query.fasta'
    require(p.exists(),'Query absente : lancer 01, 02 puis 03.')
    approval(out,'consensus_approved',sha(p))
    approved=json.loads((out/'consensus_approved.json').read_text())
    require(approved['config_sha256']==config_sha(c),
            'Configuration modifiée après revue du consensus : relancer 01, 02 et 03.')
    require(approved['candidate_sha256']==sha(out/'02_consensus'/'candidate.fasta'),
            'Consensus candidat modifié : relancer 02 et 03.')
    provenance=json.loads((out/'01_qc'/'provenance.json').read_text())
    for role,key in [('forward','forward_ab1'),('reverse','reverse_ab1')]:
        require(sha(Path(c['raw_dir'])/c[key])==provenance['raw_sha256'][role],
                'AB1 source modifié après consensus : reprendre depuis 01.')
    return p
