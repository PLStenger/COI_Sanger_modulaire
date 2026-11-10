from common import *
import io, re, time, urllib.request, urllib.parse
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def read_panel(path):
    with path.open(newline='') as h:
        first=h.readline();h.seek(0)
        sep='\t' if '\t' in first else ';'
        reader=csv.DictReader(h,delimiter=sep)
        require(reader.fieldnames==['accession','expected_species','role','notes'],
                'En-tête : accession;expected_species;role;notes (séparateur ; ou tabulation).')
        rows=[{k:(v or '').strip() for k,v in row.items() if k is not None} for row in reader if any(row.values())]
    require(len(rows)>=4,'Ajouter au moins quatre références dans le panel (dont un outgroup).')
    acc=[r['accession'] for r in rows]
    require(len(acc)==len(set(acc)),'Accessions dupliquées.')
    require(all(re.fullmatch(r'[A-Z]{1,6}_?\d+\.\d+',x) for x in acc),'Utiliser des accessions.version réelles, pas des exemples.')
    require(all(r['role'] in ['reference','outgroup'] for r in rows),'role = reference ou outgroup.')
    require(any(r['role']=='outgroup' for r in rows),'Choisir au moins un outgroup pertinent.')
    return rows


def fetch(accession,path,email):
    if path.exists():
        record=SeqIO.read(path,'genbank')
        require(record.id==accession,f'Cache incohérent : {record.id} != {accession}')
        return record
    params=urllib.parse.urlencode(dict(db='nuccore',id=accession,rettype='gbwithparts',retmode='text',tool='coi_sanger_modulaire',email=email))
    url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'+params
    last=None
    for attempt in range(5):
        try:
            req=urllib.request.Request(url,headers={'User-Agent':'COI-Sanger-modulaire/1.0'})
            with urllib.request.urlopen(req,timeout=90) as res:text=res.read().decode()
            record=SeqIO.read(io.StringIO(text),'genbank')
            require(record.id==accession,f'GenBank a renvoyé {record.id} au lieu de {accession}')
            path.write_text(text);time.sleep(.4);return record
        except Exception as e:last=e;time.sleep(2**attempt)
    raise RuntimeError(f'{accession} : téléchargement impossible : {last}')


def coi_feature(f):
    if f.type != 'CDS':return False
    genes={re.sub('[^a-z0-9]','',g.lower()) for g in f.qualifiers.get('gene',[])}
    product=' '.join(f.qualifiers.get('product',[])).lower()
    return bool(genes & {'cox1','coi','co1','coxi'}) or bool(
        'cytochrome' in product and 'oxidase' in product and
        re.search(r'subunit\s+(?:i|1)\b',product))


def extract_coi(record,query):
    a=PairwiseAligner();a.mode='local';a.match_score=2;a.mismatch_score=-3;a.open_gap_score=-10;a.extend_gap_score=-1
    features=[f for f in record.features if coi_feature(f)]
    candidates=[]
    if features:
        for index,f in enumerate(features):
            try:s=str(f.extract(record.seq)).upper()
            except Exception:continue
            candidates.append((s,f'annotated_CDS_{index+1}',str(f.location)))
    else:
        require(len(record.seq)<=3000,
                f'{record.id} : longue séquence sans CDS COI reconnue ; extraire COI manuellement, pas de mitogénome entier dans MAFFT.')
        candidates=[(str(record.seq).upper(),'unannotated_short_record','full_record')]
    scored=[]
    for s,method,location in candidates:
        for orient,seq in [('forward',s),('reverse_complement',revcomp(s))]:
            try:al=a.align(query,seq)[0]
            except IndexError:continue
            scored.append((al.score,seq,al,method,location,orient))
    require(bool(scored),f'{record.id} : aucun alignement COI détecté.')
    score,s,al,method,location,orientation=max(scored,key=lambda x:x[0])
    qb,rb=al.aligned
    pairs=[(query[i],s[j]) for (i0,i1),(j0,j1) in zip(qb,rb) for i,j in zip(range(i0,i1),range(j0,j1)) if query[i] in 'ACGT' and s[j] in 'ACGT']
    coverage=len(pairs)/sum(b in 'ACGT' for b in query)
    start=max(0,int(rb[0,0])-15);end=min(len(s),int(rb[-1,1])+15)
    return s[start:end],dict(extraction=method,feature_location=location,orientation=orientation,
                           local_comparable_bp=len(pairs),local_coverage=round(coverage,5),
                           local_identity=round(sum(x==y for x,y in pairs)/len(pairs),5) if pairs else 0,
                           extracted_bp=end-start,alignment_score=score)


def main():
    c,out=load();query_path=validated_query(c,out)
    require(bool(c.get('email','').strip()),'Renseigner email dans la configuration pour les appels NCBI.')
    panel=ROOT/c['references_file'];require(panel.exists(),f'Panel absent : {panel}')
    rows=read_panel(panel);query=str(SeqIO.read(query_path,'fasta').seq).upper()
    records=[SeqRecord(Seq(query),id='QUERY_'+c['sample'],description='')];audit=[];metadata={}
    for row in rows:
        acc=row['accession'];rec=fetch(acc,out/'04_references'/f'{acc}.gb',c['email'])
        source=next((f for f in rec.features if f.type=='source'),None)
        quals=source.qualifiers if source else {}
        organism=rec.annotations.get('organism','')
        fragment,ex=extract_coi(rec,query)
        retained=ex['local_comparable_bp']>=c['min_reference_comparable_bp'] and ex['local_coverage']>=c['min_reference_coverage']
        # L'identifiant reste indépendant de l'orthographe du nom d'espèce.
        name='REF_'+acc
        species_match=(row['expected_species'].casefold()==organism.casefold()) if row['expected_species'] else None
        item=dict(id=name,accession=acc,expected_species=row['expected_species'],organism=organism,
                  species_match=species_match,role=row['role'],notes=row['notes'],definition=rec.description,
                  voucher=';'.join(quals.get('specimen_voucher',[])),bio_material=';'.join(quals.get('bio_material',[])),isolate=';'.join(quals.get('isolate',[])),
                  country=';'.join(quals.get('geo_loc_name',quals.get('country',[]))),
                  record_sha256=sha(out/'04_references'/f'{acc}.gb'),retained=retained,**ex)
        audit.append(item)
        if retained:
            records.append(SeqRecord(Seq(fragment),id=name,description=''));metadata[name]=item
        print(acc,organism,'retenue' if retained else 'ECARTEE',ex['orientation'],ex['local_comparable_bp'])
    tsv(out/'04_references'/'audit.tsv',audit,list(audit[0]))
    dump(out/'04_references'/'metadata.json',metadata)
    target=out/'04_references'/'oriented_COI.fasta';SeqIO.write(records,target,'fasta')
    dump(out/'04_references'/'provenance.json',dict(query_sha256=sha(query_path),panel_sha256=sha(panel),
         fasta_sha256=sha(target),metadata_sha256=sha(out/'04_references'/'metadata.json'),config_sha256=config_sha(c)))
    require(len(records)>=5,'Moins de quatre références conservées ; enrichir/corriger le panel.')
    require(any(x['role']=='outgroup' for x in metadata.values()),'Aucun outgroup conservé. Choisir un groupe externe plus proche et couvrant COI.')
    print('Audit à examiner : organism, expected_species, voucher, couverture, exclusions. Aucun nom déclaré ne remplace le nom NCBI.')

if __name__=='__main__':main()
