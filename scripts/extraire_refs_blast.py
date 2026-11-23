import re
import sys
from pathlib import Path

if len(sys.argv) != 4:
    raise SystemExit("Usage : python scripts/extraire_refs_blast.py RAPPORT.txt RID ECHANTILLON")

source = Path(sys.argv[1])
rid = sys.argv[2]
sample = sys.argv[3]
destination = source.parent / (sample + "_references_depuis_BLAST.tsv")

rows = []
seen = set()
in_table = False

for line in source.read_text().splitlines():
    if "Sequences producing significant alignments:" in line:
        in_table = True
        continue
    if not in_table:
        continue
    if line.lstrip().startswith(">"):
        break
    if re.match(r"^\s*(?:Alignments\b|Score\s*=|Query\s+\d+)", line):
        break
    match = re.search(r"\b([A-Z]{1,6}_?\d+\.\d+)\s*$", line)
    if not match:
        continue
    accession = match.group(1)
    if accession in seen:
        continue
    name = re.match(r"^\s*([A-Z][a-z]+)\s+([a-z][a-z-]+)\b", line)
    species = " ".join(name.groups()) if name else ""
    note = "BLAST " + sample + " ; RID " + rid + " ; nom a verifier"
    rows.append("\t".join([accession, species, "reference", note]))
    seen.add(accession)

if not rows:
    raise SystemExit("Aucune accession extraite : verifier le rapport.")

header = "accession\texpected_species\trole\tnotes"
destination.write_text(header + "\n" + "\n".join(rows) + "\n")
print(len(rows), "references extraites dans", destination)
