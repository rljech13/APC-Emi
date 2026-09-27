"""Assemble a one-sequence-per-species Emi2/FBXO43 orthologue set for alignment."""
import sys, re, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from seqs import read_fa

HERE = os.path.dirname(os.path.abspath(__file__))
LOCAL = os.path.join(os.path.dirname(__file__), "..", "Fasta", "Emi2") + os.sep


def read_ncbi(path):
    """Return list of (accession, species, sequence)."""
    recs, acc, sp, buf = [], None, None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith('>'):
            if acc:
                recs.append((acc, sp, ''.join(buf)))
            acc = line[1:].split()[0]
            m = re.search(r'\[([^\]]+)\]', line)
            sp = m.group(1) if m else 'unknown'
            sp = sp.split('(')[0].strip()
            buf = []
        elif line:
            buf.append(line.strip())
    if acc:
        recs.append((acc, sp, ''.join(buf)))
    return recs


records = []
records += read_ncbi(os.path.join(HERE, 'squamata_fbxo43_raw.fasta'))
records += read_ncbi(os.path.join(HERE, 'outgroups.fasta'))

# keep the longest isoform per species
best = {}
for acc, sp, s in records:
    if len(s) < 350:          # drop fragments / mis-annotations
        continue
    if sp not in best or len(s) > len(best[sp][1]):
        best[sp] = (acc, s)

# reviewed UniProt entries
uni = read_fa(os.path.join(HERE, 'Emi2_reviewed_all.fasta'))
uni_map = {
    'sp|Q4G163|FBX43_HUMAN': 'Homo sapiens',
    'sp|Q8CDI2|FBX43_MOUSE': 'Mus musculus',
    'sp|Q8AXF4|FBX43_XENLA': 'Xenopus laevis',
    'sp|Q66H04|FBX43_RAT': 'Rattus norvegicus',
}
for k, sp in uni_map.items():
    if k in uni:
        best[sp] = (k.split('|')[1], uni[k])

# the study's own sequences
for fn, sp in [('Emi2_Dval_AAseq.fasta', 'Darevskia valentini'),
               ('Emi2_Dnai_AAseq.fasta', 'Darevskia raddei nairensis'),
               ('Emi2_Pmur_AAseq.fasta', 'Podarcis muralis_study'),
               ('Emi2_Emac_AAseq.fasta', 'Eublepharis macularius_study')]:
    d = read_fa(LOCAL + fn)
    key = list(d)[0]
    best[sp] = (key, d[key])

# emit
sp_names = sorted(best)
with open(os.path.join(HERE, 'emi2_orthologues.fasta'), 'w') as fh:
    for sp in sp_names:
        acc, s = best[sp]
        tag = sp.replace(' ', '_')
        fh.write(f'>{tag}|{acc}|{len(s)}\n')
        for i in range(0, len(s), 60):
            fh.write(s[i:i + 60] + '\n')

print(f'{len(sp_names)} species written')
for sp in sp_names:
    print(f'  {sp:42s} {best[sp][0]:18s} {len(best[sp][1])}')
