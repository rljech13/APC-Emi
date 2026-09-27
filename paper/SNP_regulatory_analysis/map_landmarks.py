"""Project literature-defined Emi2 regulatory landmarks onto Darevskia numbering
via the MAFFT L-INS-i orthologue alignment, and tabulate cross-vertebrate
residue distributions at the ten Dval/Dnai polymorphic positions."""
import os, sys, json
from collections import Counter
from Bio import AlignIO

HERE = os.path.dirname(os.path.abspath(__file__))
ALN = AlignIO.read(os.path.join(HERE, 'emi2_orthologues_linsi.aln.fasta'), 'fasta')

seqs = {r.id.split('|')[0]: str(r.seq) for r in ALN}
REF = 'Darevskia_valentini'


def col_index(tag):
    """ungapped residue number (1-based) -> alignment column, and inverse"""
    s = seqs[tag]
    res2col, col2res = {}, {}
    n = 0
    for c, ch in enumerate(s):
        if ch != '-':
            n += 1
            res2col[n] = c
            col2res[c] = n
    return res2col, col2res


IDX = {t: col_index(t) for t in seqs}


def project(src_tag, src_pos, dst_tag=REF):
    """map residue src_pos in src_tag onto dst_tag; returns (pos_or_None, src_aa, dst_aa)"""
    col = IDX[src_tag][0].get(src_pos)
    if col is None:
        return None, None, None
    src_aa = seqs[src_tag][col]
    dst_aa = seqs[dst_tag][col]
    dst_pos = IDX[dst_tag][1].get(col)
    return dst_pos, src_aa, dst_aa


# --------------------------------------------------------------------------
# Literature landmarks. (species tag, position, label, source)
# --------------------------------------------------------------------------
LANDMARKS = [
    # --- beta-TrCP phosphodegron, Plk1-phosphorylated ---
    ('Xenopus_laevis',  33, 'DSGx3S degron pSer-1 (Plk1)', 'Hansen 2006; Rauh 2005'),
    ('Xenopus_laevis',  38, 'DSGx3S degron pSer-2 (Plk1)', 'Hansen 2006; Rauh 2005'),
    # --- Plk1 polo-box docking phosphothreonines (Jia 2015, mouse numbering) ---
    ('Mus_musculus',   152, 'Plk1-PBD docking pThr-1',     'Jia 2015 (5DMV)'),
    ('Mus_musculus',   169, 'Plk1-PBD Tyr-pocket anchor Phe', 'Jia 2015 (5DMV)'),
    ('Mus_musculus',   176, 'Plk1-PBD docking pThr-2',     'Jia 2015 (5DMV/5DMZ)'),
    ('Mus_musculus',   157, 'Plk1-PBD acidic-patch contact D', 'Jia 2015'),
    ('Mus_musculus',   158, 'Plk1-PBD acidic-patch contact V', 'Jia 2015'),
    ('Mus_musculus',   159, 'Plk1-PBD acidic-patch contact V', 'Jia 2015'),
    # --- CaMKII / Plk1 priming site in Xenopus ---
    ('Xenopus_laevis', 170, 'Xenopus Thr170 (Plk1-PBD site)', 'Isoda 2011'),
    ('Xenopus_laevis', 195, 'Xenopus Thr195 (CaMKII / Plk1-PBD site)', 'Rauh 2005; Hansen 2006; Isoda 2011'),
    # --- p90rsk / Mos-MAPK sites ---
    ('Xenopus_laevis', 335, 'p90rsk Ser335', 'Inoue 2007'),
    ('Xenopus_laevis', 336, 'p90rsk Thr336', 'Nishiyama 2007; Inoue 2007'),
    ('Xenopus_laevis', 342, 'p90rsk Ser342', 'Nishiyama 2007'),
    ('Xenopus_laevis', 344, 'p90rsk Ser344', 'Nishiyama 2007'),
    # --- C-terminal Cdk1 sites ---
    ('Xenopus_laevis', 545, 'Cdk1 Thr545', 'Wu 2007; Isoda 2011'),
    ('Xenopus_laevis', 551, 'Cdk1 Thr551', 'Wu 2007; Isoda 2011'),
    # --- APC/C-binding elements ---
    ('Xenopus_laevis', 529, 'D-box Arg', 'Ohe 2010'),
    ('Xenopus_laevis', 532, 'D-box Leu', 'Ohe 2010'),
    ('Xenopus_laevis', 583, 'ZBR essential Cys583', 'Ohe 2010'),
    ('Xenopus_laevis', 629, 'RL tail start', 'Ohe 2010'),
    ('Xenopus_laevis', 651, 'RL tail end (Leu)', 'Ohe 2010'),
]

print('=' * 100)
print('LANDMARK PROJECTION ONTO DAREVSKIA (Dval) NUMBERING')
print('=' * 100)
print(f'{"source":18s} {"pos":>5s} {"aa":>3s}  ->  {"Dval pos":>8s} {"Dval aa":>7s}   label')
proj = {}
for tag, pos, label, src in LANDMARKS:
    dp, sa, da = project(tag, pos)
    proj[label] = (dp, da)
    flag = '' if sa == da else '   [not identical]'
    print(f'{tag[:18]:18s} {pos:5d} {sa!s:>3s}  ->  {str(dp):>8s} {da!s:>7s}   {label}{flag}')

# --------------------------------------------------------------------------
# Cross-vertebrate residue distribution at the ten polymorphic positions
# --------------------------------------------------------------------------
SUBS = [(76, 'Q', 'H'), (92, 'V', 'I'), (102, 'A', 'T'), (171, 'K', 'R'),
        (219, 'E', 'D'), (391, 'I', 'V'), (395, 'D', 'N'), (425, 'L', 'F'),
        (549, 'F', 'I'), (590, 'V', 'M')]

EXCLUDE = {'Darevskia_valentini', 'Darevskia_raddei', 'Podarcis_muralis_study',
           'Eublepharis_macularius_study'}
# handle the truncated tag for D. raddei nairensis
EXCLUDE |= {t for t in seqs if t.startswith('Darevskia')}

SQUAMATES_ONLY = {t for t in seqs if t not in EXCLUDE and t not in
                  {'Homo_sapiens', 'Mus_musculus', 'Rattus_norvegicus',
                   'Xenopus_laevis', 'Gallus_gallus', 'Danio_rerio',
                   'Chelonia_mydas', 'Alligator_mississippiensis'}}
OUTGROUP = {'Homo_sapiens', 'Mus_musculus', 'Rattus_norvegicus', 'Xenopus_laevis',
            'Gallus_gallus', 'Danio_rerio', 'Chelonia_mydas',
            'Alligator_mississippiensis'}

print('\n' + '=' * 100)
print('CROSS-VERTEBRATE RESIDUE DISTRIBUTION AT THE TEN POLYMORPHIC POSITIONS')
print(f'(n squamates = {len(SQUAMATES_ONLY)}, n non-squamate outgroups = {len(OUTGROUP)})')
print('=' * 100)

REPORT = {}
for pos, aval, anai in SUBS:
    col = IDX[REF][0][pos]
    sq = Counter(seqs[t][col] for t in SQUAMATES_ONLY)
    og = Counter(seqs[t][col] for t in OUTGROUP)
    og_detail = {t.split('_')[0][:4] + '.' + t.split('_')[-1][:3]: seqs[t][col] for t in sorted(OUTGROUP)}
    tot = sq + og
    n = sum(tot.values())
    modal, modaln = tot.most_common(1)[0]
    REPORT[pos] = dict(dval=aval, dnai=anai, squamate=dict(sq), outgroup=dict(og),
                       outgroup_detail=og_detail,
                       modal=modal, modal_frac=round(modaln / n, 3), n=n)
    print(f'\n  {aval}{pos}{anai}   (Dval={aval}, Dnai={anai})')
    print(f'    squamates (n={sum(sq.values())}): ' +
          ', '.join(f'{a}:{c}' for a, c in sq.most_common()))
    print(f'    outgroups (n={sum(og.values())}): ' +
          ', '.join(f'{a}:{c}' for a, c in og.most_common()))
    print(f'    per-outgroup: ' + ', '.join(f'{k}={v}' for k, v in og_detail.items()))
    print(f'    modal residue over all {n}: {modal} ({REPORT[pos]["modal_frac"]:.0%})')

json.dump({'landmarks': {k: v[0] for k, v in proj.items()},
           'positions': REPORT}, open(os.path.join(HERE, 'mapping_results.json'), 'w'),
          indent=2)
print('\nwrote mapping_results.json')
