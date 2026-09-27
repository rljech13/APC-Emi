"""Per-substitution local context: domain boundaries, distance to nearest mapped
regulatory element, and kinase/degron consensus motifs created or destroyed."""
import os, sys, re, math
from Bio import AlignIO

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from seqs import read_fa

B = os.path.join(os.path.dirname(__file__), "..", "Fasta", "Emi2") + os.sep
dval = read_fa(B + 'Emi2_Dval_AAseq.fasta')['Emi2_Dval_AAseq']
dnai = read_fa(B + 'Emi2_Dnai_AAseq.fasta')['Emi2_Dnai_AAseq']

ALN = AlignIO.read(os.path.join(HERE, 'emi2_orthologues_linsi.aln.fasta'), 'fasta')
S = {r.id.split('|')[0]: str(r.seq) for r in ALN}


def index(tag):
    r2c, c2r, n = {}, {}, 0
    for c, ch in enumerate(S[tag]):
        if ch != '-':
            n += 1
            r2c[n] = c
            c2r[c] = n
    return r2c, c2r


I = {t: index(t) for t in S}
REF = 'Darevskia_valentini'


def proj(tag, pos):
    col = I[tag][0].get(pos)
    if col is None:
        return None
    # walk to nearest column that is not a gap in REF
    for d in range(0, 12):
        for c in (col - d, col + d):
            if 0 <= c < len(S[REF]) and S[REF][c] != '-':
                return I[REF][1][c]
    return None


print('=' * 92)
print('DOMAIN BOUNDARIES PROJECTED ONTO DAREVSKIA NUMBERING (UniProt annotations)')
print('=' * 92)
DOMS = [('Homo_sapiens', 490, 547, 'F-box (human 490-547)'),
        ('Xenopus_laevis', 424, 499, 'F-box (Xenopus 424-499)'),
        ('Mus_musculus', 423, 480, 'F-box (mouse 423-480)'),
        ('Homo_sapiens', 636, 684, 'ZBR (human 636-684)'),
        ('Xenopus_laevis', 579, 627, 'ZBR (Xenopus 579-627)'),
        ('Homo_sapiens', 320, 426, 'UniProt "Disordered" (human 320-426)')]
for tag, s, e, lab in DOMS:
    print(f'  {lab:38s} -> Dval {proj(tag,s)}-{proj(tag,e)}')

# ------------------------------------------------------------------ elements
# label: (start, end) in Dval numbering, established above
ELEMENTS = {
    'beta-TrCP DSGx3S phosphodegron (Plk1 sites S41/S46)': (40, 46),
    'Plk1-PBD docking module (pT182 / F199 / pT206)': (176, 207),
    'PP2A-B56 LxxIxE motif L358-x-x-L361-x-E363 + p90rsk pS359/pT360': (358, 363),
    'p90rsk cluster pS366 / pS368': (366, 368),
    'F-box domain': (452, 537),
    'D-box R548-x-x-L551': (548, 551),
    'Cdk1 site T574-P575': (574, 575),
    'ZBR (8 Zn-coordinating Cys 607-648)': (603, 651),
    'RL tail (APC/C docking)': (662, 675),
}

SUBS = [(76, 'Q', 'H'), (92, 'V', 'I'), (102, 'A', 'T'), (171, 'K', 'R'),
        (219, 'E', 'D'), (391, 'I', 'V'), (395, 'D', 'N'), (425, 'L', 'F'),
        (549, 'F', 'I'), (590, 'V', 'M')]

print('\n' + '=' * 92)
print('DISTANCE FROM EACH SUBSTITUTION TO THE NEAREST MAPPED REGULATORY ELEMENT')
print('=' * 92)
for p, a, b in SUBS:
    rows = []
    for lab, (s, e) in ELEMENTS.items():
        d = 0 if s <= p <= e else (s - p if p < s else p - e)
        rows.append((d, lab, s, e))
    rows.sort()
    d, lab, s, e = rows[0]
    inside = 'INSIDE' if d == 0 else f'{d} aa away'
    print(f'  {a}{p}{b:<2s}  {inside:>12s}  {lab} [{s}-{e}]')
    for d2, lab2, s2, e2 in rows[1:3]:
        print(f'{"":18s}  {d2:>9d} aa  {lab2} [{s2}-{e2}]')

# ------------------------------------------------------------------ motifs
print('\n' + '=' * 92)
print('LOCAL SEQUENCE CONTEXT AND KINASE / DEGRON CONSENSUS SCAN (+/- 8 residues)')
print('=' * 92)

MOTIFS = [
    ('Cdk1/MAPK proline-directed S/T-P', r'[ST]P'),
    ('full Cdk1 consensus S/T-P-x-K/R', r'[ST]P.[KR]'),
    ('basophilic R/K-x-x-S/T (PKA/CaMKII/RSK)', r'[RK]..[ST]'),
    ('CaMKII R-x-x-S/T', r'R..[ST]'),
    ('polo-box docking S-pS/pT-P', r'S[ST]P'),
    ('polo-box docking S-pS/pT (minimal)', r'S[ST]'),
    ('Plk1 kinase consensus D/E-x-S/T', r'[DE].[ST]'),
    ('beta-TrCP degron DSGxxxS', r'DSG...S'),
    ('beta-TrCP degron D-S-G (core)', r'DSG'),
    ('PP2A-B56 LxxIxE', r'[LMFI]..[ILV].E'),
    ('D-box R-x-x-L', r'R..L'),
    ('CK1 consensus pS-x-x-S (acidic n-3)', r'[DES]..[ST]'),
]

HYDRO = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4,
             H=-3.2, I=4.5, L=3.8, K=-3.9, M=1.9, F=2.8, P=-1.6, S=-0.8,
             T=-0.7, W=-0.9, Y=-1.3, V=4.2)
CHARGE = dict(D=-1, E=-1, K=1, R=1, H=0.1)
VOL = dict(A=88.6, R=173.4, N=114.1, D=111.1, C=108.5, Q=143.8, E=138.4,
           G=60.1, H=153.2, I=166.7, L=166.7, K=168.6, M=162.9, F=189.9,
           P=112.7, S=89.0, T=116.1, W=227.8, Y=193.6, V=140.0)

for p, a, b in SUBS:
    lo, hi = p - 9, p + 8
    wv, wn = dval[lo:hi], dnai[lo:hi]
    print(f'\n  --- {a}{p}{b} ---')
    print(f'    Dval {lo+1}-{hi}: {wv[:8]}[{wv[8]}]{wv[9:]}')
    print(f'    Dnai {lo+1}-{hi}: {wn[:8]}[{wn[8]}]{wn[9:]}')
    dh = HYDRO[b] - HYDRO[a]
    dc = CHARGE.get(b, 0) - CHARGE.get(a, 0)
    dv = VOL[b] - VOL[a]
    print(f'    delta Kyte-Doolittle {dh:+.1f} | delta charge {dc:+.1f} | delta volume {dv:+.1f} A^3')
    for lab, rx in MOTIFS:
        hv = {(m.start() + lo + 1, m.group()) for m in re.finditer(f'(?=({rx}))', wv)
              for m in [m]} if False else set()
        hv = {(m.start() + lo + 1, m.group(1)) for m in re.finditer(f'(?=({rx}))', wv)}
        hn = {(m.start() + lo + 1, m.group(1)) for m in re.finditer(f'(?=({rx}))', wn)}
        gained = sorted(hn - {(q, s) for q, s in hv} , key=lambda z: z[0])
        lost = sorted(hv - {(q, s) for q, s in hn}, key=lambda z: z[0])
        # only report changes that actually involve the substituted position
        gained = [g for g in gained if g[0] <= p <= g[0] + len(g[1]) - 1]
        lost = [g for g in lost if g[0] <= p <= g[0] + len(g[1]) - 1]
        if gained or lost:
            msg = []
            if gained:
                msg.append('CREATED in Dnai: ' + ', '.join(f'{s}@{q}' for q, s in gained))
            if lost:
                msg.append('DESTROYED in Dnai: ' + ', '.join(f'{s}@{q}' for q, s in lost))
            print(f'      * {lab}: ' + ' | '.join(msg))
