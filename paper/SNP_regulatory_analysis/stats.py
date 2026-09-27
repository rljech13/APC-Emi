"""Quantitative tests:
  (1) Are the ten Dval/Dnai substitutions enriched in mapped functional elements?
  (2) Are they enriched in ordered (vs disordered) regions?
  (3) Column-wise conservation percentile of each of the ten positions.
  (4) Sequence-quality cross-check: study's own Podarcis muralis vs RefSeq.
"""
import os, sys, math, random
from collections import Counter
from Bio import AlignIO

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from seqs import read_fa

ALN = AlignIO.read(os.path.join(HERE, 'emi2_orthologues_linsi.aln.fasta'), 'fasta')
S = {r.id.split('|')[0]: str(r.seq) for r in ALN}
REF = 'Darevskia_valentini'
L = 675
SUBS = [76, 92, 102, 171, 219, 391, 395, 425, 549, 590]


def index(tag):
    r2c, c2r, n = {}, {}, 0
    for c, ch in enumerate(S[tag]):
        if ch != '-':
            n += 1
            r2c[n] = c
            c2r[c] = n
    return r2c, c2r


I = {t: index(t) for t in S}

# ---------------------------------------------------------------- (1) elements
ELEMENTS = [(40, 46, 'beta-TrCP DSGx3S degron'),
            (176, 207, 'Plk1-PBD docking module'),
            (358, 368, 'PP2A-B56 LxxIxE + p90rsk cluster'),
            # the annotated folded domain, not the 41-aa Skp1-binding motif 479-519
            # nested inside it: the footprint is a functional-element footprint
            (452, 537, 'F-box domain'),
            (548, 551, 'D-box'),
            (574, 575, 'Cdk1 T574-P'),
            (600, 648, 'ZBR'),
            # post-ZBR and the RL tail are merged: under the current boundaries the
            # RL tail starts at 658 and would otherwise overlap post-ZBR 652-661,
            # which would double-count those residues in the footprint
            (652, 675, 'post-ZBR + RL tail (Apc2-binding)')]
covered = set()
for s, e, _ in ELEMENTS:
    covered |= set(range(s, e + 1))
obs = sum(1 for p in SUBS if p in covered)
exp = len(SUBS) * len(covered) / L
# exact hypergeometric-style permutation
random.seed(0)
N = 200000
ge = sum(1 for _ in range(N)
         if sum(1 for p in random.sample(range(1, L + 1), len(SUBS)) if p in covered) >= obs)
print('=' * 88)
print('(1) ENRICHMENT OF SUBSTITUTIONS IN MAPPED FUNCTIONAL ELEMENTS')
print('=' * 88)
for s, e, lab in ELEMENTS:
    inside = [p for p in SUBS if s <= p <= e]
    print(f'  {lab:36s} {s:>3}-{e:<3} ({e-s+1:>3} aa)   substitutions inside: {inside}')
print(f'\n  functional-element footprint: {len(covered)}/{L} residues = {len(covered)/L:.1%}')
print(f'  observed substitutions inside: {obs}   expected by chance: {exp:.2f}')
print(f'  permutation P(>= {obs} by chance) = {ge/N:.3f}   -> NO enrichment' if ge / N > 0.05
      else f'  permutation P = {ge/N:.4f}')

# ---------------------------------------------------------------- (2) disorder
from disorder import foldindex  # noqa: E402
dval = read_fa(os.path.join(os.path.dirname(__file__), "..", "Fasta", "Emi2", "Emi2_Dval_AAseq.fasta"))['Emi2_Dval_AAseq']

# ---------------------------------------------------------------- (3) conservation
print('\n' + '=' * 88)
print('(3) COLUMN CONSERVATION PERCENTILE OF THE TEN POSITIONS')
print('=' * 88)
EXCL = {t for t in S if t.startswith('Darevskia') or t.endswith('_study')}
TAXA = [t for t in S if t not in EXCL]


def col_identity(col):
    """fraction of non-gap taxa carrying the modal residue"""
    c = Counter(S[t][col] for t in TAXA if S[t][col] != '-')
    if not c:
        return 0.0, 0
    return c.most_common(1)[0][1] / sum(c.values()), sum(c.values())


# score every Dval-aligned column
scores = {}
for p in range(1, L + 1):
    col = I[REF][0][p]
    frac, n = col_identity(col)
    if n >= len(TAXA) * 0.7:
        scores[p] = frac
allvals = sorted(scores.values())


def pct(v):
    import bisect
    return bisect.bisect_left(allvals, v) / len(allvals) * 100


print(f'  (scored {len(scores)} of {L} Dval positions with >=70% taxon occupancy; n taxa = {len(TAXA)})')
print(f'  median column identity across Emi2 = {allvals[len(allvals)//2]:.2f}')
print(f'\n  {"pos":>5} {"modal-residue identity":>22} {"conservation percentile":>24}')
for p in SUBS:
    v = scores.get(p)
    if v is None:
        print(f'  {p:>5} {"low occupancy":>22}')
    else:
        print(f'  {p:>5} {v:>21.2f} {pct(v):>23.0f}')

# ---------------------------------------------------------------- (4) QC
print('\n' + '=' * 88)
print('(4) SEQUENCE-QUALITY CROSS-CHECK: study Podarcis muralis vs NCBI RefSeq')
print('=' * 88)
a, b = S['Podarcis_muralis_study'], S['Podarcis_muralis']
d = [(i, x, y) for i, (x, y) in enumerate(zip(a, b)) if x != y]
print(f'  aligned differences: {len(d)}')
for i, x, y in d[:20]:
    dp = I[REF][1].get(i)
    print(f'    aln col {i} (Dval pos {dp}): study={x} RefSeq={y}')

a2, b2 = S['Eublepharis_macularius_study'], S['Eublepharis_macularius']
d2 = [(i, x, y) for i, (x, y) in enumerate(zip(a2, b2)) if x != y]
print(f'\n  study Eublepharis macularius vs RefSeq: {len(d2)} differences')
