"""Disorder / order assessment for Darevskia Emi2.

Three independent lines:
  1. FoldIndex (Prilusky et al. 2005, Bioinformatics 21:3435) computed directly
     on the Darevskia valentini sequence.
  2. Per-residue pLDDT of the full-length AlphaFold-DB models of the human
     (Q4G163) and Xenopus laevis (Q8AXF4) orthologues, projected onto
     Darevskia numbering through the MAFFT L-INS-i alignment.
  3. Per-residue pLDDT of the study's own AF3 Emi2(400-675) models.
"""
import os, sys
from Bio import AlignIO

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from seqs import read_fa

# ---------------------------------------------------------------- FoldIndex
KD = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4,
          H=-3.2, I=4.5, L=3.8, K=-3.9, M=1.9, F=2.8, P=-1.6, S=-0.8,
          T=-0.7, W=-0.9, Y=-1.3, V=4.2)
CHG = dict(D=-1, E=-1, K=1, R=1)


def foldindex(seq, window=51):
    """FoldIndex = 2.785<H> - |<R>| - 1.151 ; <0 => predicted unfolded."""
    h = [(KD[a] + 4.5) / 9.0 for a in seq]          # scaled to [0,1]
    c = [CHG.get(a, 0) for a in seq]
    half, n, out = window // 2, len(seq), []
    for i in range(n):
        lo, hi = max(0, i - half), min(n, i + half + 1)
        w = hi - lo
        out.append(2.785 * (sum(h[lo:hi]) / w) - abs(sum(c[lo:hi]) / w) - 1.151)
    return out


# ---------------------------------------------------------------- pLDDT
def plddt_from_pdb(path):
    """{resnum: bfactor} using CA atoms."""
    d = {}
    for line in open(path):
        if line.startswith('ATOM') and line[12:16].strip() == 'CA':
            d[int(line[22:26])] = float(line[60:66])
    return d


# ---------------------------------------------------------------- mapping
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


def project_track(track, src_tag):
    """{src_resnum: value} -> {dval_resnum: value}"""
    out = {}
    for rp, v in track.items():
        col = I[src_tag][0].get(rp)
        if col is None:
            continue
        dp = I[REF][1].get(col)
        if dp:
            out[dp] = v
    return out


dval = read_fa(os.path.join(os.path.dirname(__file__), "..", "Fasta", "Emi2", "Emi2_Dval_AAseq.fasta"))['Emi2_Dval_AAseq']
fi = foldindex(dval)

hum = project_track(plddt_from_pdb(os.path.join(HERE, 'AF_Q4G163_v6.pdb')), 'Homo_sapiens')
xen = project_track(plddt_from_pdb(os.path.join(HERE, 'AF_Q8AXF4_v6.pdb')), 'Xenopus_laevis')

MODELS = os.path.join(os.path.dirname(__file__), "..", "Fasta", "Models") + os.sep
af3v = plddt_from_pdb(MODELS + 'Emi2_Dval_interface400-675_AF3_model0.pdb')
af3n = plddt_from_pdb(MODELS + 'Emi2_Dnai_interface400-675_AF3_model0.pdb')

SUBS = [(76, 'Q', 'H'), (92, 'V', 'I'), (102, 'A', 'T'), (171, 'K', 'R'),
        (219, 'E', 'D'), (391, 'I', 'V'), (395, 'D', 'N'), (425, 'L', 'F'),
        (549, 'F', 'I'), (590, 'V', 'M')]

print('=' * 96)
print('DISORDER / ORDER AT THE TEN POLYMORPHIC POSITIONS')
print('=' * 96)
print(f'{"pos":>5} {"sub":>6} {"FoldIdx":>8} {"call":>10} {"pLDDT_hum":>10} {"pLDDT_xen":>10} {"AF3_Dval":>9} {"AF3_Dnai":>9}')
for p, a, b in SUBS:
    f = fi[p - 1]
    call = 'disordered' if f < 0 else 'folded'
    h = hum.get(p)
    x = xen.get(p)
    v = af3v.get(p)
    n = af3n.get(p)
    fmt = lambda z: f'{z:10.1f}' if z is not None else f'{"n/a":>10}'
    fmt9 = lambda z: f'{z:9.1f}' if z is not None else f'{"n/a":>9}'
    print(f'{p:>5} {a+str(p)+b:>6} {f:8.3f} {call:>10} {fmt(h)} {fmt(x)} {fmt9(v)} {fmt9(n)}')

# ordered segments over the whole protein (from AF3 model, 400-675)
print('\n' + '=' * 96)
print('AF3 Dval Emi2(400-675) pLDDT profile: contiguous segments with pLDDT >= 70')
print('=' * 96)
res = sorted(af3v)
seg, segs = [], []
for r in res:
    if af3v[r] >= 70:
        seg.append(r)
    else:
        if len(seg) >= 4:
            segs.append((seg[0], seg[-1]))
        seg = []
if len(seg) >= 4:
    segs.append((seg[0], seg[-1]))
for s, e in segs:
    mean = sum(af3v[i] for i in range(s, e + 1)) / (e - s + 1)
    print(f'  {s}-{e}  (len {e-s+1}, mean pLDDT {mean:.1f})')
lo = [r for r in res if af3v[r] < 70]
print(f'\n  residues with pLDDT < 70: n={len(lo)} of {len(res)}')

print('\n' + '=' * 96)
print('FoldIndex profile, 25-residue blocks (negative = predicted disordered)')
print('=' * 96)
for s in range(0, len(dval), 25):
    blk = fi[s:s + 25]
    m = sum(blk) / len(blk)
    bar = '#' * max(0, int((m + 0.4) * 40))
    print(f'  {s+1:4d}-{min(s+25,len(dval)):4d}  {m:+.3f} {bar}')
