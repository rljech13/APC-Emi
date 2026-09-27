"""Structural context of Dval L425 / Dnai F425 (= human Leu449, Xenopus Met403).

Uses (a) the study's AF3 Emi2(400-675) monomer-in-complex models and
(b) the full-length AlphaFold-DB model of human FBXO43 (Q4G163),
computing per-residue pLDDT, Shrake-Rupley relative SASA, secondary structure
(backbone i->i+4 H-bond / CA geometry) and packing neighbours.
"""
import os, sys, math
from Bio.PDB import PDBParser, ShrakeRupley, DSSP
import warnings
warnings.filterwarnings('ignore')

HERE = os.path.dirname(os.path.abspath(__file__))
MODELS = os.path.join(os.path.dirname(__file__), "..", "Fasta", "Models") + os.sep

# Tien et al. 2013 theoretical maximum SASA
MAXASA = dict(A=129, R=274, N=195, D=193, C=167, E=223, Q=225, G=104, H=224,
              I=197, L=201, K=236, M=224, F=240, P=159, S=155, T=172, W=285,
              Y=263, V=174)
THREE2ONE = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLU': 'E',
             'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
             'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
             'TYR': 'Y', 'VAL': 'V'}

P = PDBParser(QUIET=True)


def analyse(path, label, targets, chain_id=None):
    st = P.get_structure('m', path)
    model = st[0]
    ch = model[chain_id] if chain_id else list(model)[0]
    sr = ShrakeRupley()
    sr.compute(model, level='R')
    print(f'\n===== {label} =====')
    print(f'  chains: {[c.id for c in model]}')
    resmap = {r.id[1]: r for r in ch if r.id[0] == ' '}
    print(f'  chain {ch.id} residue range: {min(resmap)}-{max(resmap)}')
    for t in targets:
        r = resmap.get(t)
        if r is None:
            print(f'  residue {t}: NOT PRESENT')
            continue
        aa = THREE2ONE.get(r.get_resname(), 'X')
        plddt = r['CA'].get_bfactor()
        rsa = r.sasa / MAXASA.get(aa, 200) * 100
        # packing neighbours: heavy atoms within 5 A of any side-chain atom
        sc = [a for a in r if a.get_id() not in ('N', 'CA', 'C', 'O')]
        near = set()
        for a in sc:
            for r2 in ch:
                if r2.id[0] != ' ' or r2.id[1] == t:
                    continue
                for a2 in r2:
                    if (a - a2) < 5.0:
                        near.add((r2.id[1], THREE2ONE.get(r2.get_resname(), 'X')))
                        break
        near = sorted(near)
        print(f'  {aa}{t}: pLDDT={plddt:.1f}  SASA={r.sasa:.1f} A^2  relSASA={rsa:.1f}%  '
              f'n_contacts(<5A)={len(near)}')
        print(f'        neighbours: ' + ', '.join(f'{b}{a}' for a, b in near))
    return ch, resmap


def helix_scan(ch, lo, hi):
    """crude alpha-helix detection: CA(i)-CA(i+4) distance ~6.2 A"""
    cas = {r.id[1]: r['CA'].coord for r in ch if r.id[0] == ' ' and 'CA' in r}
    out = []
    for i in range(lo, hi + 1):
        if i in cas and i + 4 in cas:
            d = float(((cas[i] - cas[i + 4]) ** 2).sum() ** 0.5)
            out.append((i, d, 'H' if 5.0 < d < 7.2 else '.'))
    return out


TARGETS = [425, 549, 590, 607]
ch_v, rv = analyse(MODELS + 'Emi2_Dval_interface400-675_AF3_model0.pdb',
                   'AF3 Dval Emi2(400-675) + Cdc20 + Apc10 (chain A)', TARGETS, 'A')
ch_n, rn = analyse(MODELS + 'Emi2_Dnai_interface400-675_AF3_model0.pdb',
                   'AF3 Dnai Emi2(400-675) + Cdc20 + Apc10 (chain A)', TARGETS, 'A')

print('\n===== helix scan around 425 (AF3 Dval, CA(i)-CA(i+4) distance) =====')
for i, d, c in helix_scan(ch_v, 408, 450):
    mark = '  <== 425' if i == 425 else ''
    print(f'   {i:4d} {d:5.2f} {c}{mark}')

# ---- human AFDB
print('\n===== AlphaFold-DB human FBXO43 Q4G163 (full length) =====')
st = P.get_structure('h', os.path.join(HERE, 'AF_Q4G163_v6.pdb'))
hm = st[0]
hch = list(hm)[0]
sr = ShrakeRupley()
sr.compute(hm, level='R')
hres = {r.id[1]: r for r in hch}
for t in [449]:
    r = hres[t]
    aa = THREE2ONE[r.get_resname()]
    rsa = r.sasa / MAXASA[aa] * 100
    print(f'  {aa}{t}: pLDDT={r["CA"].get_bfactor():.1f}  relSASA={rsa:.1f}%')
    sc = [a for a in r if a.get_id() not in ('N', 'CA', 'C', 'O')]
    near = set()
    for a in sc:
        for r2 in hch:
            if r2.id[1] == t:
                continue
            for a2 in r2:
                if (a - a2) < 5.0:
                    near.add((r2.id[1], THREE2ONE.get(r2.get_resname(), 'X')))
                    break
    print('        neighbours: ' + ', '.join(f'{b}{a}' for a, b in sorted(near)))

print('\n  human pLDDT profile 425-500 (10-residue means):')
for s in range(425, 500, 10):
    vals = [hres[i]['CA'].get_bfactor() for i in range(s, min(s + 10, 501)) if i in hres]
    print(f'    {s}-{s+len(vals)-1}: {sum(vals)/len(vals):.1f}')

print('\n===== helix scan around human 449 (AFDB) =====')
for i, d, c in helix_scan(hch, 430, 475):
    mark = '  <== 449' if i == 449 else ''
    print(f'   {i:4d} {d:5.2f} {c}{mark}')
