"""Optional local torsion expansion for empirical glycan conformers.

Sampling rotates the glycan around the ASN ND2--C1 attachment bond and around
identified glycosidic linkage bonds. The operation changes torsional orientation
while preserving bond lengths and bond angles.
"""
from __future__ import annotations
from typing import Dict, List, Set, Tuple
import numpy as np
from .geometry import rotate_points

def _residue_atom_indices(atom_names, residue_numbers):
    out={}
    for i,r in enumerate(residue_numbers):
        out.setdefault(int(r),[]).append(i)
    return out

def _descendants(edges, root_child):
    children={}
    for p,c,pa,ca,*_ in edges:
        children.setdefault(int(p),[]).append(int(c))
    seen=set(); stack=[int(root_child)]
    while stack:
        r=stack.pop()
        if r in seen: continue
        seen.add(r); stack.extend(children.get(r,[]))
    return seen

def expand_conformer(coords, atom_names, residue_numbers, edges,
                     canonical_asn_frame, n_samples, rng,
                     attachment_sigma_deg=20.0, glycosidic_sigma_deg=20.0):
    """Generate local torsional variants while preserving bond lengths/angles.

    The empirical conformer itself is not returned here; callers keep it separately.
    Rotations are around the ASN ND2--root C1 bond and carbohydrate linkage bonds.
    """
    coords=np.asarray(coords,float)
    names=np.asarray(atom_names); resnums=np.asarray(residue_numbers,int)
    res_atoms=_residue_atom_indices(names,resnums)

    def atom_idx(resid,name):
        hits=np.where((resnums==int(resid)) & (names==name))[0]
        if not len(hits): raise ValueError(f"missing atom {resid}:{name}")
        return int(hits[0])

    root=1
    c1=atom_idx(root,"C1")
    nd2=np.asarray(canonical_asn_frame,float)[2]
    out=[]
    for _ in range(n_samples):
        x=coords.copy()
        # Root attachment torsion: rotate the whole glycan except C1 around ND2-C1.
        delta=float(rng.normal(0,attachment_sigma_deg))
        movers=[i for i in range(len(x)) if i != c1]
        x=rotate_points(x,movers,nd2,x[c1],delta)

        # Each glycosidic LINK gives a physically interpretable bond axis.
        for p,c,parent_atom,child_atom,*_ in edges:
            if child_atom != "C1":
                continue
            try:
                ia=atom_idx(p,parent_atom); ib=atom_idx(c,child_atom)
            except ValueError:
                continue
            downstream=_descendants(edges,c)
            movers=[]
            for rr in downstream:
                movers.extend(res_atoms.get(rr,[]))
            movers=[i for i in movers if i not in (ia,ib)]
            if movers:
                delta=float(rng.normal(0,glycosidic_sigma_deg))
                x=rotate_points(x,movers,x[ia],x[ib],delta)
        out.append(x)
    return out
