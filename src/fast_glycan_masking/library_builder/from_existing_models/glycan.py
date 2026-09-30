"""Parse Rosetta glycoprotein PDBs and recover ASN-linked glycan trees.

The routines in this module identify N-glycans from PDB LINK records, retain the
attachment ASN frame (CB/CG/ND2), and canonicalize atom/residue ordering so
conformers from different source models can coexist in one NPZ library.
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple, Set
import numpy as np

AtomKey = Tuple[str,int,str]
ResidueKey = Tuple[str,int]

@dataclass(frozen=True)
class Atom:
    serial: int
    name: str
    resname: str
    chain: str
    resid: int
    element: str
    xyz: Tuple[float,float,float]

@dataclass(frozen=True)
class Link:
    atom1: AtomKey
    resname1: str
    atom2: AtomKey
    resname2: str

@dataclass
class GlycanTree:
    attachment_chain: str
    attachment_resid: int
    root: ResidueKey
    residues: Set[ResidueKey]
    links: List[Link]
    atoms: List[Atom]
    asn_frame: np.ndarray

def parse_atom(line):
    name=line[12:16].strip(); resname=line[17:20].strip()
    chain=line[21:22].strip() or "A"; resid=int(line[22:26])
    element=line[76:78].strip().upper()
    if not element:
        element = next((c for c in name if c.isalpha()), "C").upper()
    return Atom(int(line[6:11]),name,resname,chain,resid,element,
                (float(line[30:38]),float(line[38:46]),float(line[46:54])))

def parse_link(line):
    return Link(
        (line[21:22].strip() or "A", int(line[22:26]), line[12:16].strip()),
        line[17:20].strip(),
        (line[51:52].strip() or "A", int(line[52:56]), line[42:46].strip()),
        line[47:50].strip(),
    )

def read_pdb(path: Path):
    atoms=[]; links=[]
    with Path(path).open() as fh:
        for line in fh:
            if line.startswith(("ATOM  ","HETATM")):
                atoms.append(parse_atom(line))
            elif line.startswith("LINK"):
                links.append(parse_link(line))
    return atoms, links

def _is_asn_nd2(key, resname):
    return key[2] == "ND2" and resname.upper() == "ASN"

def find_n_glycans(path: Path) -> List[GlycanTree]:
    """Find ASN-linked glycan trees and their local CB/CG/ND2 attachment frames."""
    atoms, links = read_pdb(path)
    atom_map={(a.chain,a.resid,a.name):a for a in atoms}
    residue_links: Dict[ResidueKey,List[Tuple[ResidueKey,Link]]] = {}
    attachments=[]
    for lk in links:
        r1=(lk.atom1[0],lk.atom1[1]); r2=(lk.atom2[0],lk.atom2[1])
        residue_links.setdefault(r1,[]).append((r2,lk))
        residue_links.setdefault(r2,[]).append((r1,lk))
        if _is_asn_nd2(lk.atom1,lk.resname1) and lk.atom2[2]=="C1":
            attachments.append((lk.atom1[0],lk.atom1[1],r2,lk))
        elif _is_asn_nd2(lk.atom2,lk.resname2) and lk.atom1[2]=="C1":
            attachments.append((lk.atom2[0],lk.atom2[1],r1,lk))

    trees=[]
    for chain,resid,root,attach in attachments:
        # Traverse carbohydrate LINK graph without walking back into protein.
        seen=set(); stack=[root]
        while stack:
            r=stack.pop()
            if r in seen: continue
            seen.add(r)
            for nbr,lk in residue_links.get(r,[]):
                if nbr == (chain,resid):
                    continue
                # Protein residues are represented by ATOM records; glycan residues
                # in these Rosetta models are HETATM and connected through C1/Ox links.
                if any(a.chain==nbr[0] and a.resid==nbr[1] and a.resname=="ASN"
                       for a in atoms):
                    continue
                if nbr not in seen:
                    stack.append(nbr)
        gatoms=[a for a in atoms if (a.chain,a.resid) in seen]
        needed=[(chain,resid,"CB"),(chain,resid,"CG"),(chain,resid,"ND2")]
        try:
            frame=np.array([atom_map[k].xyz for k in needed],float)
        except KeyError as e:
            raise ValueError(f"{path}: ASN {chain}:{resid} missing {e.args[0][2]}")
        glinks=[lk for lk in links
                if ((lk.atom1[0],lk.atom1[1]) in seen and
                    (lk.atom2[0],lk.atom2[1]) in seen)]
        trees.append(GlycanTree(chain,resid,root,seen,glinks,gatoms,frame))
    return trees

def canonicalize_tree(tree: GlycanTree):
    """Return a stable glycan atom order/signature suitable for cross-model libraries."""
    """Return deterministic atom order/signature independent of PDB residue numbers."""
    adj: Dict[ResidueKey,List[Tuple[ResidueKey,Link]]] = {}
    for lk in tree.links:
        r1=(lk.atom1[0],lk.atom1[1]); r2=(lk.atom2[0],lk.atom2[1])
        adj.setdefault(r1,[]).append((r2,lk))
        adj.setdefault(r2,[]).append((r1,lk))

    def acceptor_position(lk, parent):
        atom = lk.atom1[2] if (lk.atom1[0],lk.atom1[1])==parent else lk.atom2[2]
        if atom.startswith("O") and atom[1:].isdigit():
            return int(atom[1:])
        return 99

    parent={tree.root:None}; order=[]; queue=[tree.root]
    while queue:
        r=queue.pop(0); order.append(r)
        kids=[]
        for nbr,lk in adj.get(r,[]):
            if nbr == parent.get(r): continue
            if nbr in parent: continue
            parent[nbr]=r; kids.append((acceptor_position(lk,r),nbr))
        kids.sort(key=lambda x:(x[0],x[1][0],x[1][1]))
        queue.extend([x[1] for x in kids])

    rid={r:i+1 for i,r in enumerate(order)}
    byres={}
    for a in tree.atoms: byres.setdefault((a.chain,a.resid),[]).append(a)
    canon_atoms=[]; source_atoms=[]; signature=[]
    serial=1
    for r in order:
        for a in sorted(byres[r], key=lambda x:x.name):
            source_atoms.append(a)
            canon_atoms.append(Atom(serial,a.name,a.resname,"G",rid[r],a.element,a.xyz))
            signature.append((rid[r],a.resname,a.name,a.element))
            serial+=1

    edges=[]
    for lk in tree.links:
        r1=(lk.atom1[0],lk.atom1[1]); r2=(lk.atom2[0],lk.atom2[1])
        if r1 in rid and r2 in rid:
            edges.append((rid[r1],rid[r2],lk.atom1[2],lk.atom2[2],
                          lk.resname1,lk.resname2))
    edges=tuple(sorted(edges))
    return source_atoms, canon_atoms, tuple(signature), edges
