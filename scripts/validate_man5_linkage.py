#!/usr/bin/env python3
"""Analyze ASN–first-GlcNAc linkage geometry in a Man5 conformer NPZ.

Assumes attachment_frame rows correspond to ASN CB, CG, ND2 and glycan
coordinates are expressed in the same coordinate system. Check builder source.
"""
import argparse
import csv
from pathlib import Path
import numpy as np


def dihedral_batch(p0, p1, p2, p3):
    b0 = -(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    b1 = b1 / np.linalg.norm(b1, axis=-1, keepdims=True)
    v = b0 - np.sum(b0 * b1, axis=-1, keepdims=True) * b1
    w = b2 - np.sum(b2 * b1, axis=-1, keepdims=True) * b1
    return np.degrees(np.arctan2(np.sum(np.cross(b1, v) * w, axis=-1), np.sum(v * w, axis=-1)))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--library', required=True)
    ap.add_argument('--out-prefix', required=True)
    ap.add_argument('--chunk-size', type=int, default=5000)
    args = ap.parse_args()
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    out = Path(args.out_prefix)
    out.parent.mkdir(parents=True, exist_ok=True)
    with np.load(args.library, allow_pickle=False) as z:
        coords = z['coords']
        names = z['atom_names'].astype(str)
        residues = z['residue_numbers']
        frame = z['attachment_frame'].astype(float)
        assert frame.shape == (3, 3), 'Expected 3x3 CB/CG/ND2 attachment frame'
        def idx(resid, atom):
            matches = np.flatnonzero((residues == resid) & (names == atom))
            if len(matches) != 1:
                raise ValueError(f'Expected exactly one residue {resid} atom {atom}; found {len(matches)}')
            return int(matches[0])
        c1_i, o5_i = idx(1, 'C1'), idx(1, 'O5')
        cb, cg, nd2 = frame
        n = len(coords)
        phi = np.empty(n)
        psi = np.empty(n)
        bond = np.empty(n)
        for start in range(0, n, args.chunk_size):
            stop = min(start + args.chunk_size, n)
            c1 = coords[start:stop, c1_i, :]
            o5 = coords[start:stop, o5_i, :]
            m = stop-start
            Cb = np.broadcast_to(cb, (m,3))
            Cg = np.broadcast_to(cg, (m,3))
            Nd2 = np.broadcast_to(nd2, (m,3))
            phi[start:stop] = dihedral_batch(o5, c1, Nd2, Cg)
            psi[start:stop] = dihedral_batch(c1, Nd2, Cg, Cb)
            bond[start:stop] = np.linalg.norm(c1-Nd2, axis=1)
        extras = {k:z[k] for k in ('source_position','source_pdb','glycan_role','conformer_origin') if k in z.files}
    with open(f'{out}_linkage.csv', 'w', newline='') as f:
        fields = ['index','phi_deg','psi_deg','nd2_c1_angstrom'] + list(extras)
        writer = csv.writer(f)
        writer.writerow(fields)
        for i in range(n):
            writer.writerow([i,phi[i],psi[i],bond[i]] + [extras[k][i] for k in extras])
    for key, arr, xlabel, bins, limits in [
        ('phi',phi,'Phi (degrees)',72,(-180,180)),
        ('psi',psi,'Psi (degrees)',72,(-180,180)),
        ('bond',bond,'ND2–C1 distance (Å)',70,None)]:
        fig, ax = plt.subplots(figsize=(7,4.5))
        ax.hist(arr, bins=bins, range=limits)
        ax.set_xlabel(xlabel)
        ax.set_ylabel('Conformer count')
        fig.tight_layout()
        fig.savefig(f'{out}_{key}.png',dpi=180)
        plt.close(fig)
    fig, ax = plt.subplots(figsize=(7,6))
    h = ax.hist2d(phi,psi,bins=72,range=[[-180,180],[-180,180]],cmap='viridis')
    ax.set_xlabel('Phi (degrees): O5–C1–ND2–CG')
    ax.set_ylabel('Psi (degrees): C1–ND2–CG–CB')
    fig.colorbar(h[3],ax=ax,label='Conformer count')
    fig.tight_layout()
    fig.savefig(f'{out}_phi_psi.png',dpi=180)
    plt.close(fig)
    for key, arr in [('Phi (degrees)',phi),('Psi (degrees)',psi),('ND2–C1 (Å)',bond)]:
        print(f'{key}: n={len(arr)}, min={np.min(arr):.3f}, median={np.median(arr):.3f}, max={np.max(arr):.3f}')
    print('NOTE: Validate attachment_frame row order and shared coordinate frame against builder source.')
    print(f'Wrote {out}_linkage.csv and four plots.')

if __name__ == '__main__':
    main()
