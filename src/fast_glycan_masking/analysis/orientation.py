#!/usr/bin/env python3
"""Analyze glycan orientation diversity in a library or placer output."""

from __future__ import annotations
import argparse, csv, json, re
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


def parse_site(s):
    m = re.fullmatch(r"\s*([^:\s]+)\s*:\s*(-?\d+)\s*", s)
    if not m:
        raise argparse.ArgumentTypeError("Use CHAIN:RESID, e.g. A:581")
    return m.group(1), int(m.group(2))


def nd2_from_pdb(path, chain, resid):
    with open(path) as fh:
        for line in fh:
            if line[:6].strip() not in {"ATOM", "HETATM"}: continue
            if (line[21:22].strip() or "A") != chain: continue
            try: r = int(line[22:26])
            except ValueError: continue
            if r == resid and line[12:16].strip() == "ND2":
                return np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    raise ValueError(f"ND2 not found for {chain}:{resid}")


def normalize(v):
    """Normalize an array of 3D vectors to unit length."""
    n = np.linalg.norm(v, axis=1, keepdims=True)
    if np.any(n < 1e-12): raise ValueError("Zero-length orientation vector")
    return v/n


def angles(u):
    """Convert unit vectors to azimuth and polar angles in degrees."""
    az = np.mod(np.degrees(np.arctan2(u[:,1], u[:,0])), 360.0)
    polar = np.degrees(np.arccos(np.clip(u[:,2], -1, 1)))
    return az, polar


def occupancy(az, polar, naz=12, nmu=6):
    mu = np.cos(np.radians(polar))
    h, _, _ = np.histogram2d(az, mu,
        bins=[np.linspace(0,360,naz+1), np.linspace(-1,1,nmu+1)])
    f = h.ravel()
    occupied = int(np.count_nonzero(f))
    p = f[f>0]/f.sum()
    ent = float(-np.sum(p*np.log(p))/np.log(len(f))) if len(p) else 0.0
    return h, {"occupied_bins": occupied, "total_bins": len(f),
               "occupancy_fraction": occupied/len(f),
               "normalized_entropy": ent}


def angular_stats(u):
    n=len(u)
    if n < 2: return {"n": n}
    # exact nearest-neighbor angles, blockwise so 60k libraries remain practical
    nn=np.full(n, np.inf)
    block=256
    for a in range(0,n,block):
        b=min(a+block,n)
        d=np.clip(u[a:b] @ u.T,-1,1)
        ang=np.degrees(np.arccos(d))
        ang[np.arange(b-a), np.arange(a,b)] = np.inf
        nn[a:b]=ang.min(axis=1)
    # reproducible random sample for global pairwise distribution
    rng=np.random.default_rng(2026)
    k=min(200000,max(10000,n*10))
    i=rng.integers(0,n,k); j=rng.integers(0,n,k)
    good=i!=j
    pair=np.degrees(np.arccos(np.clip(np.sum(u[i[good]]*u[j[good]],axis=1),-1,1)))
    return {
        "n": n,
        "nearest_neighbor_min_deg": float(nn.min()),
        "nearest_neighbor_mean_deg": float(nn.mean()),
        "nearest_neighbor_median_deg": float(np.median(nn)),
        "nearest_neighbor_max_deg": float(nn.max()),
        "sampled_pairwise_mean_deg": float(pair.mean()),
        "sampled_pairwise_median_deg": float(np.median(pair)),
    }


def save_outputs(vectors, labels, prefix, tag):
    u=normalize(vectors); az,polar=angles(u)
    p=Path(f"{prefix}_{tag}_orientations.csv")
    with open(p,"w",newline="") as f:
        w=csv.writer(f); w.writerow(["model","ux","uy","uz","azimuth_deg","polar_deg"])
        for lab,x,a,b in zip(labels,u,az,polar): w.writerow([lab,*x,float(a),float(b)])

    fig,ax=plt.subplots(figsize=(8,5)); ax.hist(az,bins=np.linspace(0,360,37))
    ax.set(xlabel="Azimuth (degrees)",ylabel="Count",title=f"{tag}: azimuth")
    fig.tight_layout(); fig.savefig(f"{prefix}_{tag}_azimuth.png",dpi=200); plt.close(fig)

    fig,ax=plt.subplots(figsize=(8,5)); ax.hist(polar,bins=np.linspace(0,180,19))
    ax.set(xlabel="Polar angle (degrees)",ylabel="Count",title=f"{tag}: polar angle")
    fig.tight_layout(); fig.savefig(f"{prefix}_{tag}_polar.png",dpi=200); plt.close(fig)

    h,occ=occupancy(az,polar)
    fig,ax=plt.subplots(figsize=(9,5))
    im=ax.imshow(h.T,origin="lower",aspect="auto",extent=[0,360,-1,1])
    ax.set(xlabel="Azimuth (degrees)",ylabel="cos(polar angle)",
           title=f"{tag}: equal-area spherical occupancy")
    fig.colorbar(im,ax=ax,label="Count")
    fig.tight_layout(); fig.savefig(f"{prefix}_{tag}_sphere.png",dpi=200); plt.close(fig)
    return {**angular_stats(u),**occ}


def analyze_library(path,prefix):
    """Analyze root and whole-glycan orientation distributions in an NPZ library."""
    with np.load(path,allow_pickle=False) as z:
        c=np.asarray(z["coords"],float)
        frame=np.asarray(z["attachment_frame"],float)
        names=np.asarray(z["atom_names"]).astype(str)
        resids=np.asarray(z["residue_numbers"],int)
    if c.ndim!=3: raise ValueError(f"Library coords should be (N,A,3), got {c.shape}")
    nd2=frame[2]
    result={"input_type":"library","library":str(Path(path).resolve())}
    result["whole_glycan"]=save_outputs(c.mean(1)-nd2,np.arange(len(c)),prefix,"center")
    root=resids[0]
    hits=np.where((resids==root)&(names=="C1"))[0]
    if len(hits):
        result["root_c1"]=save_outputs(c[:,hits[0],:]-nd2,np.arange(len(c)),prefix,"root")
    return result


def analyze_placed(path,protein,site,prefix):
    """Analyze directional diversity of glycans written by the placement module."""
    chain,resid=site; nd2=nd2_from_pdb(protein,chain,resid)
    with np.load(path,allow_pickle=False) as z:
        c=np.asarray(z["coords"],float)
        labels=np.asarray(z["site_labels"]).astype(str)
        selected=np.asarray(z["selected_library_indices"],int) if "selected_library_indices" in z.files else None
    if c.ndim==3: c=c[:,None,:,:]
    wanted=f"{chain}{resid}"
    hit=np.where(labels==wanted)[0]
    if len(hit)!=1: raise ValueError(f"{wanted} not uniquely found; sites={labels.tolist()}")
    si=int(hit[0]); sc=c[:,si,:,:]
    model_labels=selected[:,si] if selected is not None else np.arange(len(sc))
    return {"input_type":"placed","placed":str(Path(path).resolve()),"site":f"{chain}:{resid}",
            "whole_glycan":save_outputs(sc.mean(1)-nd2,model_labels,prefix,"placed")}


def main():
    ap=argparse.ArgumentParser()
    g=ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--library"); g.add_argument("--placed")
    ap.add_argument("--protein"); ap.add_argument("--site",type=parse_site)
    ap.add_argument("--out-prefix",required=True)
    a=ap.parse_args()
    Path(a.out_prefix).parent.mkdir(parents=True,exist_ok=True)
    if a.library: summary=analyze_library(a.library,a.out_prefix)
    else:
        if not a.protein or not a.site: ap.error("--placed requires --protein and --site")
        summary=analyze_placed(a.placed,a.protein,a.site,a.out_prefix)
    out=f"{a.out_prefix}_orientation_summary.json"
    Path(out).write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2)); print(f"[done] {out}")


if __name__=="__main__": main()
