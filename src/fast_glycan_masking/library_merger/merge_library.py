#!/usr/bin/env python3
"""Merge two Man5 NPZ libraries, normalizing heavy atoms and ASN attachment frames.

Requires numpy and scipy. Supports Rosetta empirical builder and torsion-library
formats. Preserves per-conformer provenance, with missing values filled; stores
original static metadata in a separate JSON report.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

ATOM = {'coords','atom_names','residue_numbers','residue_names','elements','chains','attachment_frame','linkage_edges','library_kind'}

def load(path):
    with np.load(path,allow_pickle=False) as z: d={k:z[k] for k in z.files}
    for k in ('coords','atom_names','residue_numbers','residue_names','elements','attachment_frame'):
        if k not in d: raise ValueError(f'{path}: missing {k}')
    if d['coords'].ndim!=3 or d['coords'].shape[2]!=3 or d['attachment_frame'].shape!=(3,3):
        raise ValueError(f'{path}: invalid coords/frame shapes')
    n,a,_=d['coords'].shape
    if n==0 or any(len(d[k])!=a for k in ('atom_names','residue_numbers','residue_names','elements')):
        raise ValueError(f'{path}: empty or inconsistent atom arrays')
    if not np.isfinite(d['coords']).all() or not np.isfinite(d['attachment_frame']).all():
        raise ValueError(f'{path}: nonfinite coordinates')
    return d

def normalized_atoms(d):
    r=d['residue_numbers']; uniq=list(dict.fromkeys(map(int,r)))
    order={v:i+1 for i,v in enumerate(uniq)}
    mask=np.array([str(e).upper() not in ('H','D') for e in d['elements']])
    ix=np.where(mask)[0]
    keys=[(order[int(r[i])],str(d['residue_names'][i]),str(d['atom_names'][i]),str(d['elements'][i]).upper()) for i in ix]
    if len(set(keys))!=len(keys): raise ValueError('Ambiguous repeated heavy-atom keys')
    return ix,keys

def edges(d):
    if 'linkage_edges' in d:
        e=d['linkage_edges'].astype(str)
        if e.ndim!=2 or e.shape[1]!=6: raise ValueError('Invalid linkage_edges')
        return sorted(tuple(x) for x in e)
    required=['linkage_parent_resid','linkage_child_resid','linkage_parent_atom','linkage_child_atom']
    if not all(k in d for k in required): raise ValueError('Missing glycan linkage metadata')
    residues=list(dict.fromkeys(map(int,d['residue_numbers'])))
    renum={v:i+1 for i,v in enumerate(residues)}
    names={int(r):str(d['residue_names'][np.where(d['residue_numbers']==r)[0][0]]) for r in residues}
    out=[]
    for p,c,pa,ca in zip(*(d[k] for k in required)):
        p,c=int(p),int(c)
        out.append((str(renum[p]),str(renum[c]),str(pa),str(ca),names[p],names[c]))
    return sorted(out)

def kabsch(source,target):
    sc,tc=source.mean(0),target.mean(0)
    u,_,vt=np.linalg.svd((source-sc).T@(target-tc))
    r=vt.T@np.diag([1,1,np.sign(np.linalg.det(vt.T@u.T))])@u.T
    return r,tc-r@sc

def deduplicate(x,cutoff,batch=128):
    """First-kept RMSD filter, KD-tree for prior batches; exact within batch."""
    if cutoff is None:return np.arange(len(x))
    if cutoff<=0:raise ValueError('dedup cutoff must be > 0')
    flat=np.ascontiguousarray(x.reshape(len(x),-1))
    radius=cutoff*np.sqrt(x.shape[1])
    keep=[]
    for start in range(0,len(flat),batch):
        tree=cKDTree(flat[keep]) if keep else None
        local=[]
        for i in range(start,min(start+batch,len(flat))):
            if tree is not None and tree.query_ball_point(flat[i],radius):continue
            if local and np.any(np.linalg.norm(flat[local]-flat[i],axis=1)<=radius):continue
            local.append(i)
        keep.extend(local)
        if start%4096==0:print(f'[dedup] processed {start}/{len(flat)}; retained {len(keep)}',flush=True)
    return np.array(keep,dtype=int)

def per_conformer(d):
    n=len(d['coords'])
    return {k:v for k,v in d.items() if k not in ATOM and v.ndim>=1 and len(v)==n}

def to_json(v):
    a=np.asarray(v)
    return a.item() if a.ndim==0 else a.tolist()

def merge(a_path,b_path,out,cutoff):
    a,b=load(a_path),load(b_path)
    ia,ka=normalized_atoms(a);ib,kb=normalized_atoms(b)
    if set(ka)!=set(kb):
        raise ValueError(f'Heavy atom sets differ: only A={list(set(ka)-set(kb))[:5]}, only B={list(set(kb)-set(ka))[:5]}')
    if edges(a)!=edges(b):raise ValueError('Linkage topologies differ')
    lookup={k:i for i,k in enumerate(kb)}
    perm=np.array([lookup[k] for k in ka])
    r,t=kabsch(np.asarray(b['attachment_frame'],float),np.asarray(a['attachment_frame'],float))
    alignment=float(np.sqrt(np.mean(np.sum((b['attachment_frame']@r.T+t-a['attachment_frame'])**2,axis=1))))
    if alignment>0.3:raise ValueError(f'ASN attachment frames mismatch (RMSD {alignment:.3f} Å > 0.3 Å)')
    x=np.concatenate([a['coords'][:,ia,:],b['coords'][:,ib[perm],:]@r.T+t])
    print(f'[merge] {len(a["coords"])} + {len(b["coords"])} conformers; {len(ia)} heavy atoms; frame RMSD {alignment:.4f} Å',flush=True)
    keep=deduplicate(x,cutoff)
    na,nb=len(a['coords']),len(b['coords'])
    pa,pb=per_conformer(a),per_conformer(b)
    outdata={
        'coords':x[keep],
        'atom_names':np.asarray([k[2] for k in ka]),
        'residue_numbers':np.asarray([k[0] for k in ka],dtype=int),
        'residue_names':np.asarray([k[1] for k in ka]),
        'elements':np.asarray([k[3] for k in ka]),
        'chains':np.full(len(ka),'G'),
        'attachment_frame':a['attachment_frame'],
        'linkage_edges':np.asarray(edges(a)),
        'library_kind':np.array('merged_man5_heavy_atom'),
        'source_library':np.where(keep<na,str(Path(a_path).resolve()),str(Path(b_path).resolve())),
        'original_library_index':np.where(keep<na,keep,keep-na),
    }
    for k in sorted(set(pa)|set(pb)):
        va,vb=pa.get(k),pb.get(k)
        sample=va if va is not None else vb
        shape=sample.shape[1:]
        if va is not None and va.shape[1:]!=shape or vb is not None and vb.shape[1:]!=shape:
            raise ValueError(f'Incompatible provenance dimensions: {k}')
        def fill(n):
            if sample.dtype.kind in 'f':return np.full((n,)+shape,np.nan)
            if sample.dtype.kind in 'iu':return np.full((n,)+shape,-1,dtype=np.int64)
            if sample.dtype.kind=='b':return np.zeros((n,)+shape,dtype=bool)
            return np.full((n,)+shape,'',dtype='<U1')
        if va is None:va=fill(na)
        if vb is None:vb=fill(nb)
        if va.dtype.kind in 'US' or vb.dtype.kind in 'US':va,vb=va.astype(str),vb.astype(str)
        outdata[k]=np.concatenate([va,vb])[keep]
    if 'parent_empirical_index' in outdata:
        old=outdata['parent_empirical_index'].copy()
        remap={int(v):i for i,v in enumerate(keep)}
        outdata['original_parent_empirical_index']=old.copy()
        outdata['parent_empirical_index']=np.array([remap.get(int(p)+(0 if int(i)<na else na),-1) if int(p)>=0 else -1 for i,p in zip(keep,old)],dtype=int)
    out=Path(out);out.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(out,**outdata)
    static={}
    for label,d in [('a',a),('b',b)]:
        static[label]={k:to_json(v) for k,v in d.items() if k not in ATOM and k not in per_conformer(d)}
    summary={'library_a':str(Path(a_path).resolve()),'library_b':str(Path(b_path).resolve()),'output':str(out.resolve()),'input_a':na,'input_b':nb,'total':na+nb,'retained':len(keep),'duplicates_removed':na+nb-len(keep),'retained_a':int(np.sum(keep<na)),'retained_b':int(np.sum(keep>=na)),'heavy_atoms':len(ka),'dedup_rmsd_angstrom':cutoff,'asn_frame_alignment_rmsd_angstrom':alignment,'input_atom_counts':[a['coords'].shape[1],b['coords'].shape[1]],'original_static_metadata':static,'note':'Output is heavy-atom only. Missing provenance filled with NaN/-1/empty strings. Parent indices -1 when parent is absent.'}
    report=out.with_name(out.stem+'_merge_summary.json');report.write_text(json.dumps(summary,indent=2)+'\n')
    return summary,report

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--library-a',required=True);p.add_argument('--library-b',required=True)
    p.add_argument('--out',required=True)
    g=p.add_mutually_exclusive_group();g.add_argument('--dedup-rmsd',type=float,default=0.5);g.add_argument('--no-dedup',action='store_true')
    a=p.parse_args();summary,report=merge(a.library_a,a.library_b,a.out,None if a.no_dedup else a.dedup_rmsd)
    print(json.dumps({k:v for k,v in summary.items() if k!='original_static_metadata'},indent=2));print('[done]',report)
if __name__=='__main__':main()
