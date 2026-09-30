"""Rigid-body and rotation geometry used by the empirical library builder."""
from __future__ import annotations
import numpy as np

def kabsch(source: np.ndarray, target: np.ndarray):
    """Return R,t such that source @ R.T + t best overlays target."""
    source = np.asarray(source, float)
    target = np.asarray(target, float)
    sc = source.mean(axis=0)
    tc = target.mean(axis=0)
    h = (source - sc).T @ (target - tc)
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    r = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    t = tc - r @ sc
    return r, t

def transform(coords: np.ndarray, source_frame: np.ndarray, target_frame: np.ndarray):
    r, t = kabsch(source_frame, target_frame)
    return np.asarray(coords, float) @ r.T + t

def dihedral_deg(p0, p1, p2, p3):
    p0,p1,p2,p3 = map(lambda x: np.asarray(x,float), (p0,p1,p2,p3))
    b0 = -(p1-p0); b1 = p2-p1; b2 = p3-p2
    b1 /= np.linalg.norm(b1)
    v = b0 - np.dot(b0,b1)*b1
    w = b2 - np.dot(b2,b1)*b1
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1,v),w), np.dot(v,w))))

def rotation_matrix(axis, angle_deg):
    axis = np.asarray(axis,float)
    axis /= np.linalg.norm(axis)
    x,y,z = axis
    a = np.radians(angle_deg)
    c,s,C = np.cos(a),np.sin(a),1-np.cos(a)
    return np.array([
        [c+x*x*C, x*y*C-z*s, x*z*C+y*s],
        [y*x*C+z*s, c+y*y*C, y*z*C-x*s],
        [z*x*C-y*s, z*y*C+x*s, c+z*z*C],
    ])

def rotate_points(coords, indices, axis_a, axis_b, delta_deg):
    out = np.asarray(coords,float).copy()
    a,b = np.asarray(axis_a,float),np.asarray(axis_b,float)
    r = rotation_matrix(b-a, delta_deg)
    idx = np.asarray(sorted(set(indices)), dtype=int)
    out[idx] = (out[idx]-a) @ r.T + a
    return out
