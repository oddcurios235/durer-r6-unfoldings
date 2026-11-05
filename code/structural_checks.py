#!/usr/bin/env python3
"""Exact structural checks for the disjoint-projection R6 witness."""
from fractions import Fraction as Q
from pathlib import Path
import json,itertools
import prismatoid_checker as C
import r6_main_checks
ROOT=Path(__file__).resolve().parents[1]
def run(check=False):
    data=r6_main_checks.cases()[-1];V=C.get_vertices(data);F=data['faces'];E=C.verify_hull(V,F)
    deg=[0]*len(V)
    for a,b in E:deg[a]+=1;deg[b]+=1
    normals=[C.cross(C.vec(V[f[0]],V[f[1]]),C.vec(V[f[0]],V[f[2]])) for f in F]
    parallel=[[i,j] for i,j in itertools.combinations(range(len(F)),2) if C.cross(normals[i],normals[j])==(0,0,0)]
    assert parallel==[[0,len(F)-1]]
    disjoint={str(i):next(j for j in range(len(F)) if not set(F[i])&set(F[j])) for i in range(len(F))}
    assert deg[0]>4 and max(deg[4:])==4
    c,s,b,l,t,h=[Q(data['parameters'][x]) for x in ('c','s','b','l','t','h')]
    assert c> s
    # Physical x/c equals oblique coordinate sum.
    lower_max=max(v[0]/c for v in V[:4]);upper_min=min(v[0]/c for v in V[4:])
    assert lower_max==3 and upper_min==Q(19,4)
    out=dict(name=data['name'],degrees=deg,parallel_facet_pairs=parallel,
             disjoint_facet_witnesses=disjoint,lower_label_sum_max=str(lower_max),upper_label_sum_min=str(upper_min),
             distinguished_lower_corner_dot=str(c*c-s*s))
    path=ROOT/'certificates/structural_witness.json'
    if check:assert json.loads(path.read_text())==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':
    import sys;run('--check' in sys.argv)
