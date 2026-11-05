#!/usr/bin/env python3
"""Principal finite prismatoid unfolding tests using the bundled checker.
These rational-input tests do not prove the universal assertion; for the
analytic proof, see the manuscript (main.tex / main.pdf).
"""
from fractions import Fraction as Q
from pathlib import Path
import json,argparse
import prismatoid_checker as C

def build(name,c,s,b,l,t,h,chain):
    c,s,b,l,t,h=map(Q,(c,s,b,l,t,h)); chain=[tuple(map(Q,z)) for z in chain]
    assert c*c+s*s==1 and 0<s<=c and b>0 and l>0 and h>0
    assert t>=max(b,l)
    assert chain[0]==(0,t) and chain[-1]==(t,0)
    assert all(x<X and y>Y for (x,y),(X,Y) in zip(chain,chain[1:]))
    V=[(0,0,0),(b*c,-b*s,0),((b+l)*c,(l-b)*s,0),(l*c,l*s,0)]
    V += [(c*(x+y),s*(y-x),h) for x,y in chain]+[(2*t*c,0,h)]
    u,v,w=4,len(V)-2,len(V)-1
    F=[[3,2,1,0],[0,1,v],[0,u,3],[1,2,w,v],[2,3,u,w]]
    F += [[0,i+1,i] for i in range(u,v)]+[list(range(u,w+1))]
    T=[[0,i] for i in range(u,v+1)]+[[0,1],[0,3],[2,w]]
    T += [[1,2]] if l>=b else [[2,3]]
    return dict(name=name,vertices=[list(map(str,z)) for z in V],faces=F,cut_edges=T,
                parameters=dict(zip(('c','s','b','l','t','h'),map(str,(c,s,b,l,t,h)))),
                chain=[list(map(str,z)) for z in chain],selection='N' if l>=b else 'E')

def cases():
    out=[]
    for c,s in [(Q(9999,10001),Q(200,10001)),(Q(99,101),Q(20,101)),(Q(12,13),Q(5,13)),(Q(4,5),Q(3,5)),(Q(21,29),Q(20,29))]:
      for b,l in [(Q(1),Q(1)),(Q(1,7),Q(3)),(Q(3),Q(1,7))]:
       for t in [max(b,l),max(b,l)+min(b,l)/(2*(4*c*c-1)),max(b,l)+min(b,l)/(4*c*c-1)]:
        for k in [1,2,7]:
         # Amplitude approaching t still yields a strictly convex monotone chain.
         for amp in ([Q(0)] if k==1 else [Q(1,5),Q(19,20)]):
          chain=[(t*Q(j,k),t*(1-Q(j,k))-amp*t*Q(j,k)*(1-Q(j,k))) for j in range(k+1)]
          for h in [Q(1,100000),Q(2),Q(1000)]:
           out.append(build('r6_'+str(len(out)),c,s,b,l,t,h,chain))
    # Explicit very acute disjoint-projection high-vertex case.
    c,s=Q(99,101),Q(20,101);b,l,t=Q(1),Q(2),Q(5);k=20
    chain=[(t*Q(j,k),t*(1-Q(j,k))-Q(j,k)*(1-Q(j,k))) for j in range(k+1)]
    out.append(build('r6_long_disjoint',c,s,b,l,t,Q(1,10000),chain))
    return out

def run(check=False):
    p=Path(__file__).resolve().parents[1]/'certificates/r6_main.json'
    old=json.loads(p.read_text()) if check else None
    data=cases()
    if check:assert data==old['cases']
    records=[]
    for x in data:
        r=C.verify_case(x);records.append(r)
        print(x['name'],x['selection'],r['vertices'],r['all_face_pairs_certified'],flush=True)
    result=dict(cases=data,verification=records,totals=dict(cases=len(records),pairs=sum(r['all_face_pairs_certified'] for r in records),max_vertices=max(r['vertices'] for r in records)))
    if check:assert result==old
    else:p.write_text(json.dumps(result,indent=2)+'\n')
    print('TOTAL',result['totals'])
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--check',action='store_true');a=ap.parse_args();run(a.check)
