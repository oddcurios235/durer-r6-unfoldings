#!/usr/bin/env python3
"""Exact algebraic identities and rational-geometry checks for this manuscript.
These finite tests complement the analytic arguments in main.tex.
"""
from fractions import Fraction as Q
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]

def add(z,w):return (z[0]+w[0],z[1]+w[1])
def neg(z):return (-z[0],-z[1])
def sub(z,w):return add(z,neg(w))
def dot(z,w):return z[0]*w[0]+z[1]*w[1]

def polynomial_identity():
    # Variables k,d,x. Reduce d^2=1-k^2.
    def plus(a,b):
        r=a.copy()
        for e,v in b.items():r[e]=r.get(e,0)+v
        return {e:v for e,v in r.items() if v}
    def times(a,b):
        r={}
        for ea,va in a.items():
          for eb,vb in b.items():
            e=tuple(x+y for x,y in zip(ea,eb));r[e]=r.get(e,0)+va*vb
        return {e:v for e,v in r.items() if v}
    def reduce(a):
        r={}
        for (i,j,m),v in a.items():
            if j>=2:
                r=plus(r,reduce({(i,j-2,m):v,(i+2,j-2,m):-v}))
            else:r=plus(r,{(i,j,m):v})
        return r
    one={(0,0,0):1};k={(1,0,0):1};d={(0,1,0):1};x={(0,0,1):1}
    k2=times(k,k);x2=times(x,x)
    A=plus(times(k,plus(x2,one)),times(d,x))
    left=plus(times(A,A),{e:-v for e,v in plus(k2,x2).items()})
    right=times(times(k,x),times(plus(x2,one),plus(times(k,x),{e:2*v for e,v in d.items()})))
    assert reduce(left)==reduce(right)
    return 'coefficient identity modulo d^2+k^2=1 verified'

def verify_rational_geometry_identities():
    # Exact examples with Q represented as sqrt(k^2 a^2+h^2)+ka;
    # choose Pythagorean rational triples with h/a rational.
    tests=0
    for c,s in [(Q(99,101),Q(20,101)),(Q(4,5),Q(3,5))]:
      k=2*c*s;d=c*c-s*s
      for a in [Q(0),Q(1,7),Q(3)]:
       for l in [Q(1),Q(4,3)]:
        t=a+l
        for h in ([Q(1),Q(3)] if a==0 else [k*a*Q(3,4),k*a*Q(4,3)]):
            root=h if a==0 else (k*a*Q(5,4) if h==k*a*Q(3,4) else k*a*Q(5,3))
            assert root*root==k*k*a*a+h*h
            zshift=root+k*a
            U=(t*c,t*s);Z=(l*c+zshift*s,l*s+zshift*c)
            e=sub(U,Z);Je=(-e[1],e[0]);den=a*a+h*h
            ON=sub(Z,((l*(a*e[0]+h*Je[0]))/den,(l*(a*e[1]+h*Je[1]))/den))
            g=(h*c-t*s,h*s+t*c)
            assert dot(g,Z)-h*t==zshift*(k*h+d*t)-h*a
            assert den*(dot(g,ON)-h*t)==(h*h+t*t)*(zshift*(k*h+d*a)-h*a)
            assert dot(g,Z)>h*t and dot(g,ON)>h*t
            assert dot(e,e)==den
            assert dot(sub(ON,Z),sub(ON,Z))==l*l
            tests+=1
    return tests

def run(check=False):
    out=dict(scalar_identity=polynomial_identity(),geometry_identity_tests=verify_rational_geometry_identities())
    p=ROOT/'certificates/algebra_checks.json'
    if check:assert json.loads(p.read_text())==out
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':
    import sys;run('--check' in sys.argv)
