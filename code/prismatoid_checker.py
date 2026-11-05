#!/usr/bin/env python3
"""Exact rational-geometry checker for prismatoid unfolding tests.

For supplied rational polyhedra and original-edge cut trees, check the maximal
facets, spanning-tree structure, rigid planar developments and pairwise
separation using rational predicates and outward rational interval bounds.
The functions are imported by the bundled finite test programs.
"""
import collections, itertools, math
from fractions import Fraction as Q

# Represent a number by an interval with exact integer
# endpoints x.lo/M <= x <= x.hi/M; outward integer rounding after products.
M=1<<190
class Interval:
    __slots__=('lo','hi')
    def __init__(self,lo,hi=None):
        self.lo=int(lo);self.hi=int(hi if hi is not None else lo)
        assert self.lo<=self.hi
    @staticmethod
    def exact(v):
        v=Q(v)
        lo=(v.numerator*M)//v.denominator
        hi=-((-v.numerator*M)//v.denominator)
        return Interval(lo,hi)
    @staticmethod
    def sqrt(v):
        v=Q(v);assert v>=0
        r=math.isqrt((v.numerator*M*M)//v.denominator)
        return Interval(r,r+(r*r*v.denominator!=v.numerator*M*M))
    def __add__(self,r):
        r=asiv(r);return Interval(self.lo+r.lo,self.hi+r.hi)
    __radd__=__add__
    def __neg__(self):return Interval(-self.hi,-self.lo)
    def __sub__(self,r):return self+-asiv(r)
    def __rsub__(self,r):return asiv(r)+-self
    def __mul__(self,r):
        r=asiv(r);z=(self.lo*r.lo,self.lo*r.hi,self.hi*r.lo,self.hi*r.hi)
        return Interval(min(z)//M, -((-max(z))//M))
    __rmul__=__mul__
    def __repr__(self):return f'I[{float(self.lo/M):.13g},{float(self.hi/M):.13g}]'
def asiv(x):return x if isinstance(x,Interval) else Interval.exact(x)

def sub(a,b):return tuple(x-y for x,y in zip(a,b))
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def cross(a,b):return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def det(a,b):return a[0]*b[1]-a[1]*b[0]
def left(v):return (-v[1],v[0])
def sq3(a):return dot(a,a)
def sqdist(a,b):return sq3(sub(a,b))
def vec(a,b):return sub(b,a)
def get_vertices(data):return [tuple(Q(y) for y in x) for x in data['vertices']]

class UnionFind:
    def __init__(self,n):self.p=list(range(n));self.count=n
    def rep(self,i):
        while self.p[i]!=i:
            self.p[i]=self.p[self.p[i]];i=self.p[i]
        return i
    def join(self,a,b):
        a=self.rep(a);b=self.rep(b)
        if a==b:return False
        self.p[b]=a;self.count-=1;return True

def verify_hull(vertices,faces):
    """Enumerate all 3D supporting planes through vertex triples, no supplied facet assumptions."""
    n=len(vertices)
    found=set()
    for i,j,k in itertools.combinations(range(n),3):
        normal=cross(vec(vertices[i],vertices[j]),vec(vertices[i],vertices[k]))
        if normal==(0,0,0):continue
        values=[dot(normal,vec(vertices[i],p)) for p in vertices]
        if min(values)>=0 or max(values)<=0:
            found.add(frozenset(idx for idx,d in enumerate(values) if d==0))
    asserted={frozenset(f) for f in faces}
    assert len(asserted)==len(faces),'duplicate facets'
    assert found==asserted, f'wrong maximal face list: hull-only {found-asserted}, provided-only {asserted-found}'
    edge_faces=collections.defaultdict(list)
    for idx,f in enumerate(faces):
        assert len(f)>=3 and len(f)==len(set(f))
        nrm=cross(vec(vertices[f[0]],vertices[f[1]]),vec(vertices[f[0]],vertices[f[2]]))
        assert nrm!=(0,0,0)
        assert all(dot(nrm,vec(vertices[f[0]],v))<=0 for v in vertices),'facet orientation outward incorrect'
        for i,j in zip(f,f[1:]+f[:1]):
            assert i!=j
            e=vec(vertices[i],vertices[j])
            # Strictly convex polygon, each edge line supports every other facet vertex.
            assert all(dot(cross(e,vec(vertices[i],vertices[k])),nrm)>0 for k in f if k!=i and k!=j),'invalid facet cyclic order'
            edge_faces[tuple(sorted((i,j)))].append(idx)
    assert all(len(v)==2 for v in edge_faces.values()),'not a 2-manifold polyhedron'
    assert n-len(edge_faces)+len(faces)==2,'Euler characteristic failure'
    return edge_faces

def verify_tree(n,cuts,edge_faces):
    assert all(e in edge_faces for e in cuts),'cut is not genuine edge'
    assert len(cuts)==n-1
    U=UnionFind(n)
    assert all(U.join(*e) for e in cuts),'cycles in cut tree'
    assert U.count==1

def unfold_all(vertices,faces,edge_faces,cuts):
    """Distance-only trilateration with positive oriented face side."""
    m=len(faces)
    neighbors=[[] for _ in faces]
    for edge,side in edge_faces.items():
        if edge not in cuts:
            i,j=side
            neighbors[i].append((j,edge))
            neighbors[j].append((i,edge))
    # Verify the connected uncut dual is a tree.
    assert sum(map(len,neighbors))==2*(m-1)
    cc=UnionFind(m)
    for i in range(m):
        for j,_ in neighbors[i]:
            if i<j:assert cc.join(i,j),'uncut dual cycle'
    assert cc.count==1
    # Corner-equality relation for two faces glued along a path of intact hinges.
    corner_ids={(i,v):j for j,(i,v) in enumerate((i,v) for i,f in enumerate(faces) for v in f)}
    corners=UnionFind(len(corner_ids))
    for edge,side in edge_faces.items():
        if edge not in cuts:
            i,j=side
            for v in edge: corners.join(corner_ids[(i,v)],corner_ids[(j,v)])
    developed=[None]*m
    f0=faces[0]
    v0,v1=f0[:2]
    # First polygon: its first edge placed along the x axis with original metric length.
    developed[0]=_place_face(vertices,f0,v0,v1,(asiv(0),asiv(0)),(Interval.sqrt(sqdist(vertices[v0],vertices[v1])),asiv(0)))
    todo=[0]
    while todo:
        face_idx=todo.pop()
        for j,hinge in neighbors[face_idx]:
            if developed[j] is not None:continue
            a,b=hinge
            f=faces[j]
            # Order hinge in *child* boundary orientation, so new child lies to left.
            if f[(f.index(a)+1)%len(f)]!=b:a,b=b,a
            developed[j]=_place_face(vertices,f,a,b,developed[face_idx][a],developed[face_idx][b])
            todo.append(j)
    assert all(d is not None for d in developed)
    return developed,corner_ids,corners

def _place_face(vertices,face,a,b,pa,pb):
    q={a:pa,b:pb}
    w=vec(pa,pb)
    d2=sqdist(vertices[a],vertices[b]);assert d2>0
    for i in face:
        if i in q:continue
        # Derive coordinates from *three pairwise squared lengths*, not 3D dot.
        aa=sqdist(vertices[a],vertices[i])
        bb=sqdist(vertices[b],vertices[i])
        coeff=(aa+d2-bb)/(2*d2)
        rad=aa/d2-coeff*coeff
        assert rad>0,'collinear facet vertices'
        perp=Interval.sqrt(rad)
        q[i]=(pa[0]+coeff*w[0]-perp*w[1],pa[1]+coeff*w[1]+perp*w[0])
    return q

def signed_area_cross(pa,pb,pz):return det(vec(pa,pb),vec(pa,pz))

def axis_witness(faces,positions,pair,corner_ids,corners):
    i,j=pair
    for a_idx,b_idx in ((i,j),(j,i)):
        f=faces[a_idx]
        for x,y in zip(f,f[1:]+f[:1]):
            side=[]
            for v in faces[b_idx]:
                if corners.rep(corner_ids[(a_idx,x)])==corners.rep(corner_ids[(b_idx,v)]) or corners.rep(corner_ids[(a_idx,y)])==corners.rep(corner_ids[(b_idx,v)]):
                    # Topologically identical boundary corner, rigorously zero by original hinges.
                    continue
                crossv=signed_area_cross(positions[a_idx][x],positions[a_idx][y],positions[b_idx][v])
                side.append(crossv)
            if all(t.hi<=0 for t in side):return (a_idx,b_idx,x,y,side)
    raise ValueError(f'no outward rational-interval convex separator for face pair {pair}')

def verify_case(data,print_pairs=False):
    vertices=get_vertices(data);faces=[list(f) for f in data['faces']]
    cut={tuple(sorted(e)) for e in data['cut_edges']}
    assert len(cut)==len(data['cut_edges'])
    graph=verify_hull(vertices,faces)
    verify_tree(len(vertices),cut,graph)
    positions,corner_ids,corners=unfold_all(vertices,faces,graph,cut)
    result=dict(name=data.get('name'),vertices=len(vertices),faces=len(faces),edges=len(graph),cuts=len(cut),uncut_hinges=len(graph)-len(cut))
    n=0
    for pair in itertools.combinations(range(len(faces)),2):
        wit=axis_witness(faces,positions,pair,corner_ids,corners)
        n+=1
        if print_pairs:print(' pair',pair,'axis',wit[:4])
    result['all_face_pairs_certified']=n
    return result
