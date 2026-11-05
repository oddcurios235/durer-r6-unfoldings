#!/usr/bin/env python3
"""Additional finite R6 tests for both named cut trees and asymmetric upper
chains, checked using the bundled geometric verifier.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse, copy, json
import r6_main_checks as R
import prismatoid_checker as C

def other_tree(case):
    out=copy.deepcopy(case)
    old,new=([1,2],[2,3]) if out['selection']=='N' else ([2,3],[1,2])
    assert old in out['cut_edges'] and new not in out['cut_edges']
    out['cut_edges'].remove(old); out['cut_edges'].append(new)
    out['selection']='E' if out['selection']=='N' else 'N'
    out['name'] += '_other'
    return out

def cases():
    out=[other_tree(z) for z in R.cases()]
    for c,s in [(Q(9999,10001),Q(200,10001)),(Q(12,13),Q(5,13)),(Q(21,29),Q(20,29))]:
      for b,l in [(Q(1),Q(1)),(Q(1,7),Q(3)),(Q(3),Q(1,7))]:
       for h in [Q(1,100000),Q(1000)]:
        for beta in [Q(-1,3),Q(1,3)]:
         t=max(b,l); k=7
         ts=[Q(j*j,k*k) for j in range(k+1)]
         chain=[(t*x,t*(1-x)-t*Q(1,2)*x*(1-x)*(1+beta*x)) for x in ts]
         z=R.build('r6_asymmetric_'+str(len(out)),c,s,b,l,t,h,chain)
         out.extend([z,other_tree(z)])
    return out

def run(check=False):
    p=Path(__file__).resolve().parents[1]/'certificates/r6_additional.json'
    old=json.loads(p.read_text()) if check else None
    data=cases()
    if check: assert data==old['cases']
    records=[]
    for z in data:
        r=C.verify_case(z);records.append(r)
        print(z['name'],z['selection'],r['all_face_pairs_certified'],flush=True)
    result=dict(cases=data,verification=records,totals=dict(cases=len(records),pairs=sum(r['all_face_pairs_certified'] for r in records),max_vertices=max(r['vertices'] for r in records)))
    if check: assert result==old
    else: p.write_text(json.dumps(result,indent=2)+'\n')
    print('ADDITIONAL TOTAL',result['totals'])
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--check',action='store_true');run(a.parse_args().check)
