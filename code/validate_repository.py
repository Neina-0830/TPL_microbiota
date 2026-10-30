#!/usr/bin/env python3
"""Check source syntax, primary input checksums and binary process removals."""
import ast
import csv
import hashlib
import sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code/07_sGLV'))
from sGLV_core import calibrated_parameters, process_removal_parameters, SCENARIOS

def main():
    for p in (ROOT/'code').rglob('*.py'): ast.parse(p.read_text(),filename=str(p))
    with (ROOT/'DATA_MANIFEST.csv').open() as f:
        for row in csv.DictReader(f):
            p=ROOT/row['Path']
            assert hashlib.sha256(p.read_bytes()).hexdigest()==row['SHA256'],p
    inp={'initial':np.array([[.1],[.2],[.3]]),
         'immigration_weight':np.array([[.4],[.7],[1.]]),
         'growth_weight':np.array([[1.],[3.],[5.]]),
         'interaction':np.array([[-1.,.2,-.1],[.1,-1.,.2],[.4,.3,-1.]])}
    cfg={'warmed_noise_multiplier':.9}
    base=calibrated_parameters(inp,cfg,'warming',.02,.1,.3,.05)
    changed={'immigration_removed':'immigration_rate','growth_removed':'growth_rate',
             'interaction_removed':'interaction_matrix','noise_removed':'noise_linear'}
    for s in SCENARIOS:
        p=process_removal_parameters(inp,cfg,'warming',.02,.1,.3,.05,s)
        for key in base:
            if key!=changed.get(s): assert np.array_equal(p[key],base[key]),(s,key)
        if s=='growth_removed': assert np.allclose(p['growth_rate'],np.mean(base['growth_rate']))
        if s=='interaction_removed': assert np.array_equal(p['interaction_matrix'],-np.eye(3))
        if s=='immigration_removed': assert np.all(p['immigration_rate']==0)
        if s=='noise_removed': assert p['noise_linear']==0
    print('PASS: Python syntax, data checksums and all process-removal identities.')

if __name__=='__main__': main()
