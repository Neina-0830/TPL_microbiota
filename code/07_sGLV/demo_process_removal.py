#!/usr/bin/env python3
"""Small synthetic example of process removal; does not reproduce field results."""
import argparse
import time
from pathlib import Path
import numpy as np
import pandas as pd
from sGLV_core import SCENARIOS, process_removal_parameters, simulate_sglv

def exponent(x):
    # Same eleven-snapshot window and relative-abundance normalization as R.
    a=x.iloc[:11].to_numpy(dtype=float)
    a=a/a.sum(axis=1,keepdims=True)
    mu=a.mean(axis=0);va=a.var(axis=0,ddof=1);keep=(mu>0)&(va>0)
    return np.polyfit(np.log10(mu[keep]),np.log10(va[keep]),1)[0]

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',default='results/demo_process_removal')
    a=p.parse_args();start=time.perf_counter();out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    rng=np.random.RandomState(731);n=20
    A=rng.normal(0,.005,(n,n));np.fill_diagonal(A,-1)
    inp={'initial':np.linspace(.01,.08,n).reshape(-1,1),
         'immigration_weight':np.linspace(.2,1,n).reshape(-1,1),
         'growth_weight':np.linspace(.5,1.5,n).reshape(-1,1),'interaction':A}
    cfg={'warmed_noise_multiplier':.95}
    rows=[]
    for rank in range(1,11):
        for cond in ('warming','unwarming'):
            d=.001*rank;g=.1*(1.1 if cond=='warming' else 1)
            for seed in range(3):
                for scenario in SCENARIOS:
                    params=process_removal_parameters(inp,cfg,cond,d,.005,g,.02,scenario)
                    ts=simulate_sglv(params,[f'T{i}' for i in range(n)],seed,1000,dt=.1)
                    rows.append(dict(Group='Synthetic',TPL='tTTPL',Rank=rank,Condition=cond,
                                     Seed=seed,Scenario=scenario,Beta=exponent(ts)))
    x=pd.DataFrame(rows);x.to_csv(out/'TPL_estimates.csv',index=False)
    avg=x.groupby(['Rank','Condition','Scenario']).Beta.mean().unstack('Scenario')
    effects=pd.DataFrame({s:avg['full']-avg[s] for s in SCENARIOS if s!='full'})
    importance=effects.abs().div(effects.abs().sum(axis=1),axis=0)
    assert np.allclose(importance.sum(axis=1),1)
    importance.to_csv(out/'expected_importance.csv')
    print(f'Wrote {len(x)} TPL estimates and {len(importance)} importance panels in {time.perf_counter()-start:.3f} seconds.')

if __name__=='__main__': main()
