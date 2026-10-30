#!/usr/bin/env python3
"""Simulate the calibrated model and four single-process removals."""
import argparse
import pandas as pd
from sGLV_core import (CONDITIONS, SEEDS, SCENARIOS, load_group_config,
    validate_group_names, load_condition_inputs, process_removal_parameters,
    process_filename, simulate_sglv, write_simulation, project_root, work_root)

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--tpl',choices=['tTTPL','pTTPL'],required=True)
    p.add_argument('--group',nargs='+',default=['all'])
    p.add_argument('--condition',choices=['warming','unwarming','both'],default='both')
    p.add_argument('--overwrite',action='store_true')
    p.add_argument('--dry-run',action='store_true')
    a = p.parse_args();root=project_root();work=work_root(root);configs=load_group_config()
    for group in validate_group_names(a.group,configs):
        file=work/'Calibration'/a.tpl/f'{group}_top10_parameters.csv'
        top=pd.read_csv(file)
        if len(top)!=10 or top['Rank'].duplicated().any():
            raise ValueError(f'{file}: ten distinct parameter sets are required')
        conditions=CONDITIONS if a.condition=='both' else (a.condition,)
        for condition in conditions:
            cache={}
            for row in top.sort_values('Rank').itertuples():
                vals=(row.Immigration,row.Interaction,row.Growth,row.Noise)
                for seed in SEEDS:
                    for scenario in SCENARIOS:
                        out=(work/'get_timeseries'/group/a.tpl/'process_removal_top10'/scenario/condition/
                             process_filename(row.Rank,*vals[:3],seed,vals[3],scenario))
                        if out.exists() and not a.overwrite: continue
                        if a.dry_run: print(out);continue
                        if row.Interaction not in cache:
                            cache[row.Interaction]=load_condition_inputs(root,configs[group],condition,row.Interaction)
                        inp=cache[row.Interaction]
                        params=process_removal_parameters(inp,configs[group],condition,*vals,scenario)
                        result=simulate_sglv(params,inp['taxa'],seed,configs[group]['resample_depth'])
                        write_simulation(result,out,overwrite=a.overwrite)

if __name__=='__main__': main()
