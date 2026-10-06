#!/usr/bin/env python3
"""Exploratory per-DAPI condition plot from one or more measurements_selected.csv files."""

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mannwhitneyu


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('tables',nargs='+',type=Path)
    parser.add_argument('--output',type=Path,default=Path('condition_comparison.png'))
    args=parser.parse_args()
    groups=defaultdict(list)
    for path in args.tables:
        with path.open(newline='') as file:
            for row in csv.DictReader(file):
                if row.get('auto_selected',row.get('selected')) not in ('True','true','1','yes'):
                    continue
                groups[row['condition']].append(float(row['cen_bor_um']))
    if len(groups)<2:
        raise ValueError('Supply selected-object tables containing at least two conditions')
    names=sorted(groups)
    fig,ax=plt.subplots(figsize=(max(5,1.5*len(names)+2),4.8),constrained_layout=True)
    rng=np.random.default_rng(0)
    for i,name in enumerate(names):
        values=np.asarray(groups[name])
        ax.scatter(i+rng.uniform(-.14,.14,len(values)),values,s=17,alpha=.55)
        mean=float(values.mean())
        sd=float(values.std(ddof=1)) if len(values)>1 else 0.
        ax.errorbar(i,mean,yerr=sd,fmt='o',color='black',capsize=5,markersize=6)
        print(f'{name}: n DAPI bodies={len(values)}, mean={mean:.4f} µm, SD={sd:.4f} µm')
    ax.set_xticks(range(len(names)),names)
    ax.set_ylabel('3D centroid to LMN border (µm)')
    ax.set_title('Selected DAPI bodies: dots, mean ± SD')
    fig.savefig(args.output,dpi=160)
    if len(names)==2:
        test=mannwhitneyu(groups[names[0]],groups[names[1]],alternative='two-sided')
        print(f'Two-sided Mann–Whitney U={test.statistic:.4g}; p={test.pvalue:.4g}')
    else:
        print('Mann–Whitney test applies to two conditions; no pooled test calculated.')
    print('Per-body test is exploratory: DAPI bodies within stacks are not independent biological replicates.')


if __name__=='__main__':main()
