"""
A heatmap of taxons to samples
"""


from .utils import parse_sample_list, pivot_on_sample_and_name, set_common_ancestor, unique_filename
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
import itertools
import math
import numpy as np

def run(options):
    samples = parse_sample_list(options['<sample_group_string>'])
    q1 = q3 = iqr = 0.0
    plt.figure(figsize=(len(samples)/1.2,len(samples)/1.2))
    if options['--outliers']:
        q1 = samples["rel_abund"].quantile(0.25)
        q3 = samples["rel_abund"].quantile(0.75)
        iqr = q3-q1
    if options['--filter']:
        samples = set_common_ancestor(
                samples, options['--filter'])
    samples = pivot_on_sample_and_name(samples)

    vmin = None
    vmax = None
    if options['--scale-min']:
        vmin = float(options['--scale-min'])
    if options['--scale-max']:
        vmax = float(options['--scale-max'])
    if options['--outliers']:
        # check for outliers and append to list
        outliershigh=[]
        outlierslow=[]
        values=[]
        for thing in samples.values:
            for val in thing:
                if (val > q3 + 1.5*iqr):
                    outliershigh.append(val)
                elif (val < q1-1.5*iqr):
                    outlierslow.append(val)
                else:
                    values.append(val)
        vmax=max(values)
        vmin=min(values)
        maxoutlier = max(outliershigh)

        yticks = samples.index.values
        xticks = samples.columns.values
        heatmap = sb.heatmap(
                samples,
                annot = True,
                annot_kws = {
                    'size': 10,
                    'alpha': 1.0 
                    },
                xticklabels = xticks,
                yticklabels = yticks,
                square = True,
                vmin=vmin,
                vmax=vmax,
                cbar = False,
                cmap = 'Reds')
        
        heatmap = sb.heatmap(
                samples,
                annot = True,
                annot_kws = {
                    'size': 10,
                    'alpha': 0.0 
                    },
                xticklabels = xticks,
                yticklabels = yticks,
                square = True,
                vmin=vmin,
                vmax=vmax,
                mask = samples.values > (q3 + 1.5*iqr),
                cmap = 'Blues')
    # Style
    plt.yticks(rotation=0)
    plt.xticks(rotation=30, ha='right')
    plt.axes().set_title('Relative Abundance')

    heatmap.get_figure().savefig('outputs/'+unique_filename()+'.png', bbox_inches='tight')