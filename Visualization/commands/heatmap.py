"""
A heatmap of taxons to samples
"""


from .utils import parse_sample_list, pivot_on_sample_and_name, set_common_ancestor, unique_filename
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
import itertools


def run(options):
    samples = parse_sample_list(options['<sample_group_string>'])
    plt.figure(figsize=(len(samples)/1.2,len(samples)/1.2))
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

    yticks = samples.index.values
    xticks = samples.columns.values
    heatmap = sb.heatmap(
            samples,
            annot = True,
            annot_kws = {
                'size': 10,
                'alpha': 0.8 
                },
            xticklabels = xticks,
            yticklabels = yticks,
            square = True,
            vmin=vmin,
            vmax=vmax)

    # Style
    plt.yticks(rotation=0)
    plt.xticks(rotation=30, ha='right')
    plt.axes().set_title('Relative Abundance')

    heatmap.get_figure().savefig('outputs/'+unique_filename()+'.png', bbox_inches='tight')
