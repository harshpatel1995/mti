"""
Barchart of relative abundance for each organism in each sample
"""

import itertools

import seaborn as sb
import pandas as pd
from numpy import median
import matplotlib.pyplot as plt

from .utils import parse_sample_list


def run(options):
    # Get the data
    samples = parse_sample_list(options['<sample_group_string>'], 
            options['--meta'])
    x_axis = options['--var']

    # Draw the plots
    barPlot = sb.factorplot(
            x=x_axis,
            y='rel_abund',
            hue='name',
            data=samples,
            size=6,
            kind='bar',
            legend=False,
            edgecolor='.2')

    barPlot.despine(left=True)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    barPlot.set_ylabels('relative abundance')
    barPlot.savefig('bar.png')
