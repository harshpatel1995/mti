"""
Barchart of relative abundance for each organism in each sample
"""

import itertools

import seaborn as sb
import pandas as pd
from numpy import median
import matplotlib.pyplot as plt

from .utils import parse_sample_list_string, parse_metadata, taxid_to_name


def run(options):
    # Get the data
    samples = parse_sample_list_string(options['<sample_group_string>'], False)
    samples['sample_organism'] = samples['sample'].map(str) + ', ' + samples['organism']


    # Get the metadata and add it to the samples
    metadata_filename = options['--meta']
    metadata = parse_metadata(metadata_filename)
    y_axis = options['--var']
    samples[y_axis] = samples['sample'].map(lambda x: metadata[x][y_axis])
    samples['organism'] = samples['organism'].map(lambda t: taxid_to_name(t))

    # Draw the plots
    barPlot = sb.factorplot(
            x=y_axis,
            y='rel_abund',
            hue='organism',
            data=samples,
            size=6,
            kind='bar',
            legend=False)

    barPlot.despine(left=True)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    barPlot.set_ylabels('relative abundance')
    barPlot.savefig('bar.png')
