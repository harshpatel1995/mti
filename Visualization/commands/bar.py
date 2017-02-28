"""
Barchart of relative abundance for each organism in each sample
"""

import itertools

import seaborn as sb
import pandas as pd
from numpy import median

from .utils import parse_gra, parse_metadata, taxid_to_name, parse_sample_group_string


def run(options):
    # Get the data
    samples = parse_sample_group_string(options['<sample_group_string>'], False)
    samples['sample_organism'] = samples['sample'].map(str) + ', ' + samples['organism']


    # Get the metadata and add it to the samples
    metadata_filename = options['--meta']
    metadata = parse_metadata(metadata_filename)
    y_axis = options['--var']
    samples[y_axis] = samples['sample'].map(lambda x: int(metadata[x][y_axis]))
    samples['organism'] = samples['organism'].map(lambda t: taxid_to_name(t))

    # Draw the plots
    barPlot = sb.barplot(
            x=y_axis,
            y='rel_abund',
            hue='organism',
            capsize=.1,
            data=samples)

    barPlot.set(ylabel='relative abundance')
    barPlot.get_figure().savefig('bar.png')
