"""
Barchart of relative abundance for each organism in each sample
"""

import itertools

import seaborn as sb
import pandas as pd
from numpy import median

from .utils import parse_gra, parse_metadata, taxid_to_name


def run(options):
    # Get the data
    gra_filenames = options['<sample_group_string>'].split(',')
    samples = {g: parse_gra(g) for g in gra_filenames}
    data = [[
        {
            'sample': sample,
            'organism': taxid,
            'rel_abund': vals['rel_abund'],
            'error': vals['error']
            }
        for taxid, vals in d.items()
        ]
        for sample, d in samples.items()
        ]
    data = list(itertools.chain.from_iterable(data))
    samples = pd.DataFrame(data=data)
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
            estimator=median,
            capsize=.1,
            data=samples)

    barPlot.get_figure().savefig('bar.png')
