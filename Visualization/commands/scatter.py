"""
Scatterplot of relative abundance vs. a metadata variable
"""


import itertools
import math
from textwrap import wrap

import pandas as pd
import seaborn as sb
from matplotlib import pyplot as plt

from .utils import parse_gra, taxid_to_name, parse_metadata


def run(options):
    # Get the data
    gra_filenames = options['<sample_group_string_v>'].split(',')
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

    # Get the metadata and add it to the samples
    metadata_filename = options['<metadata>']
    metadata = parse_metadata(metadata_filename)
    y_axis = options['<var>']
    samples[y_axis] = samples['sample'].map(lambda x: int(metadata[x][y_axis]))

    # Draw the plots
    fig = plt.figure()
    fig_size = math.ceil(math.sqrt(len(samples['organism'].unique())))
    for i, organism in enumerate(samples['organism'].unique()):
        org_samples = samples.loc[samples['organism'] == organism]
        ax = fig.add_subplot(fig_size, fig_size, i+1)
        ax.set_title('\n'.join(wrap(organism, 60)))
        ax.set(ylim=(-0.1, 1))
        scatter = sb.regplot(y_axis, 'rel_abund', org_samples, ax=ax)
    fig.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.5)
    fig.savefig('scatter.png', bbox_inches='tight')
