"""
Scatterplot of relative abundance vs. a metadata variable
"""


import itertools
import math
from textwrap import wrap

import pandas as pd
import seaborn as sb
from matplotlib import pyplot as plt

from .utils import parse_gra, taxid_to_name, parse_metadata, parse_sample_group_string


def run(options):
    # Get the data
    samples_y = parse_sample_group_string(options['<sample_group_string_v>'])

    if (options['<metadata>']):
        # Get the metadata and add it to the samples
        metadata_filename = options['<metadata>']
        metadata = parse_metadata(metadata_filename)
        x_axis = options['<var>']
        samples[x_axis] = samples['sample'].map(lambda x: int(metadata[x][x_axis]))

        # Draw the plots
        fig = plt.figure()
        fig_size = math.ceil(math.sqrt(len(samples['organism'].unique())))
        for i, organism in enumerate(samples['organism'].unique()):
            org_samples = samples.loc[samples['organism'] == organism]
            ax = fig.add_subplot(fig_size, fig_size, i+1)
            ax.set_title('\n'.join(wrap(organism, 60)))
            ax.set(ylim=(-0.1, 1))
            scatter = sb.regplot(x_axis, 'rel_abund', org_samples, ax=ax, ci=None)
        fig.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.5)
        fig.savefig('scatter.png', bbox_inches='tight')
    else:
        samples_x = parse_sample_group_string(options['<sample_group_string_h>'])
        data = pd.DataFrame()
        data[options['<sample_group_string_v>']] = samples_y['rel_abund']
        data[options['<sample_group_string_h>']] = samples_x['rel_abund']
        
        # Draw the plots
        scatter = sb.regplot(x=options['<sample_group_string_h>'], y=options['<sample_group_string_v>'], data=data, ci=None)
        greatest_y = data[options['<sample_group_string_v>']].max()
        plt.ylim(0, greatest_y + 0.01)
        scatter.get_figure().savefig('scatter.png')
