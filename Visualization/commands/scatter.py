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

    if options['<metadata>']:
        # Get the metadata and add it to the samples
        metadata_filename = options['<metadata>']
        metadata = parse_metadata(metadata_filename)
        x_axis = options['<var>']
        samples_y[x_axis] = samples_y['sample'].map(lambda x: int(metadata[x][x_axis]))

        # Draw the plots
        fig_size = math.ceil(math.sqrt(len(samples_y['organism'].unique())))
        fig = plt.figure(figsize=(fig_size*3, fig_size*3))
        for i, organism in enumerate(samples_y['organism'].unique()):
            org_samples = samples_y.loc[samples_y['organism'] == organism]
            ax = fig.add_subplot(fig_size, fig_size, i+1)
            name = taxid_to_name(organism)
            ax.set_title('\n'.join(wrap(name, 20)))
            greatest_y = samples_y['rel_abund'].max()
            ax.set(ylim=(0, greatest_y+0.1))
            scatter = sb.regplot(x_axis, 'rel_abund', org_samples, ax=ax, ci=None)
        fig.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.5)
        fig.savefig('scatter.png', bbox_inches='tight')
    else:
        samples_x = parse_sample_group_string(options['<sample_group_string_h>'])
        data = pd.DataFrame()
        data[options['<sample_group_string_v>']] = samples_y['rel_abund']
        data[options['<sample_group_string_h>']] = samples_x['rel_abund']
        
        # Draw the plots
        scatter = sb.regplot(x=options['<sample_group_string_h>'], y=options['<sample_group_string_v>'], data=data, ci=none)
        greatest_x = data[options['<sample_group_string_h>']].max()
        greatest_y = data[options['<sample_group_string_v>']].max()
        ax_max = max(greatest_y, greatest_x) + 0.01
        plt.ylim(0, ax_max)
        plt.xlim(0, ax_max)
        scatter.get_figure().savefig('scatter.png')
