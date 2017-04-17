"""
Scatterplot of relative abundance vs. a metadata variable
"""


import itertools
import math
from textwrap import wrap

import pandas as pd
from scipy import stats
import seaborn as sb
from matplotlib import pyplot as plt

from .utils import parse_sample_list, set_common_ancestor, unique_filename


def run(options):
    # 2 options for scatterplots
    if options['<metadata>']:
        samples_y = parse_sample_list(options['<sample_group_string_v>'],
                options['<metadata>'])
        if options['--filter']:
            samples_y = set_common_ancestor(
                    samples_y, options['--filter'])
        # Get the metadata and add it to the samples
        x_axis = options['<var>']
        # Draw the plots
        fig_size = math.ceil(math.sqrt(len(samples_y['name'].unique())))
        fig = plt.figure(figsize=(fig_size*3, fig_size*3))
        for i, organism in enumerate(samples_y['name'].unique()):
            org_samples = samples_y.loc[samples_y['name'] == organism]
            ax = fig.add_subplot(fig_size, fig_size, i+1)
            name = organism
            ax.set_title('\n'.join(wrap(name, 20)))
            greatest_y = samples_y['rel_abund'].max()
            ax.set(ylim=(0, greatest_y+0.1))
            scatter = sb.regplot(x_axis, 'rel_abund', org_samples, ax=ax, ci=None)
        fig.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.5)
        fig.savefig('scatter.png', bbox_inches='tight')

    else:
        samples_y = parse_sample_list(options['<sample_group_string_v>'])
        samples_x = parse_sample_list(options['<sample_group_string_h>'])
        if options['--filter']:
            samples_y = set_common_ancestor(
                    samples_y, options['--filter'])
            samples_x = set_common_ancestor(
                    samples_x, options['--filter'])
        data = pd.DataFrame()
        data[options['<sample_group_string_v>']] = samples_y['rel_abund']
        data[options['<sample_group_string_h>']] = samples_x['rel_abund']
        # Draw the plots
        data = samples_x.merge(samples_y, on='name', how='outer').fillna(0)
        scatter = sb.JointGrid('rel_abund_y', 'rel_abund_x', data)
        corr = options['--correlation']
        if corr == 'spearman' or corr == 'spearmanr':
            corr_func = stats.spearmanr
        elif corr  == 'pearson' or corr == 'pearsonr':
            corr_func = stats.pearsonr
        else:
            corr_func = None
        scatter.plot(sb.regplot, sb.distplot, corr_func)
        if options['--annotate']:
            for row in data.iterrows():
                annotate_point(row)
        scatter.savefig('outputs/'+unique_filename()+'.png')

def annotate_point(row):
    ind = row[0]
    r = row[1]
    plt.gca().annotate(r['name'], xy=(r['rel_abund_x'], r['rel_abund_y']), fontsize=5, xytext=(2, 2), textcoords='offset points')
