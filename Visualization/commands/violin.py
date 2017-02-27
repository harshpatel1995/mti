"""
Violin plot of organisms in each sample
"""


import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
import itertools

from .utils import parse_gra, taxid_to_name, parse_sample_group_string


def run(options):
    # Get the data
    samples = parse_sample_group_string(options['<sample_group_string>'], False)

    # Draw the plot
    yticks = [taxid_to_name(o) for o in samples['organism']]
    violinPlot = sb.violinplot(
            x='rel_abund', 
            y='organism', 
            data=samples,
            orient='h')

    # Style
    violinPlot.set_yticklabels(yticks)

    violinPlot.get_figure().savefig('violin.png', bbox_inches='tight')
