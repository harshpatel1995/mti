"""
A heatmap of taxons to samples
"""


from .utils import parse_gra, taxid_to_name
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
import itertools


def run(options):
    gra_filenames = options['<sample_group_string>'].split(',')
    samples = {g: parse_gra(g) for g in gra_filenames}
    l = [[{
        'sample': sample,
        'organism': taxid, 
        'rel_abund': vals['rel_abund'],
        'error': vals['error']
        }
        for taxid, vals in d.items()] for sample, d in samples.items()]
    data = list(itertools.chain.from_iterable(l))
    df = pd.DataFrame(data=data)
    samples = df.pivot ('sample', 'organism', 'rel_abund')

    vmin = None
    vmax = None
    if options['--scale-min']:
        vmin = float(options['--scale-min'])
    if options['--scale-max']:
        vmax = float(options['--scale-max'])

    yticks = samples.index.values
    xticks = [taxid_to_name(o) for o in samples.columns.values]
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

    heatmap.get_figure().savefig('heatmap.png', bbox_inches='tight')
