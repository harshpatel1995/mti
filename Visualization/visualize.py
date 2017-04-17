"""
MTI Visualization

Usage:
    visualize.py tree 
        [--filter=<organism>] 
        [--remove-threshold=<threshold>] 
        <sample_group_string>
    visualize.py heatmap 
        [--filter=<organism>]
        [--scale-min=<scale_min>]
        [--scale-max=<scale_max>]
        [--outliers]
        <sample_group_string>
    visualize.py violin 
        [--filter=<organism>]
        <sample_group_string>
    visualize.py scatter 
        [-n|--normalize-data]
        [-r|--rank-correlations]
        [--filter=<organism>]
        (v|vertical) <sample_group_string_v>
        (h|horizontal) 
            ((meta <metadata> var <var>)|(sample <sample_group_string_h>))
    visualize.py bar
        [--filter=<organism>]
        [--meta=<metadata>] 
        [--var=<var>]
        <sample_group_string>

Explanation of sample_group_string:
    A string representing the sample files to consider, possibly with groups.
    Of the form ([group]=[sample_file]+[sample_file]... | [sample_file]),...
    E.g., "g1=s1.gra+s3.gra+s4.gra,s2.gra" will result in two groups:
        (1) g1, which averages s1, s3, and s4 into one logical sample
        (2) s2
    Similarly, "s3.gra,s4.gra,sample_6.gra" will result in three groups:
        (1) s3
        (2) s4
        (3) sample_6
    Groups (e.g. "g1" in the above example) can have 1 or more sample files. 
    sample_group_string must have at least one sample file or group.
"""


from docopt import docopt
import importlib


# TODO this kind of defeats the purpose of lazy-loading
COMMANDS = [
    'tree',
    'heatmap',
    'violin',
    'scatter',
    'bar'
]


if __name__ == '__main__':
    options = docopt(__doc__, version='mti-vis 1.0')
    command = [key for key, val in options.items() if key in COMMANDS and val is True][0]
    i = importlib.import_module('commands.' + command)
    i.run(options)
