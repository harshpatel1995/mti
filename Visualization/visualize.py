"""
MTI Visualization

Usage:
	visualize.py tree <files>...
	visualize.py heatmap <files>...
	visualize.py violin <files>...
        visualize.py scatter meta <metadata> y <var> samples <files>...

Examples:
	visualize.py tree -s s1.gra
"""


from docopt import docopt
from inspect import getmembers, isclass
import commands


# TODO this kind of defeats the purpose of lazy-loading
COMMANDS = [
    'tree',
    'heatmap',
    'violin',
    'scatter'
]


if __name__ == '__main__':
    options = docopt(__doc__, version='mti-vis 1.0')
    print(options)
    command = [key for key, val in options.items() if key in COMMANDS and val is True][0]
    module = getattr(commands, command)
    module_class = getattr(module, command.capitalize())
    module_class(options).run()
