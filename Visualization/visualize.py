"""
MTI Visualization

Usage:
	visualize.py tree ((-s | --samples) <files>...)
	visualize.py heatmap ((-s | --samples) <files>...)
	visualize.py violin ((-s | --samples) <files>...)

Examples:
	visualize.py tree -s s1.gra
"""


from docopt import docopt
from inspect import getmembers, isclass
import commands


if __name__ == '__main__':
	options = docopt(__doc__, version='mti-vis 1.0')
	command = [key for key, val in options.items() if not key.startswith(('-', '<')) and val is True][0]
	module = getattr(commands, command)
	module_class = getattr(module, command.capitalize())
	module_class(options).run()