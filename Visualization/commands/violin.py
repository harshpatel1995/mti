"""
Violin plot of organisms in each sample
"""


from .base import Base
from .utils import parse_gra, taxid_to_name
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
import itertools


class Violin(Base):
	def run(self):
		# Get the data
		gra_filenames = self.options['<files>']
		samples = {g: parse_gra(g) for g in gra_filenames}
		l = [[{
						'sample': sample,
						'organism': taxid, 
						'rel_abund': vals['rel_abund'],
						'error': vals['error']
					}
				for taxid, vals in d.items()] for sample, d in samples.items()]
		data = list(itertools.chain.from_iterable(l))
		samples = pd.DataFrame(data=data)

		# Draw the plot
		yticks = [taxid_to_name(o) for o in samples['organism']]
		violinPlot = sb.violinplot(
			x='rel_abund', 
			y='organism', 
			data=samples,
			orient='h'
		)

		# Style
		violinPlot.set_yticklabels(yticks)

		violinPlot.get_figure().savefig('violin.png', bbox_inches='tight')