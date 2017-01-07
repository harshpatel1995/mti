"""
A heatmap of taxons to samples
"""


from .base import Base
from .utils import parse_gra, taxid_to_name
import pandas as pd
import seaborn as sb
import itertools


class Heatmap(Base):
	def run(self):
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
		df = pd.DataFrame(data=data)
		samples = df.pivot ('sample', 'organism', 'rel_abund')
		sb.heatmap(samples).get_figure().savefig('heatmap.png')