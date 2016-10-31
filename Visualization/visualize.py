"""
Visualize the inferred taxonomic structure of the metagenome
"""

import csv
import sys
from ete3 import PhyloTree, TreeStyle
from Bio import Entrez

Entrez.email = 'ballardt@knights.ucf.edu'
CELLULAR_ORGANISMS_TAXID = 131567

def parse_gra(filename, delimiter='\t'):
	"""Parse a gra file containing taxids and relative abundances

	args
		filename 	- .gra file to parse
		delimiter	- delimiter used in the .gra file. Default is tab.

	return
		dictionary with taxids as keys, relative abundances as values
	"""
	with open(filename, 'rb') as f:
		reader = csv.reader(f, delimiter=delimiter)
		l = list(reader)[:2] # only want taxids and rel. abund.
		return dict(zip(map(int, l[0]), map(float, l[1])))


def get_lineages(taxids):
	"""Generate the lineages for the supplied taxids.

	args
		taxids - taxids to get lineages for

	return
		dictionary with taxids as keys, list of lineage taxids as values from closest
		ancestor to farthest
	"""
	lineages = {}
	memo = {}
	for taxid in taxids:
		handle = Entrez.efetch(db='taxonomy', id=taxid, mode='text', rettype='xml')
		taxon = Entrez.read(handle)[0] # Can read return a >1 element list here?
		lineage = taxon['Lineage'].split('; ')
		lineage.reverse()
		lineages[taxid] = list(species_to_taxids(lineage, memo))
	return lineages


def species_to_taxids(species_names, memo={}):
	"""Given a list of species names, get their corresponding taxids.

	args
		species_names - list of species names to convert to taxids
		memo - dictionary of species calculated thus far (for DP)

	return
		list of taxids for each species
	"""
	for species in species_names:
		species = species.replace(' ', '+').strip()
		if species not in memo:
			search = Entrez.esearch(term=species, db='taxonomy', retmode='xml')
			memo[species] = int(Entrez.read(search)['IdList'][0])
		yield memo[species]


def construct_phylo_tree(taxids, rel_abunds=None):
	"""Create a phylogenetic with relative abundance info from a list of taxids

	args
		taxids - list of taxids
		rel_abunds - dictionary of relative abundances for taxids

	return
		root of the created tree
	"""
	nodes = {}
	lineages = get_lineages(taxids)
	for taxid, lineage in lineages.iteritems():
		path = [new_node(taxid, nodes)]
		child = path[0]
		for ancestor in lineage:
			path.append(new_node(ancestor, nodes))
			if child not in path[-1].children:
				path[-1].add_child(child)
			child = path[-1]
	return nodes[CELLULAR_ORGANISMS_TAXID]


def new_node(taxid, nodes={}):
	if taxid not in nodes:
		node = nodes.setdefault(taxid, PhyloTree())
		node.name = taxid
		node.taxid = taxid
	return nodes[taxid]


if __name__ == "__main__":
	gra_filename = sys.argv[1]
	rel_abunds = parse_gra(gra_filename)
	t = construct_phylo_tree(rel_abunds.keys())
	t.render('phylo.png', tree_style=TreeStyle())

	"""
	# Load a tree and link it to an alignment.
	t = PhyloTree("(((seqA,seqB),seqC),seqD);")
	t.link_to_alignment(alignment=fasta_txt, alg_format="fasta")

	t.render('phylo.png');
	"""