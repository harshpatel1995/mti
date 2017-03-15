"""
A taxonomic tree
"""


from .utils import parse_gra, taxid_to_name, name_to_taxid

import sys
import time
from ete3 import PhyloTree, TreeStyle, NodeStyle, faces
from Bio import Entrez

Entrez.email = 'ballardt@knights.ucf.edu'


def run(options):
    gra_filename = options['<sample_group_string>']
    sample_organisms = parse_gra(gra_filename)
    t = construct_taxonomic_tree(sample_organisms, get_lineage_entrez)

    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.branch_vertical_margin = 10
    ts.layout_fn = tree_layout

    t.render('tree.png', tree_style=ts)


def construct_taxonomic_tree(sample_organisms, lineage_query_func):
    """
    TODO update this documentation
    Create a taxonomic tree from a list of taxids

    Args:
            taxids ([string]): The list of taxids that represent leaf nodes
            lineage_query_func (func): A function that produces a list of lineage taxids
                            given a taxid

    Returns:
            The root of the created tree
    """

    # Keeping track of nodes made thus far allows us to create branches on the
    # tree and memoize
    nodes = {}

    # The element that will be returned; the root of the tree
    root = None

    # Go through each organism in the sample, fetching its lineage and creating a
    # node for it
    for taxid, data in sample_organisms.items():
        rel_abund = data['rel_abund']
        lineage = lineage_query_func(taxid)
        child = tax_node(taxid, rel_abund)

        # For each ancestor of this sample organism that we haven't seen, create a
        # node for it and add the previous node as a child. This creates a path from
        # the leaf to the most ancestral node not yet created.
        for ancestor in (l for l in lineage if l not in nodes):
            nodes[ancestor] = tax_node(ancestor)
            nodes[ancestor].add_child(child)
            nodes[ancestor].rel_abund += child.rel_abund
            child = nodes[ancestor]

        # If we have a new root, set it. Otherwise, connect the most recently 
        # created node to the rest of the tree by attaching it to its parent.
        # TODO this doesn't seem like a great way of doing it
        if child.taxid == lineage[-1]:
            root = child
        else:
            if not child.taxid in lineage:
                lineage.insert(0, child.taxid)
            parentIndex = lineage.index(child.taxid) + 1
            parent = nodes[lineage[parentIndex]]
            parent.add_child(child)
            for l in lineage[1:]:
                nodes[l].rel_abund += child.rel_abund

    return root


def tax_node(taxid, rel_abund=0):
    """
    TODO update this documentation
    Create a node on the taxonomic tree

    Args:
            taxid (string): The taxid of the node

    Returns:
            The newly created node
    """
    node = PhyloTree()
    node.name = taxid
    node.taxid = taxid
    node.rel_abund = rel_abund
    return node


def get_lineage_entrez(taxid):
    """
    Fetch a taxid's lineage on-demand via Entrez

    Args:
            taxid (string): taxid of the organism to get the lineage of

    Returns:
            A list with the lineage taxids ([string])
    """
    time.sleep(1) # To avoid Entrez from resetting our connection
    handle = Entrez.efetch(db='taxonomy', id=[taxid], mode='text', rettype='xml')
    [taxon] = Entrez.read(handle)
    lineage = taxon['Lineage'].split('; ')
    lineage.reverse()
    return [name_to_taxid(s) for s in lineage]


def tree_layout(node):
    """
    TODO
    """
    scientificName = taxid_to_name(node.taxid)

    if node.is_leaf():
        nameSize = 14
        nameColor = '#009000'
        relAbundSize = 10
    else:
        nameSize = 10
        nameColor = '#303030'
        relAbundSize = 8

    nameFace = faces.TextFace(scientificName, fsize=nameSize, fgcolor=nameColor)
    nameFace.margin_bottom = 5
    nameFace.margin_right = 10
    nameFace.margin_left = 5
    faces.add_face_to_node(nameFace, node, column=0, position='branch-top')

    relAbundFace = faces.TextFace(node.rel_abund, fsize=relAbundSize, fgcolor='#2148c8')
    relAbundFace.margin_top = 5
    relAbundFace.margin_left = 5
    faces.add_face_to_node(relAbundFace, node, column=0, position='branch-bottom')

    taxidFace = faces.TextFace(node.taxid, fsize=8, fgcolor='#303030')
    taxidFace.margin_top = 5
    taxidFace.margin_left = 5
    faces.add_face_to_node(taxidFace, node, column=0, position='branch-bottom')
