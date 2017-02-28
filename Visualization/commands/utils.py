"""
Utility functions for visualization, typically anything shared among modules
"""

import csv
import itertools
import pandas as pd
from Bio import Entrez

Entrez.email = 'ballardt@knights.ucf.edu'


def parse_gra(filenames, allow_grouping=True, delimiter='\t'):
    """
    Parse a gra file containing taxids, relative abundances, and errors

    Args:
        filenames (string): The name(s) of the .gra file to parse.
        delimiter (string): The delimiter used in the .gra file. Default is tab.

    Returns:
        A dictionary of the form:
        {
            <<taxid (string)>>: {
                rel_abund: <<rel_abund (float)>>, 
                error: <<error (float)>>
            },
            ...
        }
    """
    gra_dict = {}
    if '+' in filenames:
        if allow_grouping:
            means = {}
            filename_list = filenames.split('+')
            for filename in filename_list:
                with open(filename, 'r') as f:
                    reader = csv.reader(f, delimiter=delimiter)
                    l = list(reader)
                    taxids = l[0]
                    rel_abunds = map(float, l[1])
                    errors = map(float, l[2])
                    data = [{'rel_abund': r, 'error': e} for r, e in zip(rel_abunds, errors)]
                    sample = dict(zip(taxids, data))
                    for taxid, val in sample.items():
                        t = means.get(taxid, {'rel_abund': 0, 'error': 0, 'count': 0})
                        t['rel_abund'] += val['rel_abund']
                        t['error'] += val['error']
                        t['count'] += 1
                        means[taxid] = t
            
            for taxid, val in means.items():
                val['rel_abund'] /= val['count']
                val['error'] /= val['count']
                
            gra_dict = means
        else:
            print('Error: Cannot group samples for this kind of visualization')
            exit(1)
    else:
        with open(filenames, 'r') as f:
            reader = csv.reader(f, delimiter=delimiter)
            l = list(reader)
            taxids = l[0]
            rel_abunds = map(float, l[1])
            errors = map(float, l[2])
            data = [{'rel_abund': r, 'error': e} for r, e in zip(rel_abunds, errors)]
            gra_dict = dict(zip(taxids, data))

    return gra_dict


def parse_metadata(filename):
    """
    Parse a metadata CSV file

    Args:
        filename (string): The name of the metadata file.

    Returns:
        A dictionary of the form:
        {
            <<sample name (string)>>: {
                <<column name (string)>>: <<column value (string)>> 
            },
            ...
        }
    """
    with open(filename, 'r') as f:
        reader = csv.DictReader(f)
        metadata = {}
        for row in reader:
            key = row.pop('sample_filename')
            if key in metadata:
                print('WARNING: duplicate sample_id "' + key + '" found in "' + filename + '". The previous value(s) will be overwritten.')
            metadata[key] = row
        return metadata

def name_to_taxid(name):
    """
    Convert a species name to a taxid

    Args:
        name (string): The name of the organism in NCBI's taxonomy database

    Returns:
        The taxid of the organism (string)
    """
    name = name.replace(' ', '+').replace('[','').replace(']','').strip()
    search = Entrez.esearch(term=name, db='taxonomy', retmode='xml')
    return Entrez.read(search)['IdList'][0]


def taxid_to_name(taxid):
    """
    TODO
    """
    handle = Entrez.efetch(id=[taxid], db='taxonomy', mode='text', rettype='xml')
    [taxon] = Entrez.read(handle)
    return taxon['ScientificName']

def parse_sample_group_string(sample_group_string, allow_grouping=True):
    """
    Given a sample group string, return a corresponding pandas dataframe
    """
    gra_filenames = sample_group_string.split(',')

    samples = {g: parse_gra(g, allow_grouping=allow_grouping) for g in gra_filenames}
    data = [[
        {
            'sample': sample,
            'organism': taxid,
            'rel_abund': vals['rel_abund'],
            'error': vals['error']
            }
        for taxid, vals in d.items()
        ]
        for sample, d in samples.items()
        ]
    data = list(itertools.chain.from_iterable(data))
    return pd.DataFrame(data=data)
