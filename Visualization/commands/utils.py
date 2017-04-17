"""
Utility functions for visualization, typically anything shared among modules
"""

import csv
import itertools
import pandas as pd


def parse_sample_list(sample_list, metadata_filename=None):
    # Get metadata if we have it
    if metadata_filename is not None:
        metadata = parse_metadata(metadata_filename)
    # Get our groups and an accumulator DataFrame to return
    sample_groups = pd.DataFrame()
    groups = sample_list.split(',')
    # Get the DataFrame for each group
    for group in groups:
        group_df = pd.DataFrame()
        samples = group.split('+')
        # Get the DataFrame for each individual sample file
        for sample in samples:
            sample_df = pd.read_csv(sample, sep='\t', header=None,
                    names=['name', 'lineage', 'rel_abund'])
            # TODO: Add metadata here
            group_df = group_df.append(sample_df, ignore_index=True)
        # Average the relative abundance for each organism in the group
        group_df = group_df.groupby(['name', 'lineage']).agg({
                'rel_abund': lambda x: sum(x) / len(samples) }).reset_index()
        # Attach the name of this sample group and append it to the rest of the
        sample_concat_name = ','.join([s.replace(' ', '').split('/')[-1] for s in samples])
        group_df['sample'] = sample_concat_name
        sample_groups = sample_groups.append(group_df, ignore_index=True)
    return sample_groups


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

def pivot_on_sample_and_name(samples):
    """
    TODO
    """
    return samples.pivot(index='sample', columns='name', values='rel_abund')

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
