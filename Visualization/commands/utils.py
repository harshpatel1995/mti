"""
Utility functions for visualization, typically anything shared among modules
"""

import csv
import itertools
import pandas as pd
import uuid
import math

def unique_filename():
    return str(uuid.uuid1())

def parse_sample_list(sample_list, metadata_filename=None):
    group_list = ['name', 'lineage']
    # Get metadata if we have it
    if metadata_filename is not None:
        meta_df = parse_metadata(metadata_filename)
        group_list.extend(filter(lambda c : c != 'sample', 
                list(meta_df.columns.values)))
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
            if metadata_filename is not None:
                meta_row = meta_df.loc[meta_df['sample'] == clean_filename(sample)]
                if not meta_row.empty:
                    for c in meta_row.columns:
                        if c == 'sample':
                            continue
                        sample_df[c] = list(meta_row[c])[0]
            group_df = group_df.append(sample_df, ignore_index=True)
        # Average the relative abundance for each organism in the group
        group_df = group_df.groupby(group_list).agg({
                'rel_abund': lambda x: sum(x) / len(samples) }).reset_index()
        # Attach the name of this sample group and append it to the rest of the
        sample_concat_name = ','.join([clean_filename(s) for s in samples])
        group_df['sample'] = sample_concat_name
        sample_groups = sample_groups.append(group_df, ignore_index=True)
    return sample_groups

def parse_metadata(filename):
    """
    Parse a metadata CSV file

    Args:
        filename (string): The name of the metadata file.

    Returns:
        A pandas dataframe of the metadata
    """
    with open(filename, 'r') as f:
        meta_df = pd.read_csv(filename)
        reserved_words = ['name', 'lineage', 'rel_abund']
        rename_dict = {w: '_'+w for w in reserved_words if w in meta_df.columns}
        meta_df = meta_df.rename(columns=rename_dict)
        return meta_df

def clean_filename(filename):
    """
    todo: document
    """
    return filename.replace(' ', '').split('/')[-1]

def pivot_on_sample_and_name(samples):
    """
    TODO: document
    """
    return samples.pivot(index='sample', columns='name', values='rel_abund')

def set_common_ancestor(df, ancestor):
    """
    TODO: document
    """
    df = df.loc[df['lineage'].str.contains(r'\s?'+ancestor+'(;|$)', regex=True)]
    df_ra = df['rel_abund'] / df['rel_abund'].sum()
    df['lineage'] = df['lineage'].apply(
            lambda s: s.split(ancestor+';')[0] + ancestor)
    df = df.assign(rel_abund=df_ra)
    return df