"""
Utility functions for visualization, typically anything shared among modules
"""

import csv
from Bio import Entrez

Entrez.email = 'ballardt@knights.ucf.edu'


def parse_gra(filename, delimiter='\t'):
  """
  Parse a gra file containing taxids, relative abundances, and errors

  Args:
    filename (string): The name of the .gra file to parse.
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
  with open(filename, 'r') as f:
    reader = csv.reader(f, delimiter=delimiter)
    l = list(reader)
    taxids = l[0]
    rel_abunds = map(float, l[1])
    errors = map(float, l[2])
    data = [{'rel_abund': r, 'error': e} for r, e in zip(rel_abunds, errors)]
    return dict(zip(taxids, data))


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
  name = name.replace(' ', '+').strip()
  search = Entrez.esearch(term=name, db='taxonomy', retmode='xml')
  return Entrez.read(search)['IdList'][0]


def taxid_to_name(taxid):
  """
  TODO
  """
  handle = Entrez.efetch(id=[taxid], db='taxonomy', mode='text', rettype='xml')
  [taxon] = Entrez.read(handle)
  return taxon['ScientificName']
