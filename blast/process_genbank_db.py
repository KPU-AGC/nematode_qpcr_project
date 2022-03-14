"""GenBank Processing Scripts

This script performs two functions: 
1) Creates .fasta file containing sequences from an input GenBank file
2) Generates a metadata file describing each accession of the GenBank file. 

This tool accepts the path to a single GenBank file. 

Requires SeqIO from BioPython ('biopython') and 'pandas'.
"""

import pathlib
import argparse
import pandas as pd
from Bio import SeqIO

def output_fasta(list_fasta, path):
    """Create .fasta files
    Parameters: 
    list_fasta - list - list containing fasta strings
    path - pathLib.path - path without suffix
    Returns: 
    None
    """
    output_file = open(path, 'w')
    for fasta in list_fasta: 
        output_file.write(fasta)
    output_file.close()
def _retrieve_metadata(gb_entry):
    """Parse GenBank entry to obtain relevant metadata information
    Parameters: 
    gb_entry - SeqRecord object
    Returns: 
    metadata - tuple of metadata
    """ 
    def get_source_index(features): 
        """Get the index of the 'source' entry for GenBank entry
        Parameters: 
        features - feature table of a GenBank entry
        Returns: 
        index - int - position of 'source' entry in feature table
        """
        index = None
        for feature in features: 
            if feature.type == "source":
                index = features.index(feature)
        return index
    def get_source_qualifiers(source_qual): 
        """Get qualifiers of interest for source
        Parameters: 
        source_qual - qualifiers in source
        Returns: 
        qual - list - contains qualifiers of interest
        """
        quals = []
        interest = ('isolate', 'host', 'country')
        for qual in interest: 
            if qual in source_qual:
                if len(source_qual[qual]) == 1: 
                    quals.append(source_qual[qual][0])
                else: 
                    quals.append(','.join(source_qual[qual]))
            else: 
                quals.append(None)
        return quals
    #Get qualifiers from 'source' feature
    source_index = get_source_index(gb_entry.features)
    source_qualifiers = get_source_qualifiers(gb_entry.features[source_index].qualifiers)
    #Order: accession id, sequence length, description, organism, isolate, host, country
    metadata = (
        gb_entry.id, 
        len(gb_entry.seq), 
        gb_entry.description, 
        gb_entry.annotations['source'],
        source_qualifiers[0],
        source_qualifiers[1],
        source_qualifiers[2],
    )
    return metadata
def parse_args(): 
    parser = argparse.ArgumentParser(description='Generate .fasta and metadata files from a GenBank file.')
    parser.add_argument(
        'd', 
        metavar='input', 
        type=pathlib.Path, 
        help = 'path to a GenBank file'
    )
    parser.add_argument('-o',
        dest='output_path', 
        action='store', 
        type=pathlib.Path,
        default=None,
        help='output path'
    )
    args = parser.parse_args()
    output_path = args.d.parents
    if args.output_path:
        output_path = args.output_path
    return (args.d, output_path)

def main():

    gb_path, output_path = parse_args()

    if gb_path.is_file() is True:
        #Continue program
        data = SeqIO.parse(gb_path, 'gb')
        #fasta
        list_fasta = []
        #metadata
        metadata = {
            'acc':[],
            'seq_len':[],
            'desc':[],
            'organism':[],
            'isolate':[],
            'host':[],
            'country':[],
        }
        for entry in data: 
            list_fasta.append(entry.format('fasta'))
            meta = _retrieve_metadata(entry)
            #Add to metadata
            metadata['acc'].append(meta[0])
            metadata['seq_len'].append(meta[1])
            metadata['desc'].append(meta[2])
            metadata['organism'].append(meta[3])
            metadata['isolate'].append(meta[4])
            metadata['host'].append(meta[5])
            metadata['country'].append(meta[6])
        
        #DF
        print(len(metadata['acc']))
        print(len(metadata['seq_len']))
        print(len(metadata['desc']))
        print(len(metadata['organism']))
        print(len(metadata['isolate']))
        print(len(metadata['host']))
        print(len(metadata['country']))
        meta_df = pd.DataFrame(data=metadata)

        #Output
        output_fasta(list_fasta, output_path.joinpath(gb_path.stem + '.fasta'))
        meta_df.to_csv(output_path.joinpath(gb_path.stem + '_metadata.csv'), index=False)
    else: 
        print("Not a GenBank file.") 
        
if __name__ == "__main__": 
    main()