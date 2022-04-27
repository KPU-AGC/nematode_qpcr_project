import argparse
import pathlib
from Bio import AlignIO

def parse_args(): 
    parser = argparse.ArgumentParser('Program to retrieve accessions used to construct an alignment')
    parser.add_argument('alignment_path', 
        action='store',
        type=pathlib.Path,
        help='Path to alignment'
    )
    parser.add_argument('-o', 
        dest='output_path', 
        type=pathlib.Path, 
        default=None,
        help='Output path'
    )
    args = parser.parse_args()
    output_path = args.alignment_path.parent
    if args.output_path: 
        output_path = args.output_path
    return (args.alignment_path, output_path)

def main(): 
    alignment_path, output_path = parse_args()

    alignment = AlignIO.read(alignment_path, 'fasta')

    list_id = []
    for seq in alignment: 
        list_id.append(seq.id.split('.')[0])

    #Output
    output_file = open(output_path.joinpath('list_accessions.txt'), 'w')
    for id in list_id: 
        output_file.write(f'{id}\n')
    output_file.close()

if __name__ == '__main__':
    main()