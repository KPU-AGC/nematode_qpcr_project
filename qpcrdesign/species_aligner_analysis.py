from Bio import SeqIO, AlignIO, Seq
import argparse
import pathlib
import csv
import subprocess

class genbankHandler(): 
    def __init__(self, gb_path, output_path): 
        def parse_gb(gb_path): 
            gb_data = []
            data = SeqIO.parse(gb_path, 'gb')
            for item in data: 
                gb_data.append(item)
            return gb_data
        #data - contains raw gb file
        self.data = parse_gb(gb_path)
        #id - comes from gb file name
        self.id = gb_path.stem
        self.path = output_path
        self.fasta_paths = []
        self.alignment_paths = []
        self.consensus_path = None
        #gb - dictionary of species names mapped to a list of data corresponding to each sequence belonging to that species
        self.gb = dict()
        self.species_keys = self.gb.keys()
        self.consensus_sequences = None
        self._split_gb_by_species()

    def _split_gb_by_species(self):
        for gb_data in self.data: 
            organism = gb_data.annotations['organism']
            species_name = organism.split(' ')[1]
            #If there's no clear species designation, ignore
            if (species_name != 'sp.') and (species_name != 'cf.') and (species_name !='aff.'):
                if organism not in self.gb.keys(): 
                    self.gb[organism] = [gb_data]
                else: 
                    self.gb[organism].append(gb_data)

    def generate_metadata(self): 
        #We want to know the following: 
        #1) Number of species present
        #2) Number of species that did not have a species call
        #3) Number of sequences belonging to each species
        no_species_call = []
        #Identifying which accessions do not have a species call...
        for gb_data in self.data: 
            organism = gb_data.annotations['organism']
            accession= gb_data.id
            species_name = organism.split(' ')[1]
            if (species_name == 'sp.') or (species_name == 'cf.') or (species_name =='aff.'): 
                no_species_call.append((organism, accession))
        #Get number of species
        num_species = len(self.gb.keys())

        #Open the file
        metadata_file_path = self.path.joinpath(f'{self.id}_metadata.csv')
        metadata_file = open(metadata_file_path, 'w')
        metadata_file.write(f'Species ({str(num_species)}):\n')
        #Write the species and the number of sequences available for each one
        for species in self.gb: 
            metadata_file.write(f'{species},{len(self.gb[species])}\n')
        #write the no species
        metadata_file.write(f'No species call ({str(len(no_species_call))}):\n')
        for sequence in no_species_call: 
            metadata_file.write(f'{sequence[0]},{sequence[1]}\n')
        metadata_file.close()
        
    def output_species_gb(self): 
        for organism in self.gb.keys(): 
            species_path = self.path.joinpath(f'{organism.replace(" ","-")}_data.gb')
            SeqIO.write(self.gb[organism], species_path, 'gb')

    def output_species_fasta(self): 
        for organism in self.gb.keys(): 
            species_path = self.path.joinpath(f'{organism.replace(" ","-")}.fasta')
            SeqIO.write(self.gb[organism], species_path, 'fasta')
            self.fasta_paths.append(pathlib.Path(species_path))

    def output_species_metadata(self): 
        #Loop through each organism
        for organism in self.gb.keys(): 
            species_path = self.path.joinpath(f'{organism.replace(" ","-")}_metadata.csv')
            all_gb_data = []
            #Loop through each Genbank file
            for gb in self.gb[organism]: 
                gb_metadata = [
                    gb.id, 
                    len(gb.seq), 
                    gb.description,
                    gb.annotations['source'], 
                ]
                #Identify source feature, then store data
                for feature in gb.features: 
                    if feature.type == 'source': 
                        if 'isolate' in feature.qualifiers: 
                            gb_metadata.append(feature.qualifiers['isolate'])
                        else: 
                            gb_metadata.append('')
                        if 'host' in feature.qualifiers: 
                            gb_metadata.append(feature.qualifiers['host'])
                        else: 
                            gb_metadata.append('')
                        if 'country' in feature.qualifiers: 
                            gb_metadata.append(feature.qualifiers['country'])
                        break
                all_gb_data.append(gb_metadata)
            metadata_output = open(species_path, 'w', newline='')
            csv_writer = csv.writer(metadata_output)
            csv_writer.writerow(
                (
                    'accession_id',
                    'seq_len',
                    'description',
                    'organism',
                    'isolate',
                    'host',
                    'country'
                )
            )
            csv_writer.writerows(all_gb_data)
            metadata_output.close()

    def generate_species_alignment(self): 
        for fasta_path in self.fasta_paths: 
            aligned_path = fasta_path.with_name(fasta_path.stem + '_aligned.fasta')
            args =[
                'mafft',
                '--auto',
                str(fasta_path), 
            ]
            result = subprocess.run(args, capture_output=True)
            decoded = result.stdout.decode('utf-8')
            output_file = open(aligned_path, 'w')
            output_file.write(decoded)
            output_file.close()
            self.alignment_paths.append(pathlib.Path(aligned_path))

    def generate_consensus(self): 
        def get_seq_position_info(alignment): 
            """
            Get information for each position of the alignment.
            Output is in list format, with each index corresponding to the 0-based position
            of the alignment. 
            Each index contains dictionary with following information: 
            pos - position number
            accessions - the accessions represented at that position
            bases - base calls in corresponding order of accessions
            seqs_rep - number of accessions represented at that position
            p_seq_rep - percentage of accessions represented at that position
            """
            def get_sequence_regions(alignment): 
                """
                Identify the regions that each accession spans on the alignment.
                Algorithm: 
                1) From beginning of sequence, keep going until position is not '-',
                then record the position as the start of the region
                2) From the end of the sequence, go in reverse until position is not '-',
                then record the position as the end of the region
                3) Return tuple --> (accession, start, end)
                """
                def get_sequence_region(accession): 
                    #Variables
                    f_start = False
                    f_end = False
                    start = 0
                    end = 0
                    #Data from SeqRecord object 
                    sequence = accession.seq
                    acc = (accession.id).split('.')[0]
                    #From beginning of the sequence, keep going until position is not '-'
                    i = 0
                    while (f_start is False): 
                        if sequence[i] != '-': 
                            start = i
                            f_start = True
                        i = i + 1
                    #From ending of the sequence, keep going until position is not '-'
                    i = -1
                    while (f_end is False): 
                        if sequence[i] != '-':
                            end = i%len(sequence) #Translate negative index to positive index
                            f_end = True
                        i = i - 1
                    return ((acc, start, end))
                sequence_regions = []
                for accession in alignment: 
                    sequence_regions.append(get_sequence_region(accession))
                return sequence_regions
            def get_contributing_seqs(position, list_seq_region):
                """
                Given a position and the list of sequence regions for each accession,
                identify which accessions are represented
                """
                list_accession = []
                list_indices = []
                for sequence in list_seq_region:
                    if (position >= sequence[1]) and (position <= sequence[2]): 
                        list_accession.append(sequence[0])
                        list_indices.append(list_seq_region.index(sequence))
                return list_accession, list_indices
            list_seq_position_info = []
            #Get sequence regions
            list_seq_regions = get_sequence_regions(alignment)
            num_accessions = len(alignment)
            for position in range(alignment.get_alignment_length()):
                list_accessions, list_indices = get_contributing_seqs(position, list_seq_regions)
                list_bases = []
                #Append the basecalls of every contributing sequence at specified position
                for index in list_indices: 
                    list_bases.append(alignment[index].seq[position])
                position_info = {
                    'pos':position,
                    'accessions':tuple(list_accessions),
                    'bases':tuple(list_bases),
                    'seqs_rep':len(list_accessions),
                    'p_seq_rep':len(list_accessions)/num_accessions,
                }
                list_seq_position_info.append(position_info)
            return list_seq_position_info
        def get_consensus_seq(seq_info): 
            """
            Determine the consensus sequence of an alignment, and create position matrix
            Definition of consensus: most common base represented at that position. 
            """
            consensus_sequence = []
            for position in seq_info: 
                #Ignore any ambiguous basecalls - accept A, T, C, G, and 'gap'
                base_counts = {
                    'a':position['bases'].count('a')+position['bases'].count('A'),
                    't':position['bases'].count('t')+position['bases'].count('T'),
                    'c':position['bases'].count('c')+position['bases'].count('C'),
                    'g':position['bases'].count('g')+position['bases'].count('G'),
                    '-':position['bases'].count('-'),
                }
                #print(base_counts)
                max_basecalls = [key for key, count in base_counts.items() if count == max(base_counts.values())]
                if len(max_basecalls) == 1: 
                    consensus_sequence.append(max_basecalls[0])
                else: 
                    consensus_sequence.append('n')
            return str(Seq.Seq(''.join(consensus_sequence)).ungap(gap='-'))
        def create_fasta(cons_seqs, output_path): 
            output_file = open(output_path, 'w')
            for key in cons_seqs: 
                output_file.write(f'>{key}\n')
                output_file.write(cons_seqs[key])
                output_file.write('\n')
            output_file.close()
        consensus_seqs = dict()
        for alignment_path in self.alignment_paths: 
            alignment = AlignIO.read(alignment_path, 'fasta')
            alignment_id = alignment_path.name.split('_')[0]
            seq_info = get_seq_position_info(alignment)
            consensus_seq = get_consensus_seq(seq_info)
            consensus_seqs[alignment_id] = consensus_seq
        self.consensus_sequences = consensus_seqs
        #Output
        self.consensus_path = self.path.joinpath(self.id + '_consensus.fasta')
        create_fasta(self.consensus_sequences, self.consensus_path)
    
    def generate_consensus_alignment(self): 
        consensus_alignment_path = self.path.joinpath(self.id + '_consensus_aligned.fasta')
        args =[
                'mafft',
                '--auto',
                str(self.consensus_path), 
            ]
        result = subprocess.run(args, capture_output=True)
        decoded = result.stdout.decode('utf-8')
        output_file = open(consensus_alignment_path, 'w')
        output_file.write(decoded)
        output_file.close()

def parse_args(): 
    parser = argparse.ArgumentParser('Analyze and align multiple sequences for genus-level analysis')
    parser.add_argument(
        'genus_gb', 
        action='store',
        type=pathlib.Path,
        help='path to genbank file containing all data'
    )
    parser.add_argument(
        '--output',
        '-o',
        action='store',
        dest='output_path',
        default=None,
        type=pathlib.Path,
        help='Output destination path',
    )
    args = parser.parse_args()
    if args.output_path is None: 
        args.output_path = args.genus_gb.parent
    return (args.genus_gb, args.output_path)

def main(): 
    """
    Algorithm: 
    1) Open genus genbank file --> this is generated using Entrez query on NCBI Nucleotide
    2) Split gb into several species genbank files
    3) Align .fasta files for each species

    """
    genus_gb_path, output_path = parse_args()
    genus_data = genbankHandler(genus_gb_path, output_path)
    genus_data.output_species_gb()
    genus_data.output_species_fasta()
    genus_data.output_species_metadata()
    genus_data.generate_species_alignment()
    genus_data.generate_consensus()
    genus_data.generate_consensus_alignment()
    genus_data.generate_metadata()

if __name__ == '__main__': 
    main()