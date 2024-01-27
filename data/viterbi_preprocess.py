from pathlib import Path
from itertools import combinations, permutations
from tqdm import tqdm

input_filepath = Path('data/test/HIV1_ALL_2021_genome_DNA_reduced.fasta')
output_filepath = Path('data/test/HIV1_ALL_2021_genome_DNA_reduced_preprocessed.fasta')

with open(input_filepath, "r") as input_file:
    sequences = input_file.read().split('>')
    

sequences = sequences[1:]
characters = {'A', 'C', 'G', 'T', '-'}

for i in range(len(sequences)):
    first_newline_idx = sequences[i].find('\n')
    sequences[i] = sequences[i][first_newline_idx + 1:]
    sequences[i] = sequences[i].replace('\n', '')
    sequences[i] = sequences[i].replace('-', '')

def parse_sequences(x, y):

    x, y = list(x), list(y)
    for i in range(len(x)):
        if x[i] not in characters:
            x[i] = ''

    for i in range(len(y)):
        if y[i] not in characters:
            y[i] = ''
    return ''.join(x), ''.join(y)

pairs = [(sequences[i], sequences[i+1]) for i in range(0, len(sequences)-1, 2)]

# print(f'NUM SEQS: {len(sequences)}')
# print(f'NUM PAIRS: {len(list(pairs))}')

with open(output_filepath, 'w') as output_file:
    for sequence_x, sequence_y in tqdm(list(pairs)):
        parsed_x, parsed_y = parse_sequences(sequence_x, sequence_y)
        output_file.write(f'{parsed_x} {parsed_y}')
        output_file.write('\n')