from pathlib import Path
from itertools import combinations, permutations
from tqdm import tqdm

input_filepath = Path('data/train/HIV1_ALL_2021_genome_DNA_reduced.fasta')
output_filepath = Path('data/train/HIV1_ALL_2021_genome_DNA_reduced_preprocessed.fasta')

with open(input_filepath, "r") as input_file:
    sequences = input_file.read().split('>')

sequences = sequences[1:]

for i in range(len(sequences)):
    first_newline_idx = sequences[i].find('\n')
    sequences[i] = sequences[i][first_newline_idx + 1:]
    sequences[i] = sequences[i].replace('\n', '')


def parse_sequences(x, y):
    assert len(x) == len(y)

    ix, iy = 0, 0
    x, y = list(x), list(y)
    for i in range(len(x)):
        should_delete = True
        if x[i] == '-':
            ix = i
            if y[i] == '-':
               iy = i
            elif i != 0 and y[i-1] == '-':
                iy = i-1
            elif i != len(x)-1 and y[i+1] == '-':
                iy = i+1
            else:
                should_delete = False
            if should_delete:
                x[ix] = ''
                y[iy] = ''
    return ''.join(x), ''.join(y)


pairs = permutations(sequences, 2)

# print(f'NUM SEQS: {len(sequences)}')
# print(f'NUM PAIRS: {len(list(pairs))}')

with open(output_filepath, 'w') as output_file:
    for sequence_x, sequence_y in tqdm(list(pairs)):
        parsed_x, parsed_y = parse_sequences(sequence_x, sequence_y)
        output_file.write(f'{parsed_x[:150]} {parsed_y[:150]}')
        output_file.write('\n')