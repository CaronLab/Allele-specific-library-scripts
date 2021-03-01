#!/bin/python3

from sys import argv

if len(argv) < 3:
    print('Usage: >python3 build_spectrast_incl_lists.py /path/to/annotated_file.tsv percent_rank_cutoff\n'
        'Example: >python3 ./build_spectrast_incl_lists.py ./interact.iproph_annotated.tsv 0.5')
    exit(0)

pep_file = argv[1]
cutoff = float(argv[2])

with open(pep_file, 'r') as f:
    header = f.readline().strip().split()
    contents = [x.strip().split() for x in f.readlines()]

label_index = header.index('Label')
rank_headers = [x for x in header if '_rank' in x]
rank_indices = [header.index(x) for x in rank_headers]
allele_by_index = {x: header[x].replace('_rank', '') for x in rank_indices}

pep_index = header.index('SpectraST_Peptide')

for allele_index in rank_indices:
    with open(f'{allele_by_index[allele_index]}_inclusion_list.tsv', 'w') as f:
        for line in contents:
            if line[label_index].lower() != 'target':
                continue
            if float(line[allele_index]) <= cutoff:
                f.write(line[pep_index] + '\n')
