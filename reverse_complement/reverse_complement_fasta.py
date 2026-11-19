#!/usr/bin/env python3
"""
Reverse complement antisense sequences in a FASTA file.
Sequences with '::-' in header are antisense and get reverse complemented.
Sequences with '::+' in header are sense and kept as-is.
"""

import sys
import os


def reverse_complement(seq: str) -> str:
    """Return reverse complement of DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
                  'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
                  'N': 'N', 'n': 'n'}
    return ''.join(complement.get(base, base) for base in reversed(seq))


def process_fasta(input_file: str, output_file: str):
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        header = None
        sequence_lines = []
        
        for line in f_in:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if header is not None:
                    sequence = ''.join(sequence_lines)
                    if '::-' in header:
                        sequence = reverse_complement(sequence)
                        header = header.replace('::-', '::+', 1)
                    f_out.write(f'{header}\n{sequence}\n')
                header = line
                sequence_lines = []
            else:
                sequence_lines.append(line)
        
        if header is not None:
            sequence = ''.join(sequence_lines)
            if '::-' in header:
                sequence = reverse_complement(sequence)
                header = header.replace('::-', '::+', 1)
            f_out.write(f'{header}\n{sequence}\n')


if __name__ == '__main__':
    input_file = sys.argv[1]

    base, ext = os.path.splitext(input_file)
    output_file = f'{base}_rc{ext}'
    
    process_fasta(input_file, output_file)
    print(f'Processed {input_file} -> {output_file}')