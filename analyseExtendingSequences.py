#!/usr/bin/env python3

import os
import csv

from argparse import ArgumentParser
from collections import defaultdict

from Bio import SeqIO

import matplotlib.pyplot as plt
import numpy as np

def extract_longest_run(sequence):
    '''
    Extract the longest run of consecuative extending sequneces e.g.

    If you have the sequences:

    ACTATATAGCGCGC
    ACTATATAGCGCG
    ACTATATA
    ACTATAT
    ACTATA
    ACT

    The input array to the function would be [14, 13, 8, 7, 6, 3]
    And the output would be [8, 7, 6] as that is the longest run of
    extending sequences
    '''
    last = sequence[0]

    longest_run = []
    current_run = []
    for item in sequence[1:]:
        current_run.append(last)
        if item == last - 1:
            last = item
        else:
            if len(current_run) > 1 and len(current_run) > len(longest_run):
                longest_run = current_run
            
            current_run = []
            last = item

    current_run.append(last)
    if len(current_run) > 1 and len(current_run) > len(longest_run):
        longest_run = current_run

    return longest_run

def to_int_array(string_array):
    '''
    Cast all items in an array from strings to integers
    '''
    result = []
    for i in string_array:
        result.append(int(i))

    return result

def to_str_array(int_array):
    '''
    Cast all items in an array from integers to strings
    '''
    result = []
    for i in int_array:
        result.append(str(i))

    return result

def extract_extending_sequences(seqs):
    '''
    Rename sequences to show the longest run of extending sequences.
    Yields the sequence and an array of subsequence lengths that are
    small RNAs
    '''
    for seq in seqs:
        # extract sequence lengths from the FASTA header
        lengths = to_int_array(seq.id.split('|')[1:])

        longest_run = extract_longest_run(lengths)
        if len(longest_run) > 1:
            new_id = seq.id.split('|')[0] + '|' + '|'.join(to_str_array(longest_run))
            seq = seq[:longest_run[0]]
            seq.id = new_id

            yield seq, longest_run

if __name__ == '__main__':
    parser = ArgumentParser(description='Build graphs from a FASTQ file of collapsed sequences')

    parser.add_argument('collapsed_file', help='Path to the collapsed sequence FASTQ')

    args = parser.parse_args()

    length_map = {x: 0 for x in range(18, 51)}
    result_map = {x: 0 for x in range(18, 51)}
    result_seq_map = defaultdict(lambda: [])

    # Count how common each number of consecutive extending sequences are
    frequencys = []
    for sequence, lengths in extract_extending_sequences(SeqIO.parse(args.collapsed_file, 'fastq')):
        for l in lengths[1:]:
            length_map[l] += 1

        result_map[len(sequence)] +=  1
        result_seq_map[len(lengths)].append(sequence)
        frequencys.append(len(lengths))

    
    # Draw plots based on this
    xs = list(range(18, 51))
    ys = []
    for i in xs:
        ys.append(length_map[i])

    plt.bar(xs, ys)
    plt.xlabel('Length of sequences')
    plt.ylabel('Number of sequences')
    plt.title('Length of sequences that collapse into longer sequences in the dataset')

    plt.savefig('lengthsCollpasedInto.png')
    plt.close()

    with open('lengthsCollpasedInto.csv', 'w') as f:
        writer = csv.writer(f)
        writer.writerow(['Length of sequences'] + xs)
        writer.writerow(['Number of sequence'] + ys)

    xs = list(range(18, 51))
    ys = []
    for i in xs:
        ys.append(result_map[i])

    plt.bar(xs, ys)
    plt.xlabel('Length of sequences')
    plt.ylabel('Number of sequences')
    plt.title('Length of the longest sequence that has at least two sequences collapsed into it')

    plt.savefig('lengthsOfLongest.png')
    plt.close()

    with open('lengthsOfLongest.csv', 'w') as f:
        writer = csv.writer(f)
        writer.writerow(['Length of sequences'] + xs)
        writer.writerow(['Number of sequence'] + ys)

    plt.hist(frequencys, bins=np.arange(min(frequencys), max(frequencys)+1))
    plt.xlabel('Number of sequences')
    plt.ylabel('Frequency')
    plt.title('Number of sequences that extending sequences are collapsed from')

    plt.savefig('numberCollapedInto.png') 

    with open('numberCollaspedInto.csv', 'w') as f:
        writer = csv.writer(f)
        writer.writerow(frequencys)

    os.mkdir('Sequneces')
    for seqFreq in result_seq_map.keys():
        SeqIO.write(result_seq_map[seqFreq], f'Sequneces/frequency{seqFreq:03}.fastq', 'fastq')