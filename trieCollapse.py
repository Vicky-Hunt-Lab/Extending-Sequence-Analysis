#!/usr/bin/env python3

from math import inf
from argparse import ArgumentParser

from anytree import Node
from Bio import SeqIO, Seq, SeqRecord

class Trie(Node):
    '''
    Class for building a simple inefficent trie, inserts each character as a
    separate node. Also acts as an anytree Node.
    '''

    def __init__(self):
        super(Node, self).__init__()
        self.name = ''

    def add_string(self, string):
        '''
        Insert a string into the trie
        '''
        current_node = self
        string = string.upper() + '$'

        for character in string:
            next_node = Trie.child_has_name(current_node, character)

            if next_node is not None:
                current_node = next_node
            else:
                current_node = Node(character, parent=current_node)

    @staticmethod
    def child_has_name(node, name):
        '''
        Check if a particular child of a node has a certain name
        '''
        for child in node.children:
                if child.name == name:
                    return child

        return None
            
def retrive_longest_keys(trie, only_collapsed=False):
    '''
    Scan the trie for the longest keys and yield them
    '''
    found_count = 0

    for leaf in trie.leaves:
        if len(leaf.parent.children) == 1:
            label = f'SmallRNA{found_count:07}'
            reconstrcuted_key = ''
            
            collapsed_count = 0
            for node in leaf.parent.iter_path_reverse():
                reconstrcuted_key = node.name + reconstrcuted_key

                if Trie.child_has_name(node, '$'):
                    collapsed_count += 1
                    label = label + '|' + str(len(node.ancestors))

            if not only_collapsed or collapsed_count > 1:
                yield reconstrcuted_key, label

                found_count += 1

def into_seqrecord(iterable):
    '''
    Create SeqRecord objects from an iterable that returns sequence and name
    '''
    for seq, name in iterable:
        record = SeqRecord.SeqRecord(seq=Seq.Seq(seq), id=name, description='')
        record.letter_annotations["phred_quality"] = [40] * len(record)

        yield record

if __name__ == '__main__':
    parser = ArgumentParser(description='Collapse a set of sequnces into the longest sequences')

    parser.add_argument('input_fasta', help='Fasta or Fastq file to read sequences from')

    parser.add_argument('-o', '--output', help='Fasta file to write to', default='collapsedSequences.fasta')
    parser.add_argument('-q', '--fastq', help='Output as fastq file', action='store_true')
    parser.add_argument('-c', '--only-collapsed', help='Only produce sequences collapsed from two or more', action='store_true')

    parser.add_argument('-m', '--min-length', help='Minimum length of sequence to include', type=int, default=-inf)
    parser.add_argument('-x', '--max-length', help='Maximum length of sequence to include', type=int, default=inf)

    args = parser.parse_args()

    if args.input_fasta.endswith('.fastq') or args.input_fasta.endswith('.fq'):
        filetype = 'fastq'
    else:
        filetype = 'fasta'

    trie = Trie()
    for seq in SeqIO.parse(args.input_fasta, filetype):
        if len(seq) >= args.min_length and len(seq) <= args.max_length:
            trie.add_string(str(seq.seq))

    if args.fastq:
        filetype = 'fastq'
    else:
        filetype = 'fastq'

    SeqIO.write(into_seqrecord(retrive_longest_keys(trie, only_collapsed=args.only_collapsed)), args.output, filetype)