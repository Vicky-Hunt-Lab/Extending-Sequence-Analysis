# Extending-Sequence-Analysis
Scripts for the analysis of extending sequences of small RNA

## Usage

These scripts use the python libraries `anytree`, `biopython`, `numpy` and `matplotlib`. These can be installed using pip:

```
$ pip3 install anytree biopython numpy matplotlib
```

To use you need a file of small RNA sequences in either FASTA or FASTQ format. Firstly, to collapse and analyse them run them though the `trieCollapse.py` script as follows:

```
$ python3 trieCollapse.py <input-file>
```

This will produce the file `collapsedSequences.fasta` which contains information about which sequences are shorter versions of others in the header of the file. To extract the longest chain from each one and produce various plot with the information about these sequences, run the script `analyseExtendingSequences.py` on this with the command:

```
$ python3 analyseExtendingSequences.py collapsedSequences.fasta
```

*Note `trieCollapse.py` has a few other filtering options not described here. Use `python3 trieCollapse.py --help` to see them.*
