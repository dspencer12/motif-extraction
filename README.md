# motif-extraction

## Introduction

motif-extraction is an all-in-one tool to align peptide sequences, identify overrepresented motifs in a set of peptide sequences and generate sequence logos for each of these motifs. pymotifx.py is the central component of this package and is a python implementation of the [motif-x algorithm](http://motif-x.med.harvard.edu/), based on [rmotifx](https://www.ncbi.nlm.nih.gov/pubmed/26572964).

## Using motif-extraction

To use motif-extraction, you will need to have a setup of input peptide sequences in either fasta, prealigned or unaligned (plain text) format. You will also need a background set of peptide sequences in one of these formats. The lengths of the sequences need to be identical. For fasta and unaligned formats, the tool will automatically ensure that this is the case. For prealigned formats, it is assumed that the tool should make no modifications to the input data.

Basic usage of the motif extraction tool is:
```
python3 extract_motifs.py [INPUT_DATA] [BACKGROUND_DATA] [CENTRAL_RESIDUE] [MIN_OCCURRENCES] [CUTOFF_PVAL]
```
Detailed information on each of the options can be found in the help page for ```extract_motifs.py```, accessed by running
```
python3 extract_motifs.py --help
```

## Package Requirements

motif-extraction depends on a number of additional, but common, python packages. The non-Standard Library dependencies are as follows:
- numpy
- pandas
- scipy
- weblogo

Each of these libraries can be easily installed using the Python package manager, ```pip```.


## References
- Chou MF & Schwartz D (2011). Biological sequence motif discovery using motif-x. Curr Protoc Bioinformatics. Chapter 13:Unit 13.15-24. doi:10.1002/0471250953.bi1315s35.

- Wagih O, Sugiyama N, Ishihama Y, Beltrao P. (2015) Uncovering phosphorylation-based specificities through functional interaction networks (2015). Mol. Cell. Proteomics
