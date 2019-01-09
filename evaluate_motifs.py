#! /usr/bin/env python3
'''
A script to apply the given motifs to the specified data set, calculating
the accuracy, sensitivity, specificity and Matthews' correlation coefficient.

'''
import csv
import re

import click

from utilities import extract_seqs


def get_matches(motif, seqs):
    '''
    Apply the given motif to the input sequences, returning the successful
    matches.

    Args:
        motif (str): The motif pattern to apply.
        seqs (list): The list of sequences.
        
    Returns:
        list: The sequences matching the input motif.

    '''
    regex = re.compile(motif)
    return [s for s in seqs if regex.match(s)]

        
def evaluate_motifs(motif_file, pos_seq_file, neg_seq_file, central_res,
                    seq_len, output_file):
    '''
    Evaluate the given motifs by matching against the positive and negative
    sequence files (prealigned or fasta format) and computing various
    evaluation metrics.
    
    Args:
        motif_file (str): The path to the file containing motif patterns,
                          in CSV format with a column named 'motif'.
        pos_seq_file (str): The path to the file containing positive
                            sequences, i.e. those which are known to bear
                            the modification.
        neg_seq_file (str): The path to the file containing negative
                            sequences, i.e. those which are known not to bear
                            the modification.
        central_res (str): The residue on which sequences should be centered.
        seq_len (int): The length/target length of the sequences.
        output_file (str): The path to the file to which to write the metric
                           results.
    
    '''
    with open(motif_file) as fh:
        reader = csv.DictReader(fh)
        motifs = [l['motif'] for l in reader]
        
    def _extract_seqs(ifile):
        return extract_seqs(ifile, central_res, seq_len, write=False)
    
    pos_seqs = _extract_seqs(pos_seq_file)
    neg_seqs = _extract_seqs(neg_seq_file)
    
    motif_scores = []
    for motif in motifs:
        pos_matches = get_matches(motif, pos_seqs)
        neg_matches = get_matches(motif, neg_seqs)
        tp = len(pos_matches)
        fp = len(neg_matches)
        tn = len(neg_seqs) - fp
        fn = len(pos_seqs) - tp
        
        acc = (tp + tn) / (tp + tn + fp + fn)
        sen = tp / (tp + fn)
        spec = tn / (fp + tn)
        
        mcc = ((tp * tn) - (fp * fn)) / (((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)) ** 0.5)
        motif_scores.append((motif, acc, sen, spec, mcc))
        
    with open(output_file, 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(['motif', 'accuracy', 'sensitivity', 'specificity', 'mcc'])
        writer.writerows(motif_scores)
    

@click.command()
@click.argument('motif-file')
@click.argument('pos-sequence-file')
@click.argument('neg-sequence-file')
@click.argument('central-residue')
@click.argument('seq-length', type=int)
@click.option('--output-file', '-o',
              help='The file to which to save the calculations',
              default='motif_analysis_results.csv')
def main(motif_file, pos_sequence_file, neg_sequence_file, central_residue,
         seq_length, output_file):
    evaluate_motifs(motif_file, pos_sequence_file, neg_sequence_file,
                    central_residue, seq_length, output_file)


if __name__ == '__main__':
    main()