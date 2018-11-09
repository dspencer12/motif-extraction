#! /usr/bin/env python3

import argparse
import os
import re

import pandas as pd
import weblogolib as wl


def generate_logo(seqfile, title):
    '''
    Generate the sequence logo from the specified sequences.

    Args:
        seqfile (str): The path to a sequence file.

    '''
    with open(seqfile, 'r') as fh:
        seqlen = len(fh.readline().rstrip('\n'))
        fh.seek(0)
        seqs = wl.read_seq_data(fh)

    data = wl.LogoData.from_seqs(seqs)

    options = wl.LogoOptions()
    options.title = title
    options.fineprint = ''
    #options.stack_width = 16

    options.first_index = -1 * int(seqlen / 2)

    form = wl.LogoFormat(data, options)

    eps = wl.eps_formatter(data, form)
    eps_file = seqfile[:-4] + '.eps'

    with open(eps_file, 'wb') as fh:
        fh.write(eps)


def generate_logos(motifs, seqs, output_dir):
    '''
    Given the extracted motifs, generate sequence logos.

    Args:
        motifs (pandas.DataFrame)
        seqs (list): The list of foreground sequences.
        output_dir (str): The directory to which to save logos.

    '''
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    for idx, row in motifs.iterrows():
        motif = row.loc['motif']
        regex = re.compile(motif)
        matches = [s for s in seqs if regex.match(s)]

        match_file = os.path.join(output_dir, 'motif_' + motif + '_motif.txt')
        with open(match_file, 'w') as fh:
            fh.write('\n'.join(matches))

        generate_logo(match_file, motif)


def parse_args():
    '''
    Parse the command line arguments.

    Returns:
        argparse.Namespace

    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'motif_csv',
        help='The CSV file containing motif information'
    )
    parser.add_argument(
        'seq_file',
        help='The path to the file containing the foreground sequences'
    )
    parser.add_argument(
        '-o',
        '--output-dir',
        dest='output_dir',
        default='logos',
        help='The destination directory for saving logos'
    )
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    with open(args.seq_file, 'r') as fh:
        seqs = [s.rstrip() for s in fh]
    generate_logos(pd.read_csv(args.motif_csv), seqs, args.output_dir)
