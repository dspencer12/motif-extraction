#! /usr/bin/env python3

import argparse
import os
import re

import pandas as pd
import weblogolib as wl


LOGO_DIR = 'logos'


def generate_logo(seqfile):
    '''
    Generate the sequence logo from the specified sequences.

    Args:
        seqfile (str): The path to a sequence file.

    '''
    with open(seqfile, 'r') as fh:
        seqs = wl.read_seq_data(fh)

    data = wl.LogoData.from_seqs(seqs)

    options = wl.LogoOptions()
    options.title = 'Title 1'

    form = wl.LogoFormat(data, options)

    eps = wl.eps_formatter(data, form)

    with open(seqfile + '.eps', 'wb') as fh:
        fh.write(eps)


def generate_logos(motifs, seqs):
    '''
    Given the extracted motifs, generate sequence logos.

    Args:
        motifs (pandas.DataFrame)
        seqs (list): The list of foreground sequences.

    '''
    if not os.path.isdir(LOGO_DIR):
        os.mkdir(LOGO_DIR)

    for idx, row in motifs.iterrows():
        motif = row.loc['motif']
        regex = re.compile(motif)
        matches = [s for s in seqs if regex.match(s)]

        if len(matches) != row.loc['num_fg_matches']:
            print('WARNING: fewer matches found than expected')

        match_file = os.path.join(LOGO_DIR, 'motif_' + motif + '.txt')
        with open(match_file, 'w') as fh:
            fh.write('\n'.join(matches))

        generate_logo(match_file)


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
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    with open(args.seq_file, 'r') as fh:
        seqs = [s.rstrip() for s in fh]
    generate_logos(pd.read_csv(args.motif_csv), seqs)
