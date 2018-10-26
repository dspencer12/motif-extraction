#! /usr/bin/env python3
'''
Extract statistically significant motifs from the input peptide sequences.

'''

import argparse
import enum
import sys

from align import align_sequence
from pymotifx import motifx
import sequence_logos


class Format(enum.Enum):
    FASTA = enum.auto()
    PREALIGNED = enum.auto()
    UNALIGNED = enum.auto()


def _detect_format(ifile):
    '''
    Detect the format of the input file.

    Args:
        ifile (str): The path to the input file.

    Returns:
        Format

    '''
    with open(ifile, 'r') as fh:
        first = fh.readline().rstrip()

        if first.startswith('>'):
            return Format.FASTA

        fh.seek(0)

        # Non-FASTA files
        # TODO: optimise detection of prealigned input
        length = len(first)
        res = first[int(length / 2)]
        for line in fh:
            line = line.rstrip()
            if len(line) != length or line[int(len(line) / 2)] != res:
                return Format.UNALIGNED

    return Format.PREALIGNED


def generate_bg(seqs, central_res, length):
    '''
    Generate the background sequences given full input sequences.

    '''
    hlen = int(0.5 * length)

    bg = []
    for seq in seqs:
        indices = [i for i, res in enumerate(seq) if res == central_res]
        inner_seqs = [seq[max(idx - hlen, 0):min(idx + hlen + 1, len(seq))]
                      for idx in indices]
        bg.extend([align_sequence(s, central_res, length) for s in inner_seqs])

    return bg


def extract_seqs(ifile, central_res, length, is_bg=False, write=True):
    '''
    Process the given input file to return the aligned sequences.

    Args:
        ifile (str): The path to the input file.
        central_res (str): The amino acid residue on which to align the
                           sequences.
        length (int): The total length of each amino acid sequence.
        is_bg (bool, optional): A boolean flag indicating whether the input
                                corresponds to background data.
        write (bool, optional): If True, write the aligned sequences to a file.

    Returns:
        list: The aligned sequences.

    '''
    form = _detect_format(ifile)

    needs_alignment = form == Format.FASTA or form == Format.UNALIGNED

    if length is None and needs_alignment:
        sys.exit('The length parameter is required for FASTA or unaligned '
                 'input files')

    if needs_alignment:
        if length % 2 == 0:
            sys.exit('The sequence length must be an odd number')

    with open(ifile, 'r') as fh:
        if form == Format.FASTA:
            if is_bg:
                seqs = generate_bg(
                    [s.rstrip() for s in fh if not s.startswith('>')],
                    central_res, length)
            else:
                seqs = [align_sequence(s.rstrip(), central_res, length)
                        for s in fh if not s.startswith('>')]
        elif form == Format.UNALIGNED:
            seqs = [align_sequence(s.rstrip(), central_res, length)
                    for s in fh]
        else:
            seqs = [s.rstrip() for s in fh]

    seqs = list(set(seqs))

    if write and needs_alignment:
        with open(ifile + '.aligned', 'w') as fh:
            fh.write('\n'.join(seqs))

    return seqs


def extract(fg_file, bg_file, central_res, length, min_occs, max_p, verbose=0):
    '''
    '''
    def _extract_seqs(ifile, **kwargs):
        return extract_seqs(ifile, central_res, length, **kwargs)

    fg_seqs = _extract_seqs(fg_file)
    bg_seqs = _extract_seqs(bg_file, is_bg=True)

    motifs = motifx(fg_seqs, bg_seqs, central_res, min_occs, max_p,
                    verbose=verbose)

    return (motifs, fg_seqs)


def parse_args():
    '''
    Parse the command line arguments.

    Returns:
        argparse.Namespace

    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'fg_file',
        help=('The foreground sequence file. Accepts FASTA, prealigned or '
              'unaligned formats.')
    )
    parser.add_argument(
        'bg_file',
        help=('The background sequence file. Accepts FASTA, prealigned or '
              'unaligned formats.')
    )
    parser.add_argument(
        'central_res',
        help='The central residue of the peptide sequences.'
    )
    parser.add_argument(
        'min_occs',
        type=int,
        help=('The minimum number of times a motif must occur in the '
              'foreground.')
    )
    parser.add_argument(
        'p_cutoff',
        type=float,
        help='The maximum allowed p-value.'
    )
    parser.add_argument(
        '-l',
        '--length',
        default=None,
        type=int,
        help='The desired length of the peptide sequences'
    )
    parser.add_argument(
        '-o',
        '--output',
        default='output.csv',
        help='The file to which motif results are written'
    )
    parser.add_argument(
        '-p',
        '--plot',
        action='store_true',
        default=False,
        help='Set to true to plot a sequence logo for each motif'
    )
    parser.add_argument(
        '-v',
        '--verbose',
        action='count',
        default=0,
        help='Enable verbose logging'
    )
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    motifs, fg_seqs = extract(args.fg_file, args.bg_file, args.central_res,
                              args.length, args.min_occs, args.p_cutoff,
                              verbose=args.verbose)

    if motifs is None:
        sys.exit('No motifs found')

    motifs.to_csv(args.output)

    if args.plot:
        sequence_logos.generate_logos(motifs, fg_seqs)
