#! /usr/bin/env python3
'''
A python implementation of the motif-x algorithm, based on rmotifx.

Requires numpy, pandas and scipy.

'''

import math
import re
import sys

import numpy as np
import pandas as pd
from scipy.stats import binom


RESIDUES = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M',
            'F', 'P', 'S', 'T', 'W', 'Y', 'V']


def _check_sequence_lengths(seqs, cat):
    length = len(seqs[0])
    if any([len(seqs[ii]) != length for ii in range(1, len(seqs))]):
        sys.exit('The {} sequences must all have the same length'.format(cat))


def _check_sequences_empty(fg_seqs, bg_seqs):
    if len(fg_seqs) == 0:
        sys.exit('No foreground sequences available')
    if len(bg_seqs) == 0:
        sys.exit('No background sequences available')


def _central_index(seqlen):
    return int(seqlen / 2)


def build_pfm(seqs):
    '''
    Convert the sequences to a position frequency matrix.

    Args:
        seqs (list): The list of peptide sequences.

    Returns:
        np.ndarray

    '''
    # TODO: optimise this function further
    seq_len = len(seqs[0])
    pfm = np.zeros((len(RESIDUES), seq_len))
    seq_df = pd.DataFrame([list(s) for s in seqs])
    for res_idx, res in enumerate(RESIDUES):
        pfm[res_idx, :] = (seq_df == res).sum()

    return pfm


def build_ppm(seqs):
    '''
    Convert the sequences to a position probability matrix.

    Args:
        seqs (list): The list of peptide sequences.

    Returns:
        np.ndarray

    '''
    pfm = build_pfm(seqs)
    return pfm / pfm.sum(axis=0)


def build_motif_pattern(motif, cent_res, seqlen):
    '''
    Construct a pattern for the motif.

    Args:
        motif (list of tuples): A list of tuples,
                                (amino acid residue, position).
        cent_res (str): The central residue(s).
        seqlen (int): The length of the peptide sequences.

    Returns:
        string

    '''
    pat = ['.'] * seqlen
    pat[_central_index(seqlen)] = cent_res
    for res, pos in motif:
        pat[pos] = res
    return ''.join(pat)


def find_motif(fg_seqs, bg_seqs, min_occs, max_p, central_res, verbose=0):
    '''
    Find the motifs in the sequences.

    '''
    fg_size, bg_size = len(fg_seqs), len(bg_seqs)

    seqlen = len(fg_seqs[0])
    central_idx = _central_index(seqlen)

    if verbose > 0:
        print('Step 1: Recursively build motifs')

    motif = []   # A list of index tuples
    pvals = []
    iteration = 1
    while True:
        if len(fg_seqs) == 0 or len(bg_seqs) == 0:
            # There are no sequences remaining
            break

        fg_pfm = build_pfm(fg_seqs)
        bg_ppm = build_ppm(bg_seqs)

        #binomial = 1 - binom.cdf(fg_pfm, len(fg_seqs), bg_ppm, 1)
        binomial = binom.sf(fg_pfm, len(fg_seqs), bg_ppm, 1)

        # TODO: can we revert to the commented version now that R perl
        #       has been replicated?
        #binomial[binomial == 0] = np.finfo(np.float64).tiny
        binomial[binomial < 1e-16] = 1e-16

        # Don't consider the central residue significant when building
        # the motif
        binomial[:, central_idx] = 1.

        # Don't consider the previous motifs significant
        for pos in motif:
            binomial[pos] = 1.

        # Set the p-value to 1 for anything which appears too infrequently
        binomial[fg_pfm < min_occs] = 1.

        # Find the minimum p-value
        min_p = binomial.min()

        minima = np.nonzero((binomial == min_p) & (binomial < max_p))

        if minima[0].size == 0:
            break

        if minima[0].size > 1:
            # Get the first column with a maximum frequency
            # Sort positions by column (emulate R behaviour)
            positions = sorted(zip(minima[1], minima[0]))
            max_freq = fg_pfm[minima].max()
            min_row, min_col = [(r, c) for c, r in positions
                                if fg_pfm[r, c] == max_freq][0]
        else:
            min_row, min_col = minima[0][0], minima[1][0]

        motif.append((min_row, min_col))

        res = RESIDUES[min_row]

        if verbose > 0:
            print('Iteration {}, residue={}, row={}, col={}, p-value={},'
                  'fg-size={}, bg-size={}'.format(
                      iteration, res, min_row, min_col, min_p,
                      len(fg_seqs), len(bg_seqs)))

        fg_seqs = [s for s in fg_seqs if s[min_col] == res]
        bg_seqs = [s for s in bg_seqs if s[min_col] == res]

        pvals.append(min_p)
        iteration += 1

    if not motif:
        return None

    # Convert row index to the amino acid residue
    motif = [(RESIDUES[row], col) for row, col in motif]

    score = sum([-1 * math.log10(p) for p in pvals])

    motif_pat = build_motif_pattern(motif, central_res, seqlen)
    motif_regex = re.compile(motif_pat)

    def get_matches(seqs):
        return [s for s in seqs if motif_regex.match(s)]

    return {
        'score': score,
        'motif': motif_pat,
        'num_fg_matches': len(get_matches(fg_seqs)),
        'num_bg_matches': len(get_matches(bg_seqs)),
        'fg_size': fg_size,
        'bg_size': bg_size
    }


def motifx(fg_seqs, bg_seqs, central_res, min_occs, max_p, verbose=0):
    '''
    Extract statistically significant motifs from the foreground sequences.

    Args:
        fg_seqs (list): The aligned amino acid sequences.
        bg_seqs (list): The aligned background amino acid sequences.
        central_res (str): The residue on which the sequences are centered.
        min_occs (int): The minimum number of times the extracted motifs must
                        occur.
        max_p (float): The p-value threshold for the extracted motifs.
        verbose (int, optional): The logging level.

    Returns:
        pandas.DataFrame

    '''
    # Check that the input sequences are not empty
    _check_sequences_empty(fg_seqs, bg_seqs)

    # Check that the central residue does contain amino acid residues
    if len(set(RESIDUES) & set(central_res)) == 0:
        sys.exit('The central residue sequence must contain at least one '
                 'amino acid residue')

    # Check that the input sequences are of the same length
    _check_sequence_lengths(fg_seqs, 'foreground')
    _check_sequence_lengths(bg_seqs, 'background')

    # Check that the sequence lengths are valid
    fg_len = len(fg_seqs[0])
    if (fg_len < 3 or fg_len > 35):
        sys.exit('The sequence length must be between 3 and 35')
    if fg_len % 2 == 0:
        sys.exit('The sequence length must be an odd number')
    if fg_len != len(bg_seqs[0]):
        sys.exit('The foreground and background sequence lengths must be '
                 'identical')

    # The central index and the length of each side of the sequence
    central_idx = _central_index(fg_len)

    def filter_res(seqs):
        return [s for s in seqs if s[central_idx] in central_res]

    fg_seqs = filter_res(fg_seqs)
    bg_seqs = filter_res(bg_seqs)

    _check_sequences_empty(fg_seqs, bg_seqs)

    motifs = pd.DataFrame()
    while True:
        motif = find_motif(fg_seqs, bg_seqs, min_occs, max_p, central_res,
                           verbose=verbose)

        if motif is None:
            break

        motifs = motifs.append(motif, ignore_index=True)

        if verbose > 0:
            print('Step 2: Reduce positive and negative sets')

        # Remove sequences which are already part of a motif
        motif_regex = re.compile(motif['motif'])
        fg_seqs = [s for s in fg_seqs if not motif_regex.match(s)]
        bg_seqs = [s for s in bg_seqs if not motif_regex.match(s)]

        if len(fg_seqs) < min_occs:
            break

    if verbose > 0:
        print('Converged, no more enrichments to make')

    if motifs.empty:
        return None

    for col in ['fg_size', 'bg_size', 'num_fg_matches', 'num_bg_matches']:
        motifs[col] = motifs[col].astype('int')

    motifs['fold_increase'] = ((motifs.num_fg_matches / motifs.fg_size) /
                               (motifs.num_bg_matches / motifs.bg_size))

    return motifs


if __name__ == '__main__':
    print(motifx(
        ['FAYAA', 'FEYTF', 'WVYYY', 'DAYDA', 'YAYAA', 'DFYPS',
         'NDYTA', 'NTYAE', 'NAYAA'],
        ['YBYAE', 'DAYFY', 'AYYDY', 'AEYTY'], 'Y', 1, 1, verbose=1))
