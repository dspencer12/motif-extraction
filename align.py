#! /usr/bin/env python3

import sys


def align_sequence(seq, central_res, length, extra_char='X'):
    '''
    Align the given sequence on the central residue.

    Args:
        seq (str): The sequence to be aligned.
        central_res (str): The central residue.
        length (int): The desired length of each sequence.
        extra_char (str, optional): The character to be used to make up
                                    the sequence length.

    Returns:
        str: The aligned sequence

    '''
    if length is None:
        sys.exit('Cannot align sequence without specifying the length')

    hlen = int(0.5 * length)

    res_pos = [ii for ii, res in enumerate(seq) if res == central_res]

    if hlen in res_pos:
        split_index = hlen
    else:
        # Find the instance of central_res closest to hlen
        min_dist, split_index = None, None
        for index in res_pos:
            dist = abs(hlen - index)
            if min_dist is None or dist < min_dist:
                min_dist, split_index = dist, index

    if split_index is None:
        print(seq)

    left, right = seq[:split_index], seq[split_index+1:]

    if len(left) < hlen:
        left = extra_char * (hlen - len(left)) + left
    elif len(left) > hlen:
        left = left[len(left) - hlen:]
    if len(right) < hlen:
        right += extra_char * (hlen - len(right))
    elif len(right) > hlen:
        right = right[:hlen]

    return central_res.join([left, right])
