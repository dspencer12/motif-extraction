#! /usr/bin/env python3
'''
Functions and classes for use in motif extraction.

'''
import enum
import sys


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
        length = len(first)
        res = first[int(length / 2)]
        for line in fh:
            line = line.rstrip()
            if len(line) != length or line[int(len(line) / 2)] != res:
                return Format.UNALIGNED

    return Format.PREALIGNED


def parse_fasta(fh):
    '''
    Parse the input fasta file handle to extract the sequences.

    Args:
        fh (file)

    '''
    seq = ''
    for line in fh:
        if line.startswith('>') and seq:
            yield seq
            seq = ''
        if not line.startswith('>'):
            seq += line.rstrip('\n')


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
                seqs = generate_bg(parse_fasta(fh), central_res, length)
            else:
                seqs = [align_sequence(s, central_res, length)
                        for s in parse_fasta(fh)]
        elif form == Format.UNALIGNED:
            seqs = [align_sequence(s.rstrip(), central_res, length)
                    for s in fh]
        else:
            seqs = [s.rstrip() for s in fh]

    seqs = [s for s in list(set(seqs)) if s]

    if write and needs_alignment:
        with open(ifile + '.aligned', 'w') as fh:
            fh.write('\n'.join(seqs))

    return seqs


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

    return [s for s in bg if s]


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
        return ''

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
