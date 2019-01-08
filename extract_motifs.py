#! /usr/bin/env python3
'''
Extract statistically significant motifs from the input peptide sequences.

'''
import enum
import sys

import click

from align import align_sequence
from pymotifx import motifx


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
    
    
@click.command()
@click.argument('foreground')
@click.argument('background')
@click.argument('central-residue')
@click.option('-m', '--min-occurences', 'min_occs', type=int, default=1,
              help='The minimum number of times a pattern/motif must appear '
                   'in the foreground data.', show_default=True)
@click.option('-p', '--p-cutoff', type=float, default=1e-6,
               help='The maximum allowed p-value.', show_default=True)
@click.option('-l', '--length', type=int,
               help='The target length of the peptide sequences.')
@click.option('-m', '--motif-output', default='output.csv',
               help='The file to which motif results are written.')
@click.option('-o', '--logo-output', default='logos',
              help='The directory to which to write logo information')
@click.option('--plot/--no-plot', default=False,
              help='If enabled, plot a sequence logo for each motif.')
@click.option('-v', '--verbose', count=True, help='Enable verbose logging.')
def main(foreground, background, central_residue, min_occs,
         p_cutoff, length, motif_output, logo_output, plot, verbose):
    '''
    Extract statistically significant motifs from the foreground sequences.

    Args:
        foreground (str): The foreground sequence file, in FASTA,
                          prealigned or unaligned format.

        background (str): The background sequence file, in FASTA,
                          prealigned or unaligned format.

        central-residue (str): The residue at the centre of the peptide
                               sequences.

    '''
    motifs, fg_seqs = extract(foreground, background, central_residue, length,
                              min_occs, p_cutoff, verbose=verbose)
    
    if motifs is None:
        sys.exit('No motifs found')
        
    motifs.to_csv(motif_output)
    
    if plot:
        import sequence_logos
        sequence_logos.generate_logos(motifs, fg_seqs, logo_output)


if __name__ == '__main__':
    main()
