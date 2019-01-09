#! /usr/bin/env python3
'''
Extract statistically significant motifs from the input peptide sequences.

'''
import sys

import click

from pymotifx import motifx
from utilities import extract_seqs


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
