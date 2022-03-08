import argparse

import numpy as np
from matplotlib import pyplot as plt

from ngs.io import read_fasta, read_sam, open_with_gzip


def evaluate_alignment_coverage(reference: str, alignment: str, window: int = 1000):
    contigs = {contig.identifier: contig for contig in read_fasta(open_with_gzip(reference))}
    plots = {}
    for contig in contigs.values():
        plots[contig.identifier] = np.zeros(len(contig.sequence) + 1)

    for alignment in read_sam(open_with_gzip(alignment)):
        if alignment.reference == "*" or alignment.sequence is None:
            continue
        for i in range(min(len(alignment.sequence), len(plots[alignment.reference]) - alignment.reference_position)):
            plots[alignment.reference][alignment.reference_position+i] += 1

    for plot in plots.values():
        cumsum = np.cumsum(np.insert(plot, 0, 0))
        moving_average = (cumsum[window:] - cumsum[:-window]) / window
        plt.plot(moving_average)
    plt.show()
    return {contig_name: np.mean(plot) for contig_name, plot in plots.items()}


parser = argparse.ArgumentParser(description='Evaluate alignment coverage.')
parser.add_argument('reference', type=str, help='Path to reference genome.')
parser.add_argument('alignment', type=str, help='Path to alignment.')


if __name__ == "__main__":
    args = parser.parse_args()
    evaluate_alignment_coverage(**args.__dict__)
