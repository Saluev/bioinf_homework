import argparse
import sys
from collections import defaultdict

from matplotlib import pyplot as plt

from ngs.io import read_reads


def evaluate_gc_content_distribution(filenames=None, min_base_quality="!", min_good_bases=0,
                                     min_good_bases_fraction=0.0):
    reads = read_reads(filenames)
    hist = defaultdict(int)

    good_reads = 0
    discarded_reads = 0
    for read in reads:
        gc_bases = 0
        good_bases = 0
        total_bases = len(read.sequence)
        for base, quality in zip(read.sequence, read.quality):
            if quality >= min_base_quality:
                gc_bases += (base in "GCgc")
                good_bases += 1
        if good_bases >= min_good_bases and \
                good_bases / total_bases >= min_good_bases_fraction:
            hist[round(gc_bases * 100 / good_bases)] += 1
            good_reads += 1
        else:
            discarded_reads += 1

    total_reads = good_reads + discarded_reads
    print(f"Total reads processed: {total_reads}", file=sys.stderr)
    discarded_reads_fraction = discarded_reads/total_reads
    print(f"Reads discarded due to low quality: {discarded_reads} ({discarded_reads_fraction*100:g}%)", file=sys.stderr)

    xs = range(0, 101)
    ys = [hist[x] for x in xs]
    fig, ax = plt.subplots()
    plt.plot(xs, ys)
    ax.set_xlabel('% GC')
    ax.set_ylabel('Number of reads')
    plt.show()


parser = argparse.ArgumentParser(description='Evaluate GC content distribution.')
parser.add_argument('filenames', metavar='filename', type=str, nargs='*',
                    help='List of files to process. If empty, stdin will be processed.')
parser.add_argument('--min-base-quality', type=str, default='!', metavar="!...~",
                    help='Only take into account nucleotides with quality above or equal to specified.')
parser.add_argument('--min-good-bases', type=int, default=0, metavar="N",
                    help='Only take into account reads with ≥ N good nucleotides.')
parser.add_argument('--min-good-bases-fraction', type=float, default=0, metavar="X",
                    help='Only take into account reads with ≥ X%% good nucleotides.')


if __name__ == "__main__":
    args = parser.parse_args()
    evaluate_gc_content_distribution(**args.__dict__)
