import argparse
import sys
from collections import defaultdict

from matplotlib import pyplot as plt

from ngs.io import read_reads


def evaluate_quality_distribution(filenames=None):
    reads = read_reads(filenames)
    sum_per_position = defaultdict(int)
    count_per_position = defaultdict(int)

    total_reads = 0
    for read in reads:
        for position, quality in enumerate(read.quality):
            sum_per_position[position] += ord(quality) - ord('!')
            count_per_position[position] += 1
        total_reads += 1

    print(f"Total reads processed: {total_reads}", file=sys.stderr)

    max_position = max(sum_per_position.keys())

    xs = range(0, max_position + 1)
    ys = [sum_per_position[x] / count_per_position[x] for x in xs]
    fig, ax = plt.subplots()
    plt.plot(xs, ys)
    ax.set_xlabel('Position in read')
    ax.set_ylabel('Average quality of bases')
    plt.show()


parser = argparse.ArgumentParser(description='Evaluate base quality distribution.')
parser.add_argument('filenames', metavar='filename', type=str, nargs='*',
                    help='List of files to process. If empty, stdin will be processed.')


if __name__ == "__main__":
    args = parser.parse_args()
    evaluate_quality_distribution(**args.__dict__)
