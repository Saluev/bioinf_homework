import argparse
import collections
import dataclasses
from typing import Mapping, Tuple

from ngs.io import read_fasta, read_sam, open_with_gzip


@dataclasses.dataclass
class Stats:
    replacement_count: Mapping[Tuple[str, str], int]
    avg_error_rate: float
    avg_insertion_length: float
    avg_deletion_length: float
    homopolymer_close_insertions_rate: float
    avg_insertion_quality: float
    avg_mismatch_quality: float


def evaluate_mismatch_frequency(reference: str, alignment: str, homopolymer_length_threshold=10) -> Stats:
    contigs = {contig.identifier: contig for contig in read_fasta(open_with_gzip(reference))}
    stats = collections.defaultdict(int)

    error_rates = []
    insertion_lengths = []
    deletion_lengths = []
    insertion_qualities = []
    mismatch_qualities = []

    homopolymer_close_insertions = 0
    total_insertions = 0

    def count_mismatches(fr: str, to: str, qs: str):
        nonlocal curr_wrong_positions, curr_total_positions
        for l1, l2, q in zip(fr, to, qs):
            if l1 != l2:
                curr_wrong_positions += 1
                stats[l1, l2] += 1
                mismatch_qualities.append(ord(q) - ord('!'))
        curr_total_positions += len(fr)

    def count_insertions(to: str, qs: str):
        nonlocal curr_wrong_positions, curr_total_positions
        for l2, q in zip(to, qs):
            stats["", l2] += 1
            insertion_qualities.append(ord(q) - ord('!'))
        curr_wrong_positions += len(to)
        curr_total_positions += len(to)
        insertion_lengths.append(len(to))

    def count_deletions(fr: str):
        nonlocal curr_wrong_positions, curr_total_positions
        for l1 in fr:
            stats[l1, ""] += 1
        curr_wrong_positions += len(fr)
        curr_total_positions += len(fr)
        deletion_lengths.append(len(fr))

    def homopolymer_close_to(contig_sequence: str, position: int) -> bool:
        t = homopolymer_length_threshold
        if position-t >= 0 and len(set(contig_sequence[position-t:position])) == 1:
            return True
        if position+t < len(contig_sequence) and len(set(contig_sequence[position:position+t])) == 1:
            return True
        return False

    for alignment in read_sam(open_with_gzip(alignment)):
        if alignment.reference == "*":
            continue
        contig = contigs[alignment.reference]
        alignment_offset = 0
        reference_offset = alignment.reference_position
        curr_wrong_positions = 0
        curr_total_positions = 0
        for n, letter in alignment.cigar.items:
            if letter in "MX=":
                count_mismatches(
                        contig.sequence[reference_offset:reference_offset+n],
                        alignment.sequence[alignment_offset:alignment_offset+n],
                        alignment.quality[alignment_offset:alignment_offset+n],
                )
                alignment_offset += n
                reference_offset += n
            elif letter == "I":
                count_insertions(
                    alignment.sequence[alignment_offset:alignment_offset+n],
                    alignment.quality[alignment_offset:alignment_offset+n],
                )
                alignment_offset += n

                homopolymer_close_insertions += homopolymer_close_to(contig.sequence, reference_offset)
                total_insertions += 1
            elif letter == "S":
                alignment_offset += n
            elif letter == "D":
                count_deletions(contig.sequence[reference_offset:reference_offset+n])
                reference_offset += n
            elif letter == "N":
                reference_offset += n
        if curr_total_positions > 0:
            error_rates.append(curr_wrong_positions / curr_total_positions)

    avg_error_rate = sum(error_rates) / len(error_rates)
    avg_insertion_length = sum(insertion_lengths) / len(insertion_lengths)
    avg_deletion_length = sum(deletion_lengths) / len(deletion_lengths)
    homopolymer_close_insertions_rate = homopolymer_close_insertions/total_insertions if total_insertions else 0
    return Stats(
        replacement_count=stats,
        avg_error_rate=avg_error_rate,
        avg_insertion_length=avg_insertion_length,
        avg_deletion_length=avg_deletion_length,
        homopolymer_close_insertions_rate=homopolymer_close_insertions_rate,
        avg_insertion_quality=sum(insertion_qualities) / len(insertion_qualities),
        avg_mismatch_quality=sum(mismatch_qualities) / len(mismatch_qualities),
    )


parser = argparse.ArgumentParser(description='Evaluate alignment coverage.')
parser.add_argument('reference', type=str, help='Path to reference genome.')
parser.add_argument('alignment', type=str, help='Path to alignment.')


if __name__ == "__main__":
    args = parser.parse_args()
    print(evaluate_mismatch_frequency(**args.__dict__))
