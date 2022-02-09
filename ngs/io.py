import gzip
import re
import sys
from dataclasses import dataclass
from itertools import zip_longest
from typing import Optional, TextIO, Iterable, List

from tqdm import tqdm


class Sequence(str):
    pass


class SequenceQuality(str):
    pass


@dataclass
class Read:
    identifier: str
    description: str
    sequence: Sequence
    quality: Optional[SequenceQuality]


sequence_regexp = re.compile(r'^([ATCGUN]+|[ACDEFGHIKLMNPQRSTVWYZ]+)$', re.IGNORECASE)
quality_regexp = re.compile(r'^[!-~]+$')


def read_fastq(f: TextIO) -> Iterable[Read]:
    it = iter(f)
    for identifier, sequence, delimiter, quality in zip_longest(it, it, it, it):
        if quality is None:
            raise ValueError('number of lines in the file is not divisible by 4')
        if not identifier.startswith('@'):
            raise ValueError(f'expected @READ_IDENTIFIER, got {identifier!r}')
        if not delimiter.startswith('+'):
            raise ValueError(f'expected +..., got {delimiter!r}')
        sequence = sequence.rstrip('\n\r')
        quality = quality.rstrip('\n\r')
        if not sequence_regexp.fullmatch(sequence):
            raise ValueError(f'expected nucleotide/amino acid sequence, got {sequence!r}')
        if not quality_regexp.fullmatch(quality):
            raise ValueError(f'expected quality string, got {quality!r}')
        if len(sequence) != len(quality):
            raise ValueError(f'sequence and quality strings length mismatch: {len(sequence)} != {len(quality)}')
        identifier, *description = identifier.removeprefix('@').rstrip('\n\r').split(maxsplit=1)
        yield Read(
            identifier=identifier,
            description=description[0] if description else '',
            sequence=sequence,
            quality=quality,
        )


def open_with_gzip(filename: str) -> TextIO:
    if filename.endswith(".gz"):
        return gzip.open(filename, mode="rt")
    return open(filename)


def read_reads(filenames: Optional[List[str]]) -> Iterable[Read]:
    if not filenames:
        print(f"Processing reads from stdin...", file=sys.stderr)
        yield from tqdm(read_fastq(sys.stdin), desc="Reads processed", unit="")
        print(f"Finished processing reads from stdin", file=sys.stderr)
    for filename in filenames:
        print(f"Processing reads from file {filename}...", file=sys.stderr)
        yield from tqdm(read_fastq(open_with_gzip(filename)), desc="Reads processed", unit="")
        print(f"Finished processing file {filename}.", file=sys.stderr)
