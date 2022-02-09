import gzip
import re
from dataclasses import dataclass
from itertools import zip_longest
from typing import Optional, TextIO, Iterable


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
