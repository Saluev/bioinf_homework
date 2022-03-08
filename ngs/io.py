import csv
import gzip
import re
import sys
from dataclasses import dataclass
from itertools import zip_longest
from typing import Optional, TextIO, Iterable, List, Tuple

from tqdm import tqdm


class Sequence(str):
    pass


cigar_regexp = re.compile(r"(\d+)(\w)")


class CIGAR(str):
    @property
    def items(self) -> Iterable[Tuple[int, str]]:
        for n, l in cigar_regexp.findall(self):
            yield int(n), l


class SequenceQuality(str):
    pass


@dataclass
class Read:
    identifier: str
    description: str
    sequence: Sequence
    quality: Optional[SequenceQuality] = None


@dataclass
class Alignment:
    identifier: str
    reference: str
    reference_position: int  # 0-based
    cigar: CIGAR
    sequence: Optional[Sequence]
    quality: Optional[SequenceQuality] = None


sequence_regexp = re.compile(r'^([ATCGUN]+|[ACDEFGHIKLMNPQRSTVWYZ]+)$', re.IGNORECASE)
quality_regexp = re.compile(r'^[!-~]+$')


def read_fasta(f: TextIO) -> Iterable[Read]:
    current_identifier = None
    current_description = ""
    current_sequence = []
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_sequence:
                yield Read(
                    identifier=current_identifier,
                    description=current_description,
                    sequence=Sequence("".join(current_sequence)),
                )
            current_identifier, *current_description = line.removeprefix('>').rstrip('\n\r').split(maxsplit=1)
            current_sequence = []
        else:
            if current_identifier is None:
                raise ValueError(f'expected >identifier, got {line!r}')
            if not sequence_regexp.fullmatch(line):
                raise ValueError(f'expected nucleotide/amino acid sequence, got {line!r}')
            current_sequence.append(line)
    if current_sequence:
        yield Read(
            identifier=current_identifier,
            description=current_description,
            sequence=Sequence("".join(current_sequence)),
        )


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


def read_sam(f: TextIO) -> Iterable[Alignment]:
    for i, row in enumerate(csv.reader(f, delimiter="\t", quoting=csv.QUOTE_NONE)):
        if row[0].startswith("@"):
            continue
        if len(row) < 10:
            print(row)
        sequence = Sequence(row[9])
        if sequence == "*":
            sequence = None
        elif not sequence_regexp.fullmatch(sequence):
            raise ValueError(f'expected nucleotide/amino acid sequence, got {sequence!r} in row {i}')
        quality = row[10]
        if not quality_regexp.fullmatch(quality):
            raise ValueError(f'expected quality string, got {quality!r} in row {i}')
        yield Alignment(
            identifier=row[0],
            reference=row[2],
            reference_position=int(row[3])-1,
            cigar=CIGAR(row[5]),
            sequence=sequence,
            quality=SequenceQuality(quality),
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
