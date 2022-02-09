import io
from typing import Iterable

import pytest

from ngs.io import read_fastq, Read, Sequence, SequenceQuality


def test_read_fastq():
    with pytest.raises(ValueError):
        exhaust(read_fastq(io.StringIO('wrong')))

    with pytest.raises(ValueError):
        exhaust(read_fastq(io.StringIO('@wrong\nAAAA\n+')))

    with pytest.raises(ValueError):
        exhaust(read_fastq(io.StringIO('@wrong\nACAB\n+\nAAAA')))

    assert list(read_fastq(io.StringIO('@id desc\nACAC\n+\nAAAA'))) == [
        Read(identifier='id', description='desc', sequence=Sequence('ACAC'), quality=SequenceQuality('AAAA')),
    ]

    assert list(read_fastq(io.StringIO('@id desc\nACAC\n+\nAAAA\n@di csed\nATCG\n+di\n!?#$'))) == [
        Read(identifier='id', description='desc', sequence=Sequence('ACAC'), quality=SequenceQuality('AAAA')),
        Read(identifier='di', description='csed', sequence=Sequence('ATCG'), quality=SequenceQuality('!?#$')),
    ]


def exhaust(it: Iterable):
    for _ in it:
        pass
