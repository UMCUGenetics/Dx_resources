#!/usr/bin/env python
import sys
import pytest
import re

sys.path.append("../../ExomeDepth")
from run_ExomeDepth import split_sampleid


@pytest.mark.parametrize("sample_id, expected", [
    ("U175754CM2022D12878", "2022D12878"),
    ("U175754CF2022D12878", "2022D12878"),
    ("U175754PM2022D12878", "2022D12878"),
    ("U175754PF2022D12878", "2022D12878"),
    ("U175754CO2022D12878", "2022D12878"),
    ("U175754CC2022D12878", "U175754CC2022D12878"),
    ("2022D12878", "2022D12878"),
])


def test_split_sampleid(sample_id, expected):
    sample_id = split_sampleid(sample_id) 
    assert sample_id == expected
