#!/usr/bin/env python
import pytest

from make_BEDdetail import get_familystatus_sampleid


@pytest.mark.parametrize("sample_id, expected", [
    ("U175754CM2022D12878", "child"),
    ("U175754CF2022D12878", "child"),
    ("U175754CO2022D12878", "child"),
    ("U175754PM2022D12878", "parent"),
    ("U175754PF2022D12878", "parent"),
    ("U175754PO2022D12878", "parent"),
    ("U175754CC2022D12878", "unknown"),
    ("2022D12878", "unknown"),
])
def test_familystatus_sampleid(sample_id, expected):
    sample_id = get_familystatus_sampleid(sample_id)
    assert sample_id == expected
