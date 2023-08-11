#!/usr/bin/env python
import os

import pytest

from run_ExomeDepth import get_sampleid
from run_ExomeDepth import get_vcf_stats
from run_ExomeDepth import write_log_stats


def rootdir():
    return os.path.dirname(os.path.abspath(__file__))


# Test get_sampleid function
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
    sample_id = get_sampleid(sample_id)
    assert sample_id == expected


# Test get_vcf_stats function
VCF_test_file = os.path.join(
    rootdir(), 'HC_RS-SSv7-2023-1_U224111CM2023D06051_230324_A01131_0312_AHGK2MDMXY_NICU_U224111_exome_calls.vcf'
)
assert get_vcf_stats(VCF_test_file) == (0.9923, 56.98, 86)


# Test write_log_stats function
def test_write_logstats(tmp_path):
    path = tmp_path / "logstats_folder"
    path.mkdir()
    logstats = path / "HC_TestSample_stats.log"
    message = "TestSample\tHC\tRefset1\tfemale\t0.9923\t56.98\t86\tWARNNG:ThisIsATest"

    write_log_stats(
        "_", "exome_calls", "\tWARNNG:ThisIsATest", 0.9923, 56.98,
        86, path, "HC", "TestSample", "Refset1", "female"
    )

    assert message in logstats.read_text()
