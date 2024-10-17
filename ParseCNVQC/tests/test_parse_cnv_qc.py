#!/usr/bin/env python
import pytest

from parse_cnv_qc import include_sample_number


def test_include_sample_number():
    folder1 = "./run1"
    folder2 = "./run2"
    folder3 = "./run3"
    rawfolder_run1 = {"./run1": [["run1"], 0]}
    rawfolder_run2 = {"./run2": [["run2"], 0]}
    rawfolder_run3 = {"./run3": [["run3"], 0]}
    projects1 = ["CREV4", "NICU"]
    projects2 = ["WGS"]

    assert include_sample_number(folder1, rawfolder_run1, projects1, [])[0] == {"./run1": [["run1"], 2]}
    assert include_sample_number(folder1, rawfolder_run1, projects1, [])[1] == []
    assert include_sample_number(folder2, rawfolder_run2, projects1, [])[0] == {"./run2": [["run2"], 0]}
    assert include_sample_number(folder2, rawfolder_run2, projects1, [])[1] == [
        f"no samplesheet for run {folder2}, assuming unknown number of samples in run"
    ]
    assert include_sample_number(folder3, rawfolder_run3, projects2, [])[0] == {"./run3": [["run3"], 0]}
