#!/usr/bin/env python
import pytest

from parse_cnv_qc import get_number_samples_per_run_from_samplesheet

def test_get_number_samples_per_run_from_samplesheet_withsamples():
    folder = "./run1"
    rawfolder_run = {"./run1": [["run1"], 0]}
    projects = ["CREV4", "NICU"]
    assert get_number_samples_per_run_from_samplesheet(folder, rawfolder_run, projects, [])[0] == {"./run1": [["run1"], 2]}

def test_get_number_samples_per_run_from_samplesheet_nosamples():
    folder = "./run1"
    rawfolder_run = {"./run1": [["run1"], 0]}
    projects = ["WGS", "RNASEQ"]
    assert get_number_samples_per_run_from_samplesheet(folder, rawfolder_run, projects, [])[0] == {"./run1": [["run1"], 0]}

def test_get_number_samples_per_run_from_samplesheet_nosamplesheet():
    folder = "./run2"
    rawfolder_run = {"./run2": [["run2"], 0]}
    projects = ["CREV4", "NICU"]
    assert get_number_samples_per_run_from_samplesheet(folder, rawfolder_run, projects, [])[0] == {"./run2": [["run2"], 0]}

def test_get_number_samples_per_run_from_samplesheet_nowarning():
    folder = "./run1"
    rawfolder_run = {"./run1": [["run1"], 0]}
    projects = ["CREV4", "NICU"]
    assert get_number_samples_per_run_from_samplesheet(folder, rawfolder_run, projects, [])[1] == []

def test_get_number_samples_per_run_from_samplesheet_with_warning_nosamplsheet():
    folder = "./run2"
    rawfolder_run = {"./run2": [["run2"], 0]}
    projects = ["CREV4", "NICU"]
    assert get_number_samples_per_run_from_samplesheet(folder, rawfolder_run, projects, [])[1] == [
        f"no samplesheet for run {folder}, assuming unknown number of samples in run"
    ]

