#! /usr/bin/env python3

import pysam
import glob

from database.database import connect_database
from database.models import Sample


def add_sample_flowcell_to_db(sample_id, flowcell_id, refset):
    Session = connect_database()
    with Session() as session:
        sample = session.query(Sample).filter(Sample.sample == sample_id).filter(Sample.flowcell == flowcell_id).one_or_none()
        if not sample:
            session.add(Sample(sample=sample_id, flowcell=flowcell_id, refset=refset))
            session.commit()
            return refset, True
        else:
            return sample.refset, False


def add_sample_to_db(flowcell_id, sample_id, refset):
    flowcell_id = get_flowcell_id(flowcell_id)
    refset_db, added = add_sample_flowcell_to_db(sample_id, flowcell_id, refset)
    return flowcell_id, sample_id, refset_db, added


def add_sample_to_db_and_return_refset_bam(bam, refset):
    sample_id = get_sample_id(bam)
    flowcell_id = get_flowcell_id_bam(bam)
    refset_db, added = add_sample_flowcell_to_db(sample_id, flowcell_id, refset)
    return flowcell_id, sample_id, refset_db, added


def change_refset_in_db(flowcell_id, sample_id, refset):
    Session = connect_database()
    with Session() as session:
        sample_update = (
            session.query(Sample)
            .filter(Sample.sample == sample_id)
            .filter(Sample.flowcell == flowcell_id)
            .one_or_none()
        )
        if sample_update:
            sample_update.refset = refset
            session.add(sample_update)
            session.commit()
            return flowcell_id, sample_id, refset, True

        else:
            return flowcell_id, sample_id, refset, False


def return_all_samples():
    sample_list = []
    Session = connect_database()
    with Session() as session:
        for item in session.query(Sample):
            sample_list.append("{0}\t{1}\t{2}".format(item.sample, item.flowcell, item.refset))
    return sample_list


def parse_refset(flowcell_id, sample_id):
    Session = connect_database()
    with Session() as session:
        sample = session.query(Sample).filter(Sample.sample == sample_id).filter(Sample.flowcell == flowcell_id).one_or_none()
        if sample:
            return sample.refset
        else:
            return None


def query_refset(flowcell_id, sample_id):
    flowcell_id = get_flowcell_id(flowcell_id)
    return parse_refset(flowcell_id, sample_id), flowcell_id


def query_refset_bam(bam):
    flowcell_id = get_flowcell_id_bam(bam)
    sample_id = get_sample_id(bam)
    return parse_refset(flowcell_id, sample_id), flowcell_id, sample_id


def get_flowcell_id_bam(bam):
    with pysam.AlignmentFile(bam, "rb") as workfile:
        readgroups = []
        for readgroup in workfile.header['RG']:
            if readgroup['PU'] not in readgroup:
                readgroups.append(readgroup['PU'])
    return "_".join(sorted(set(readgroups)))


def delete_sample_db(flowcell_id, sample_id):
    Session = connect_database()
    flowcell_id = get_flowcell_id(flowcell_id)
    with Session() as session:
        sample = session.query(Sample).filter(Sample.sample == sample_id).filter(Sample.flowcell == flowcell_id).one_or_none()
        if sample:
            session.delete(sample)
            session.commit()
            return True, flowcell_id
        else:
            return False, flowcell_id


def get_flowcell_id(flowcells_arg):
    return "_".join(sorted(set(flowcells_arg)))


def return_refset_bam(bam):
    Session = connect_database()
    sample_id = get_sample_id(bam)
    flowcell_id = get_flowcell_id_bam(bam)
    with Session() as session:
        sample = session.query(Sample).filter(Sample.sample == sample_id).filter(Sample.flowcell == flowcell_id).one_or_none()
        if sample:
            return sample.refset
        else:
            return "refset_unknown"


def get_sample_id(bam):
    with pysam.AlignmentFile(bam, "rb") as workfile:
        sampleid = []
        for readgroup in workfile.header['RG']:
            sampleid.append(readgroup['SM'])
        sampleid = list(set(sampleid))
        sampleid = "_".join(sampleid)
    return sampleid


def fill_database(path, conflicts):
    folders = set(glob.glob("{}*".format(path), recursive=True))
    print("Number of folders detected = {}".format(len(folders)))
    qc_files = glob.glob("{}/QC/CNV/*exomedepth_summary.txt".format(path), recursive=True)

    bam_files = []
    qc_files = {}
    for folder in folders:
        qc_file = glob.glob("{}/QC/CNV/*exomedepth_summary.txt".format(folder), recursive=True)
        if not qc_file:
            print("WARNING: CNV QC file missing folder {} Skipping all samples in folder!".format(folder))
            continue

        if len(qc_file) > 1:
            print("WARNING: Multiple QC files detected in folder {} Skipping all samples in folder!".format(folder))
            continue

        if folder not in qc_files:
            qc_files[folder] = qc_file[0]
        else:
            "WARNING: folder {} has been detected multiple times. Skipping all samples in folder!".format(folder)

        for nf_bam in glob.glob("{}/bam_files/*.bam".format(folder), recursive=True): #  Nextflow analysis
            bam_files.append(nf_bam)
        for iap_bam in glob.glob("{}/*/mapping/*realigned.bam".format(folder), recursive=True): # IAP analysis
            bam_files.append(iap_bam)

    print("Number of folders with qc_file detected = {}".format(len(qc_files)))
    print("Number of BAM files detected = {}".format(len(bam_files)))


    # Parse reference set for each sample from CNV QC summary file
    sample_refset = parse_refset_qc_files(qc_files)

    # Fill reference database for all samples (bam files) and resolve conflict
    conflicts = add_database_bam(bam_files, sample_refset)
    for conflict in conflicts:
        print(conflict)


def parse_refset_qc_files(qc_files):
    sample_refset = {}
    for qc_file in qc_files:
        refset_qc = open(qc_files[qc_file], "r").readlines()
        for line in refset_qc:
            if "REFSET" in line:
                splitline = line.rstrip().split(";")

                #  Determine index position of REFSET=
                header = []
                for item in splitline:
                    header.append(item.split("=")[0])
                refset_index = header.index('REFSET')

                #  Parse used reference set from QC file into dictionary
                if splitline[0] not in sample_refset:
                    if "WARNING" in line:
                        sample_refset[splitline[0]] = [splitline[refset_index].replace("REFSET=", ""), line.rstrip().split("\t")[1:]]
                    else:
                        sample_refset[splitline[0]] = [splitline[refset_index].replace("REFSET=", "")]

    return sample_refset


def add_database_bam(bam_files, sample_refset):
    conflicts = []
    for bam in bam_files:
        sample_id = get_sample_id(bam)
        flowcell_id = get_flowcell_id_bam(bam)
        refset_db, added = add_sample_flowcell_to_db(sample_id, flowcell_id, sample_refset[sample_id][0])
        if len(sample_refset[sample_id]) > 1:
            warning = sample_refset[sample_id][1]
        else:
            warning = "None"

        if not added:
            conflicts.append([sample_id, flowcell_id, sample_refset[sample_id], refset_db, warning])

    return conflicts

