#! /usr/bin/env python3

import pysam

from database.database import connect_database
from database.models import Sample


def add_sample_flowcell_to_db(sample_id, flowcell_id, refset, print_message=None):
    Session = connect_database()
    with Session() as session:
        sample = session.query(Sample).filter(Sample.sample == sample_id).filter(Sample.flowcell == flowcell_id).one_or_none()
        if not sample:
            session.add(Sample(sample=sample_id, flowcell=flowcell_id, refset=refset))
            session.commit()
            if print_message:
                print("## Sample {0} with flowcell_id {1} and with refset {2} added to database".format(
                   sample_id, flowcell_id, refset)
                )                
            return refset
        else:
            if print_message:
                print("## Sample {0} with flowcell_id {1} and with refset {2} already in database".format(
                   sample_id, flowcell_id, refset)
                )
            return sample.refset


def add_sample_to_db(flowcell_id, sample_id, refset, print_message=None):
    flowcell_id = get_flowcell_id(flowcell_id)
    add_sample_flowcell_to_db(sample_id, flowcell_id, refset, print_message)


def add_sample_to_db_and_return_refset_bam(bam, refset, print_refset_stdout=None):
    sample_id = get_sample_id(bam)
    flowcell_id = get_flowcell_id_bam(bam)
    refset_db = add_sample_flowcell_to_db(sample_id, flowcell_id, refset)

    if print_refset_stdout:
        print(refset_db)

    return refset_db


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
            print("## Changed refset of sample {0} with flowcell_id {1} to refset {2}".format(
                sample_id, flowcell_id, refset)
            )
        else:
            print("## Sample {0} with flowcell_id {1} not in refset database".format(
                sample_id, flowcell_id)
            )


def print_all_samples():
    Session = connect_database()
    with Session() as session:
        print("Name\tFlowcell\tRefset\tFamilyID")
        for item in session.query(Sample):
            print("{0}\t{1}\t{2}".format(item.sample, item.flowcell, item.refset))


def print_refset(flowcell_id, sample_id):
    Session = connect_database()
    with Session() as session:
        sample = session.query(Sample).filter(Sample.sample == sample_id).filter(Sample.flowcell == flowcell_id).one_or_none()
        if sample:
            print(sample.refset)
        else:
            print("## Sample {0} with flowcell_id {1} not detected in database".format(
                sample_id, flowcell_id)
            )

def query_refset(flowcell_id, sample_id):
    flowcell_id = get_flowcell_id(flowcell_id)
    print_refset(flowcell_id, sample_id)


def query_refset_bam(bam):
    flowcell_id = get_flowcell_id_bam(bam)
    sample_id = get_sample_id(bam)
    print_refset(flowcell_id, sample_id)


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
            print("## Deleted sample {0} with flowcell_id {1}.".format(sample_id, flowcell_id))
        else:
            print("## Sample {0} with flowcell_id {1} not in database.".format(sample_id, flowcell_id))


def get_flowcell_id(flowcells_arg):
    return "_".join(sorted(set(flowcells_arg)))


def return_refset_bam(bam):
    Session = connect_database()
    sample_id = get_sample_id(bam)
    flowcell_id = get_flowcell_id_bam(bam)
    with Session() as session:
        sample =  session.query(Sample).filter(Sample.sample == sample_id).filter(Sample.flowcell == flowcell_id).one_or_none()
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
