#! /usr/bin/env python3

from database.database import connect_database
from database.models import Sample

import pysam


def add_sample_to_db(flowcell_id, sample_id, refset):
    Session = connect_database()
    flowcell_id = get_flowcell_id(flowcell_id)
    with Session() as session:
        if not session.query(Sample).filter(Sample.sample == sample_id).filter(Sample.flowcell == flowcell_id).all():
            sample = Sample(sample=sample_id, flowcell=flowcell_id, refset=refset)
            session.add(sample)
            session.commit()
            print("## Sample {0} added to database with flowcell_id {1} and with refset {2}".format(
                sample_id, flowcell_id, refset)
            )


def change_refset_in_db(flowcell_id, sample_id, refset):
    Session = connect_database()
    with Session() as session:
        sample_update = (
            session.query(Sample)
            .filter(Sample.sample == sample_id)
            .filter(Sample.flowcell == flowcell_id)
            .one()
        )
        sample_update.refset = refset
        session.add(sample_update)
        session.commit()
        print("## Changed refset of sample {0} with flowcell_id {1} to refset {2}".format(
            sample_id, flowcell_id, refset)
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
        if session.query(Sample).filter(Sample.sample == sample_id).filter(Sample.flowcell == flowcell_id).all():
            print(
                session.query(Sample)
                .filter(Sample.sample == sample_id)
                .filter(Sample.flowcell == flowcell_id)
                .one()
                .refset
            )
        else:
            print("## Sample {0} with flowcell_id {1} not detected in database".format(
                sample_id, flowcell_id)
            )


def query_refset(flowcell_id, sample_id):
    flowcell_id = get_flowcell_id(flowcell_id)
    print_refset(flowcell_id, sample_id)


def query_refset_bam(bam):
    flowcell_id = get_flowcelid_bam(bam)
    sample_id = get_sample_id(bam)
    print_refset(flowcell_id, sample_id)


def get_flowcelid_bam(bam):
    with pysam.AlignmentFile(bam, "rb") as workfile:
        readgroups = []
        for readgroup in workfile.header['RG']:
            if readgroup['PU'] not in readgroup:
                readgroups.append(readgroup['PU'])
        readgroups = list(set(readgroups))
        flowcell_id = "_".join(readgroups)
    return flowcell_id


def delete_sample_db(flowcell_id, sample_id):
    Session = connect_database()
    flowcell_id = get_flowcell_id(flowcell_id)
    with Session() as session:
        samples = session.query(Sample).filter(Sample.sample == sample_id).filter(Sample.flowcell == flowcell_id)
        if samples:
            samples.delete()
            session.commit()
            print("Deleted sample {0} with flowcell_id {1}.".format(sample_id, flowcell_id))
        else:
            print("Sample {0} with flowcell_id {1} not in database.".format(sample_id, flowcell_id))


def get_flowcell_id(flowcellsarg):
    flowcells = list(set(flowcellsarg))
    flowcells.sort()
    flowcell_id = "_".join(flowcells)
    return flowcell_id


def add_sample_to_db_and_return_refset_bam(bam, refset, print_refset_stdout=None):
    Session = connect_database()
    sample_id = get_sample_id(bam)
    flowcell_id = get_flowcelid_bam(bam)
    with Session() as session:
        if not session.query(Sample).filter(Sample.sample == sample_id).filter(Sample.flowcell == flowcell_id).all():
            sample = Sample(sample=sample_id, flowcell=flowcell_id, refset=refset)
            session.add(sample)
            session.commit()
        refset_db = (
            session.query(Sample).filter(Sample.sample == sample_id).filter(Sample.flowcell == flowcell_id).one().refset
        )
        if print_refset_stdout:
            print(refset_db)
        return refset_db


def return_refset_bam(bam):
    Session = connect_database()
    sample_id = get_sample_id(bam)
    flowcell_id = get_flowcelid_bam(bam)
    with Session() as session:
        if session.query(Sample).filter(Sample.sample == sample_id).filter(Sample.flowcell == flowcell_id).all():
            return session.query(Sample).filter(Sample.sample == sample_id).filter(Sample.flowcell == flowcell_id).one().refset


def get_sample_id(bam):
    with pysam.AlignmentFile(bam, "rb") as workfile:
        sampleid = []
        for readgroup in workfile.header['RG']:
            sampleid.append(readgroup['SM'])
        sampleid = list(set(sampleid))
        sampleid = "_".join(sampleid)
    return sampleid
