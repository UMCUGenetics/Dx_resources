#! /usr/bin/env python3

import pysam

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
