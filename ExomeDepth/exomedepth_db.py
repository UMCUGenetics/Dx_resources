#!venv/bin/python
"""Exomedepth db CLI."""

import sys
import argparse
from database import connect_database
from models import Sample
import pysam
import settings

def store_sample(Session, sample):
    with Session() as session:
        session.add(sample)
        session.commit()

def create_sample(sample_id, flowcell_id, refset):
    entry = Sample(
        sample=sample_id,
        flowcell=flowcell_id,
        refset=refset
    )
    return entry

def add_sample_to_db(args):
    Session = connect_database()
    with Session() as session:
        if not session.query(Sample).filter(Sample.sample == args.sample_id).filter(Sample.flowcell == args.flowcell_id).all():
            entry = create_sample(args.sample_id, args.flowcell_id, args.refset)
            store_sample(Session, entry)
            print("## Sample {0} added to database with flowcell {1} and with refset {2}".format(args.sample_id, args.flowcell_id, args.refset))

def change_refset_in_db(args):
    Session = connect_database()
    with Session() as session:
        sample_update = session.query(Sample).filter(Sample.sample == args.sample_id).filter(Sample.flowcell == args.flowcell_id).one()
        sample_update.refset = args.refset
        session.add(sample_update)
        session.commit()
        print("## Changed refset of sample {0} with flowcell {1} to refset {2}".format(args.sample_id, args.flowcell_id, args.refset))

def print_all_samples(args):
    Session = connect_database()
    with Session() as session:
        print("Name\tFlowcell\tRefset")
        for item in session.query(Sample).all():
            print("{0}\t{1}\t{2}".format(item.sample, item.flowcell, item.refset))

def query_refset(args):
    Session = connect_database()
    with Session() as session:
        print(session.query(Sample).filter(Sample.sample == args.sample_id).filter(Sample.flowcell == args.flowcell_id).one().refset)

def get_flowcelid(bam):
    workfile = pysam.AlignmentFile(bam, "rb")
    readgroups = []
    for readgroup in workfile.header['RG']:
        if readgroup['PU'] not in readgroup:
            readgroups.append(readgroup['PU'])
    readgroups=list(set(readgroups))
    flowcell_id = "_".join(readgroups)
    return flowcell_id


def query_refset_bam(args):
    Session = connect_database()
    bam_file = []
    for bam in args.bam_id:
        if args.sample_id in bam:
            bam_file.append(bam)
    if len(bam_file) == 1:
        flowcell_id = get_flowcelid(bam_file[0])
        with Session() as session:
            print(session.query(Sample).filter(Sample.sample == args.sample_id).filter(Sample.flowcell == flowcell_id).one().refset)
    else:
        sys.exit("None or multiple BAMs files detected for same sample")

def delete_sample_db(args):
    Session = connect_database()
    with Session() as session:
        if session.query(Sample).filter(Sample.sample == args.sample_id).filter(Sample.flowcell == args.flowcell_id).all():
            session.query(Sample).filter(Sample.sample == args.sample_id).filter(Sample.flowcell == args.flowcell_id).delete()
            session.commit() 
            print("Deleted sample {0} with flowcell {1}.".format(args.sample_id, args.flowcell_id))
        else:
            print("Sample {0} with flowcell {1} not in database.".format(args.sample_id, args.flowcell_id))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparser = parser.add_subparsers()

    parser_add = subparser.add_parser('add_sample', help='add sample + flowcell to database')
    parser_add.add_argument('sample_id', help='sample id')
    parser_add.add_argument('flowcell_id', help='flowcell barcode')
    parser_add.add_argument('--refset', default = settings.refset, help='exomedepth reference set ID [default = settings.refser]')
    parser_add.set_defaults(func = add_sample_to_db)

    parser_query_refset = subparser.add_parser('query_refset', help='query refset of sample + flowcell')
    parser_query_refset.add_argument('sample_id', help='sample id')
    parser_query_refset.add_argument('flowcell_id', help='flowcell barcode')
    parser_query_refset.set_defaults(func = query_refset)

    parser_change = subparser.add_parser('change_refset', help='change refset of specific sample + flowcell')
    parser_change.add_argument('sample_id', help='sample id')
    parser_change.add_argument('flowcell_id', help='flowcell barcode')
    parser_change.add_argument('refset', help='new refset ID')
    parser_change.set_defaults(func = change_refset_in_db)

    parser_allsamples = subparser.add_parser('print_all_samples', help='Print database (tab-separated)')
    parser_allsamples.set_defaults(func = print_all_samples)

    parser_query_refset_bam = subparser.add_parser('query_refset_bam', help='query refset of sample based on BAM file')
    parser_query_refset_bam.add_argument('sample_id', help='sample id')
    parser_query_refset_bam.add_argument('bam_id', nargs='+', help='full path to BAM file(s), space seperated. Requires at least 1 BAM file')
    parser_query_refset_bam.set_defaults(func = query_refset_bam)

    parser_delete_sample = subparser.add_parser('delete_sample', help='Delete samples from database (sample + flowcell)')
    parser_delete_sample.add_argument('sample_id', help='sample id')
    parser_delete_sample.add_argument('flowcell_id', help='flowcell barcode')
    parser_delete_sample.set_defaults(func = delete_sample_db)

    args = parser.parse_args()
    args.func(args)
