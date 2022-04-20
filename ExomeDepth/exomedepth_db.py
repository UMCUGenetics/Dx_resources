#! /usr/bin/env python3

import argparse
from database.database import connect_database
from database.models import Sample

import pysam
import settings


def add_sample_to_db(args):
    Session = connect_database()
    flowcell_id = get_flowcell_id(args.flowcell_id)
    with Session() as session:
        if not session.query(Sample).filter(Sample.sample == args.sample_id).filter(Sample.flowcell == flowcell_id).all():
            sample = Sample(sample=args.sample_id, flowcell=flowcell_id, refset=args.refset)
            session.add(sample)
            session.commit()
            print("## Sample {0} added to database with flowcell {1} and with refset {2}".format(
                args.sample_id, flowcell_id, args.refset)
            )


def change_refset_in_db(args):
    Session = connect_database()
    with Session() as session:
        sample_update = (
            session.query(Sample)
            .filter(Sample.sample == args.sample_id)
            .filter(Sample.flowcell == args.flowcell_id)
            .one()
        )
        sample_update.refset = args.refset
        session.add(sample_update)
        session.commit()
        print("## Changed refset of sample {0} with flowcell {1} to refset {2}".format(
            args.sample_id, args.flowcell_id, args.refset)
        )


def print_all_samples(args):
    Session = connect_database()
    with Session() as session:
        print("Name\tFlowcell\tRefset\tFamilyID")
        for item in session.query(Sample):
            print("{0}\t{1}\t{2}".format(item.sample, item.flowcell, item.refset))


def print_refset(flowcell_id, sample_id):
    Session = connect_database()
    with Session() as session:
        print(
            session.query(Sample)
            .filter(Sample.sample == sample_id)
            .filter(Sample.flowcell == flowcell_id)
            .one()
            .refset
        )


def query_refset(args):
    flowcell_id = get_flowcell_id(args.flowcell_id)
    print_refset(flowcell_id, args.sample_id)


def query_refset_bam(args):
    flowcell_id = get_flowcelid_bam(args.bam)
    sample_id = get_sample_id(args.bam)
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


def delete_sample_db(args):
    Session = connect_database()
    flowcell_id = get_flowcell_id(args.flowcell_id)
    with Session() as session:
        samples = session.query(Sample).filter(Sample.sample == args.sample_id).filter(Sample.flowcell == flowcell_id)
        if samples:
            samples.delete()
            session.commit()
            print("Deleted sample {0} with flowcell {1}.".format(args.sample_id, flowcell_id))
        else:
            print("Sample {0} with flowcell {1} not in database.".format(args.sample_id, flowcell_id))


def get_flowcell_id(flowcellsarg):
    flowcells = list(set(flowcellsarg))
    flowcells.sort()
    flowcell_id = "_".join(flowcells)
    return flowcell_id


def add_sample_to_db_and_return_refset_bam(args):
    Session = connect_database()
    sample_id = get_sample_id(args.bam)
    flowcell_id = get_flowcelid_bam(args.bam)
    with Session() as session:
        if not session.query(Sample).filter(Sample.sample == sample_id).filter(Sample.flowcell == flowcell_id).all():
            sample = Sample(sample=sample_id, flowcell=flowcell_id, refset=args.refset)
            session.add(sample)
            session.commit()
        refset_db = (
            session.query(Sample).filter(Sample.sample == sample_id).filter(Sample.flowcell == flowcell_id).one().refset
        )
        if args.print_refset_stdout:
            print(refset_db)
        return refset_db


def get_sample_id(bam):
    with pysam.AlignmentFile(bam, "rb") as workfile:
        sampleid = []
        for readgroup in workfile.header['RG']:
            sampleid.append(readgroup['SM'])
        sampleid = list(set(sampleid))
        sampleid = "_".join(sampleid)
    return sampleid


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparser = parser.add_subparsers()

    """ Arguments add sample to database"""
    parser_add = subparser.add_parser('add_sample', help='add sample to database based on flowcell barcode')
    parser_add.add_argument('sample_id', help='sample id')
    parser_add.add_argument('flowcell_id', nargs='+', help='flowcell barcode')
    parser_add.add_argument(
        '--refset', default=settings.refset,
        help='exomedepth reference set ID [default = settings.refset]'
    )
    parser_add.set_defaults(func=add_sample_to_db)

    """ Arguments add sample to database based on BAM file and return refset"""
    parser_add_return_bam = subparser.add_parser(
        'add_sample_return_refset_bam',
        help='add sample to db based on BAM file and return refset'
    )
    parser_add_return_bam.add_argument('bam', help='full path to BAM file')
    parser_add_return_bam.add_argument(
        '--refset', default=settings.refset,
        help='exomedepth reference set ID [default = settings.refset]'
    )
    parser_add_return_bam.add_argument('--print_refset_stdout', default=True, help='print refset in stdout [default = True]')
    parser_add_return_bam.set_defaults(func=add_sample_to_db_and_return_refset_bam)

    """ Arguments query refset"""
    parser_query_refset = subparser.add_parser(
        'query_refset',
        help='return refset based on sampleid and flowcellid'
    )
    parser_query_refset.add_argument('sample_id', help='sample id')
    parser_query_refset.add_argument('flowcell_id', nargs='+', help='flowcell barcode')
    parser_query_refset.set_defaults(func=query_refset)

    """ Arguments query refset based on BAM file"""
    parser_query_refset_bam = subparser.add_parser(
        'query_refset_bam',
        help='return refset of sample based on BAM file'
    )
    parser_query_refset_bam.add_argument('bam', help='full path to BAM file')
    parser_query_refset_bam.set_defaults(func=query_refset_bam)

    """ Arguments change refset in database for sample"""
    parser_change = subparser.add_parser(
        'change_refset', help='change refset of a sample based on sampleid and flowcellid'
    )
    parser_change.add_argument('sample_id', help='sample id')
    parser_change.add_argument('flowcell_id', help='flowcell barcode (as in database)')
    parser_change.add_argument('refset', help='new refset ID')
    parser_change.set_defaults(func=change_refset_in_db)

    """ Arguments to print all entries of database """
    parser_allsamples = subparser.add_parser(
        'print_all_samples',
        help='Print full database table (tab-separated)'
    )
    parser_allsamples.set_defaults(func=print_all_samples)

    """ Arguments delete sample in database"""
    parser_delete_sample = subparser.add_parser(
        'delete_sample',
        help='Delete samples from database (sample + flowcell)'
    )
    parser_delete_sample.add_argument('sample_id', help='sample id')
    parser_delete_sample.add_argument('flowcell_id', nargs='+', help='flowcell barcode')
    parser_delete_sample.set_defaults(func=delete_sample_db)

    args = parser.parse_args()
    args.func(args)
