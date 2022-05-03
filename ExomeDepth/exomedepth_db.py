#! /usr/bin/env python3

import argparse
import database.functions

import settings


def call_add_sample_to_db(args):
    database.functions.add_sample_to_db(args.flowcell_id, args.sample_id, args.refset, args.print_message)


def call_change_refset_in_db(args):
    database.functions.change_refset_in_db(args.flowcell_id, args.sample_id, args.refset)


def call_query_refset(args):
    database.functions.query_refset(args.flowcell_id, args.sample_id)


def call_query_refset_bam(args):
    database.functions.query_refset_bam(args.bam)


def call_delete_sample_db(args):
    database.functions.delete_sample_db(args.flowcell_id, args.sample_id)


def call_add_sample_to_db_and_return_refset_bam(args):
    database.functions.add_sample_to_db_and_return_refset_bam(args.bam, args.refset, args.print_refset_stdout)


def call_print_all_samples(args):
    database.functions.print_all_samples()


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
    parser_add.add_argument('--print_message', default=True, help='print message if added to db or not [default = True]')
    parser_add.set_defaults(func=call_add_sample_to_db)

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
    parser_add_return_bam.set_defaults(func=call_add_sample_to_db_and_return_refset_bam)

    """ Arguments query refset"""
    parser_query_refset = subparser.add_parser(
        'query_refset',
        help='return refset based on sampleid and flowcellid'
    )
    parser_query_refset.add_argument('sample_id', help='sample id')
    parser_query_refset.add_argument('flowcell_id', nargs='+', help='flowcell barcode')
    parser_query_refset.set_defaults(func=call_query_refset)

    """ Arguments query refset based on BAM file"""
    parser_query_refset_bam = subparser.add_parser(
        'query_refset_bam',
        help='return refset of sample based on BAM file'
    )
    parser_query_refset_bam.add_argument('bam', help='full path to BAM file')
    parser_query_refset_bam.set_defaults(func=call_query_refset_bam)

    """ Arguments change refset in database for sample"""
    parser_change = subparser.add_parser(
        'change_refset', help='change refset of a sample based on sampleid and flowcellid'
    )
    parser_change.add_argument('sample_id', help='sample id')
    parser_change.add_argument('flowcell_id', help='flowcell barcode (as in database)')
    parser_change.add_argument('refset', help='new refset ID')
    parser_change.set_defaults(func=call_change_refset_in_db)

    """ Arguments to print all entries of database """
    parser_allsamples = subparser.add_parser(
        'print_all_samples',
        help='Print full database table (tab-separated)'
    )
    parser_allsamples.set_defaults(func=call_print_all_samples)

    """ Arguments delete sample in database"""
    parser_delete_sample = subparser.add_parser(
        'delete_sample',
        help='Delete samples from database (sample + flowcell)'
    )
    parser_delete_sample.add_argument('sample_id', help='sample id')
    parser_delete_sample.add_argument('flowcell_id', nargs='+', help='flowcell barcode')
    parser_delete_sample.set_defaults(func=call_delete_sample_db)

    args = parser.parse_args()
    args.func(args)
