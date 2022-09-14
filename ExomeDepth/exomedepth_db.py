#! /usr/bin/env python3

import argparse
import database.functions

import settings


def print_message(status, sample_id=None, flowcell_id=None, refset=None):
    if status == "not_detected_in_db":
        print("## Sample {0} with flowcell_id {1} not in refset database".format(sample_id, flowcell_id))
    elif status == "added_to_db":
        print("## Sample {0} with flowcell_id {1} and with refset {2} added to database".format(
            sample_id, flowcell_id, refset
        ))
    elif status == "present_in_db":
        print("## Sample {0} with flowcell_id {1} and with refset {2} already in database".format(
            sample_id, flowcell_id, refset
        ))
    elif status == "updated_db":
        print("## Changed refset of sample {0} with flowcell_id {1} to refset {2}".format(
           sample_id, flowcell_id, refset
        ))
    elif status == "deleted_in_db":
        print("## Deleted sample {0} with flowcell_id {1}.".format(sample_id, flowcell_id))


def print_added(flowcell_id, sample_id, refset, added):
    if added:
        print_message("added_to_db", sample_id=sample_id, flowcell_id=flowcell_id, refset=refset)
    else:
        print_message("present_in_db", sample_id=sample_id, flowcell_id=flowcell_id, refset=refset)


def add_sample_to_db(args):
    flowcell_id, sample_id, refset_db, added = database.functions.add_sample_to_db(
        args.flowcell_id, args.sample_id, args.refset
    )
    if args.print_message_stdout:
        print_added(flowcell_id, sample_id, refset_db, added)
    if args.print_refset_stdout:
        print(refset_db)


def add_sample_to_db_and_return_refset_bam(args):
    flowcell_id, sample_id, refset_db, added = database.functions.add_sample_to_db_and_return_refset_bam(
        args.bam, args.refset
    )
    if args.print_message_stdout:
        print_added(flowcell_id, sample_id, refset_db, added)
    if args.print_refset_stdout:
        print(refset_db)


def change_refset_in_db(args):
    flowcell_id, sample_id, refset, updated = database.functions.change_refset_in_db(
        args.flowcell_id, args.sample_id, args.refset
    )
    if updated:
        print_message("updated_db", sample_id=sample_id, flowcell_id=flowcell_id, refset=refset)
    else:
        print_message("not_detected_in_db", sample_id=sample_id, flowcell_id=flowcell_id)


def query_refset(args):
    refset_db, flowcell_id = database.functions.query_refset(args.flowcell_id, args.sample_id)
    if refset_db:
        print(refset_db)
    else:
        print_message("not_detected_in_db", sample_id=args.sample_id, flowcell_id=flowcell_id)


def query_refset_bam(args):
    refset_db, flowcell_id, sample_id = database.functions.query_refset_bam(args.bam)
    if refset_db:
        print(refset_db)
    else:
        print_message("not_detected_in_db", sample_id=sample_id, flowcell_id=flowcell_id)


def delete_sample_db(args):
    removed, flowcell_id = database.functions.delete_sample_db(args.flowcell_id, args.sample_id)
    if removed:
        print_message("deleted_in_db", sample_id=args.sample_id, flowcell_id=flowcell_id)
    else:
        print_message("not_detected_in_db", sample_id=args.sample_id, flowcell_id=flowcell_id)


def print_all_samples(args):
    sample_list = database.functions.return_all_samples()
    print("Name\tFlowcell\tRefset")
    for line in sample_list:
        print(line)


def fill_database(args):
    folders = database.functions.get_folders(args.path)
    conflicts = {"warning": {}, "present": {}}
    for folder in folders:
        qc_file, bam_files, warning = database.functions.get_qc_bam_files("{}/{}".format(args.path, folder))
        if warning:
            continue
        sample_refset = database.functions.parse_refset_qc_file(qc_file)
        conflicts = database.functions.add_database_bam(bam_files, sample_refset, conflicts)

    args.conflict_file.write("Conflict\tSample\tFlowcellID\tRefsetDB\tRefsetSample\tBAM\n")
    for conflict in conflicts:
        for sample in conflicts[conflict]:
            args.conflict_file.write("{}\t{}\t{}\n".format(conflict, sample, "\t".join(conflicts[conflict][sample])))


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
    parser_add.add_argument('--print_message_stdout', action='store_true', help='print message if added to db or not')
    parser_add.add_argument('--print_refset_stdout', action='store_true', help='print refset in stdout')
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
    parser_add_return_bam.add_argument(
        '--print_message_stdout', action='store_true', help='print message if added to db or not'
    )
    parser_add_return_bam.add_argument('--print_refset_stdout', action='store_true', help='print refset in stdout]')
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

    """ Arguments fill database based on cohort sampels"""
    parser_fill_database = subparser.add_parser(
        'fill_database',
        help='Fill database based on folder path and BAM files'
    )
    parser_fill_database.add_argument('path', help='path of folder to be included')
    parser_fill_database.add_argument(
        '--conflict_file',
        default="conflicts.txt",
        type=argparse.FileType('w', encoding='UTF-8'),
        help='output file name containing conflicts that must be resolved manually (default = conflicts.txt)'
    )

    parser_fill_database.set_defaults(func=fill_database)

    args = parser.parse_args()
    args.func(args)
