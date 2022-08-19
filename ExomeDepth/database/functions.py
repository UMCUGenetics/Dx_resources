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


def get_qc_bam_files(path):
    folders = set(glob.glob("{}*".format(path), recursive=True))
    print("Number of folders detected = {}".format(len(folders)))
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

    return qc_files, bam_files
 

def parse_refset_qc_files(qc_files):
    sample_refset = {}
    for qc_file in qc_files:
        print(qc_files[qc_file])
        with open(qc_files[qc_file], 'r') as refset_qc:
            for line in refset_qc.readlines():
                if "REFSET" in line:
                    splitline = line.rstrip().split(";")

                    #  Determine index position of "REFSET="
                    header = []
                    for item in splitline:
                        header.append(item.split("=")[0])
                    refset_index = header.index('REFSET')

                    #  Parse used reference set from QC file into dictionary
                    refset_sample = splitline[refset_index].replace("REFSET=", "")
                    sample_id = splitline[0] 

                    warning = "none"
                    if "WARNING" in line:
                        warning = ",".join(line.rstrip().split("\t")[1:])

                    if sample_id not in sample_refset:
                        sample_refset[sample_id] = {refset_sample : [warning]}
                    else:
                        if refset_sample not in sample_refset[sample_id]:
                            sample_refset[sample_id][refset_sample] = [warning]
                        else:
                            sample_refset[sample_id][refset_sample] += [warning]
                
    #for item in sample_refset:
    #    print("jrrp", item, sample_refset[item])

    return sample_refset


def add_database_bam(bam_files, sample_refset):
    conflicts = {"warning":{}, "present":{}, "multiple":{}}
    for bam in bam_files:
        print("Processing bam {}".format(bam))
        sample_id = get_sample_id(bam)
        flowcell_id = get_flowcell_id_bam(bam)

        if len(sample_refset[sample_id]) == 1:  # Only one refset known for sample in  all parsed samples
            refset = list(sample_refset[sample_id].keys())[0]
            if "none" not in sample_refset[sample_id][refset]:  # None needs to be present, otherwise warning only sample 
                refset_db, added = add_sample_flowcell_to_db(sample_id, flowcell_id, refset) 
                if not added and refset_db != refset:  #Already present in db with different refset
                    conflicts["present"][sample_id] = [flowcell_id, refset_db, refset, sample_refset[sample_id][refset]]
            else:  # check if warning only sample is already in db
                refset_db = parse_refset (flowcell_id, sample_id)
                if not refset_db:
                    conflicts["warning"][sample_id] = [flowcell_id, refset_db, refset, sample_refset[sample_id][refset]]
        else: # multiple refset detected
            non_warning_refsets = []
            for refset in sample_refset[sample_id]:
                if "none" in sample_refset[sample_id][refset]:  # Select refsets that have at least one analysis without a warning
                    non_warning_refsets.append(refset)
            if len(non_warning_refsets) == 1:  # Only one refset without warning, thus likely the correct one
                refset = non_warning_refsets[0]
                refset_db, added = add_sample_flowcell_to_db(sample_id, flowcell_id, refset)
                if not added and refset_db != refset:  #Already present in db with different refset
                    conflicts["warning"][sample_id] = [flowcell_id, refset_db, refset, sample_refset[sample_id][refset]]
            else:
                #conflicts = []
                #for refset in sample_refset[sample_id]:
                #    conflicts.append([refset , sample_refset[sample_id][refset]])  
                #not_added["multiple"][sample_id] = [flowcell_id, "not_present", conflicts]
                conflicts["multiple"][sample_id] = [flowcell_id, "not_present", non_warning_refsets, "none"]
                
    #for item in not_added:
    #    print(item, not_added[item])

    return conflicts


#def resolve_conflicts(not_added, sample_refset):
#    unresolved = []
#
#    for not_add in not_added:
#        print(not_add, sample_refset[not_add[0]])

        #if conflict[3][0] == conflict conflict[4]: 
        
        #for item in conflict[2]:
        #    print("jrrp",item)

        #if sample_refset[conflict[0]]:
        #    refset_list = []
        #    for item in sample_refset[conflict[0]]:
        #        ignore = False
        #        for warning in item[1]:
        #            if "DO_NOT_USE_MergeSample" in warning:
        #                ignore = True

        #else:
        #    print("Sample {} not detected in qc_files".format(conflict[0]))
        #    continue

    #for conflict in conflicts:
    #    if splitline[4] in qc_dic:
    #        refset_list = []
    #        for item in qc_dic[splitline[4]]:
    #            ignore = False
    #            ignore = 0
    #            for warning in item[1]:
    #                if "DO_NOT_USE_MergeSample" in warning:
    #                    ignore = True
    #                    ignore += 1
    #            if ignore == False:
    #                refset_list.append(item[0])
    
    #        if len(refset_list) == 1:  # Only 1 unique refset that is not a merge sample:  this must be the correct refset.
    #            action = (
    #                "source /diaggen/software/development/Dx_resources_v700/ExomeDepth/venv/bin/activate &&"
    #                "/diaggen/software/development/Dx_resources_v700/ExomeDepth/exomedepth_db.py change_refset {} {} {} "
    #            ).format(splitline[4], splitline[7], refset_list[0])
    #            os.system(action)
    #        else:  # Needs manual resolve
    #            print("Not resolved, check manually {} {} {}".format(splitline[4], splitline[7], refset_list))

    #    else: # In case sampleID is not detected in QC summary
    #        print("Sample not detected\t{}".format(splitline[4]))


    #return unresolved
