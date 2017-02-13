#! /usr/bin/env python
import sys, os, re, commands, getopt
from optparse import OptionParser
from optparse import OptionGroup
from os.path import join, isdir, isfile, split
from os import system, walk, listdir, chdir, rename, stat
from email import encoders
from email.mime.multipart import MIMEMultipart
from email.mime.base import MIMEBase
from email.mime.text import MIMEText
import smtplib
import mimetypes
import automatic_settings as s
import datetime, time
sys.dont_write_bytecode = True ## do not make .pyc of module(s)

## Author: M.G. Elferink (code borrowed from R. Ernst :) )
## Date: 25-01-2017
## Purpose: this script will be used in a cronjob to automatically fire IAP analysis for illumina NGS data.
## Usage: ./automatic_IAP.py

################
## Functions  ##
################

def send_email(sender, receivers, subject, text, attachment=None):
	mail = MIMEMultipart()
	mail['Subject'] = subject
	mail['From'] = sender
	mail['To'] = ';'.join(receivers)
	if attachment:
	       	filename = attachment.split('/')[-1]
	       	fp = open(attachment, 'rb')
	       	ctype, encoding = mimetypes.guess_type(attachment)
		if ctype is None or encoding is not None:
        		# No guess could be made, or the file is encoded (compressed), so
        		# use a generic bag-of-bits type.
        		ctype = 'application/octet-stream'
		maintype, subtype = ctype.split('/', 1)
		msg = MIMEBase(maintype, subtype)
		msg.set_payload(fp.read())
		fp.close()
		# Encode the payload using Base64
		encoders.encode_base64(msg)
		msg.add_header('Content-Disposition', 'attachment', filename=filename)
		mail.attach(msg)
	msg = MIMEText(text,'html')
	mail.attach(msg)
	m = smtplib.SMTP('smtp-open.umcutrecht.nl')
	m.sendmail(sender, receivers, mail.as_string())
	m.quit()
	
def get_sub_dir_path(directory):
	"""Return absolute path to all directories in input directory."""
	return filter(isdir,[join(directory, f) for f in listdir(directory)])

def getfile_insensitive(path):
    directory, filename = os.path.split(path)
    directory, filename = (directory or '.'), filename.lower()
    for f in os.listdir(directory):
        newpath = join(directory, f)
        if isfile(newpath) and f.lower() == filename:
            return newpath

def isfile_insensitive(path):
    return getfile_insensitive(path) is not None

def check_ped_file(run_dir,ped_file):
	correct=[]
	if(isfile_insensitive(ped_file)): 
		ped_file= getfile_insensitive(ped_file)
		lines=open(str(ped_file),"r").readlines()
		header_count=0
		line_count=0
		item_count=0
		for line in lines:
			if "Family" in line:	# check if header is still present in file
				header_count+=1
			if len(line.split()) != 6:	# check if all lines have the 6 required fields
				line_count+=1
			for item in line:	# check if all cells are filled
				if item == "":
					item_count+=1
		if header_count>0:
			correct+=["header"]
		if line_count > 0:
			correct+=["line_incorrect"]
		if item_count >0:
			correct+=["cell_empty"]
	else:
		correct=["PED_file_missing"]
		
	if len(correct)==0:
		return True,correct
	else:
		#make_mail("ped_file",run_dir,None,correct)
		return False, correct
		
def check_timestamp(run_dir):
	folderCreation = os.path.getctime(run_dir)
	date_now = time.time()
	seconds_ago = date_now - 60*60*24*14 # Number of seconds in 14 days
        if folderCreation < seconds_ago:
        	return True
       	else:
		return False

def check_raw_run(sequencer_dirs,processed,ped_folder):
	""" Check all folders in raw_data and determine if transfer is complete, pedigree is correct """
	run_dic={}
	error_dic={}
	for project in s.projects:
		for sequence_dir in sequencer_dirs:
			for run_dir in get_sub_dir_path(sequence_dir):
				try:
					project_folder=get_sub_dir_path(str(run_dir)+"/Data/Intensities/BaseCalls/")
					for folder in project_folder:
						if str(project) in str(folder.split("/")[-1]).upper():
							project_run=str(folder.split("/")[-1])
							transfer_done = '{}/{}'.format(run_dir, 'TransferDone.txt')
							ped_file='{}/{}'.format(run_dir,str(run_dir.split("/")[-1])+'.ped')
							ped_file_target='{}/{}'.format(str(ped_folder),str(run_dir.split("/")[-1])+'_'+str(project_run)+'.ped')
							run_id = '{}/{}'.format(processed,(run_dir.split("/")[-1]+"_"+str(project_run)))
							analysis_complete = '{}/{}'.format(run_id, 'processing.done')
							ped_file_format=check_ped_file(run_dir,ped_file)
							if time_check=="on":
								old=check_timestamp(run_dir)
							elif time_check=="off":
								old=False
	
							if(isfile(transfer_done) and ped_file_format[0] is True and not isfile(analysis_complete) and old is False): ## ignore runs with .done in processed data.
								run_dic[run_dir]=run_id,project_run
							if (not isfile(ped_file_target)):
								os.system("cp "+str(ped_file)+" "+str(ped_file_target))
								os.system("chmod 775 "+str(ped_file_target))
							elif(isfile(transfer_done) and ped_file_format[0] == False and not isfile(analysis_complete)):
								printline="PED file is incorrect: "
	        	                	        	for item in ped_file_format[1]:
       	        	                 		        	printline+=item + "\t"
        	        	               			printline +="Run="+str(run_id)
                                	      			error_dic[run_id]=printline
								print "PED_file incorrect or missing "+str(sequence_dir)
						###### DANGER, if processed run is deleted before raw data, the analysis is performed again! ######
						###### Write additional file in run completed folder (combined with specs?)
				except:
					pass
	return run_dic,error_dic

def html_table(text,data):
	column=len(data.split("\n")[0].split())
	length=0
	for item in data.split("\n"):
		for cell in item.split():
			if len(item) > length:
				length=len(item)
	width =(length * (column+15))
	text = '{}\n<table width={}>'.format(text,width)
	text= '{}\n<tr>'.format(text)
	for item in data.split("\n")[0].split():
        	splititem=item.split()
        	for cell in splititem:
                	text=('{}\n\n<th>'+str(cell.replace("_dedup",""))+'</th>').format(text)
	text= '{}\n</tr>'.format(text)
      	for item in data.split("\n")[1:]:
        	text= '{}\n<tr>'.format(text)
        	splititem=item.split()
        	for cell in splititem:
                	text=('{}\n\n<td><center>'+str(cell)+'</center></td>').format(text)
               	text= '{}\n</tr>'.format(text)
     	text = '{}\n</table>'.format(text)
	return text

def make_mail(state,run_id,qc=None,correct=None):
        """ Send failed/completed mail """
	if state == "started":
                subject = 'IAP started: IAP analysis for run {} started!'.format(run_id.split("/")[-1])
		text = "<html><body><p>"+'IAP analysis for {} has started'.format(run_id)+"</p></body></html>"
                email_to = [s.email_from]
        if state == "failed":
                subject = 'ERROR: IAP analysis {} failed 3 times!'.format(run_id.split("/")[-1])
		text = "<html><body><p><font color=\"red\">"+'IAP analysis of {} failed 3 times'.format(run_id)+"</font></p></body></html>"
                email_to = s.failed_mail

        if state == "completed":
		subject = 'COMPLETE: IAP analysis {} completed succesfully!'.format(run_id.split("/")[-1])
                text = "<html><body><p>"+'IAP analysis of {} completed succesfully.'.format(run_id)
              	text = '{}<br><br><u>QCStats:</u>'.format(text)
                text=html_table(text,qc[0])
                text = '{}<u>GENDER_check:</u>'.format(text)
                text=html_table(text,qc[1])
		text = '{}<u>KINSHIP_analysis:</u>'.format(text)
                text=html_table(text,qc[2])
                text = '{}</p></body></html>'.format(text)
                email_to = s.finished_mail	

	if state == "warning":
        	subject = 'WARNING: IAP analysis {} completed with WARNINGS!'.format(run_id.split("/")[-1])
		text = "<html><body><p>"+\
                        'IAP analysis of {} completed with WARNINGS.<br><br>\
                        <b><font color="red">DATA WILL NOT BE TRANSFERED TO BGARRAY!</font></b><br><br>\
                        <u>WARNINGS:</u>'.format(run_id)+"<br>"

		for item in qc[3]:
			text = '{}<br>{}\n'.format(text,item.rstrip())

		text = '{}<br><br><u>QCStats:</u>'.format(text)
		text=html_table(text,qc[0])
		text = '{}<u>GENDER_check:</u>'.format(text)
		text=html_table(text,qc[1])
		text = '{}<u>KINSHIP_analysis:</u>'.format(text)
                text=html_table(text,qc[2])
		text = '{}</p></body></html>'.format(text)
		email_to = s.finished_mail

	""" # excluded as mail will be send infinity
	if state == "ped_file":
		subject = 'PED file missing/incomplete for run {}'.format(run_id.split("/")[-1])
		text = 'PED file missing/incomplete for run {}.\n IAP not started!\n\nReason:\n'.format(run_id)
		for item in correct:
			 text = '{} {}'.format(text,item)
		email_to = s.finished_mail
	"""	

	send_email(s.email_from, email_to, subject, text)

def submit_IAP(IAP,INI,output,input,email_from,project):
	os.system("perl "+str(IAP)+"/illumina_createConfig.pl -ip "+str(INI)+ " -o " + str(output) + " -f " + str(input) + \
                  "/Data/Intensities/BaseCalls/"+ str(project) + " -m "+ str(email_from) + " -run")
	make_mail("started",output)	

def clean_up(run_id):
	os.chdir(run_id)
	os.system("python "+str(s.clean_up_folder))
	
def extract_qc(run_id,run):
        #add information of Conversion Report to trend analysis file. Q30, Cluster, PF, PF%, mismatch. Per sample/per run?

	"""QCSTATS """
	qcstats = '{}/{}'.format(run_id, '/QCStats/HSMetrics_summary.transposed.txt')
	qc_lines=open(str(qcstats),"r").readlines()
	qc=["Sample\t",qc_lines[0][1:]] # add sample name = first line in output
	try:
		len(warning)
	except:
		warning=[]
	for line in qc_lines:
		try:
			splitline=line.split()
			s.qc_metric[splitline[0]]
			qc+=[line]
			x=0
			while x< len(splitline):
				try:
					s.qc_criteria[splitline[0]] # check if first column (header) part of quality criteria
					try:
						keys= s.qc_criteria[splitline[0]].keys() #if dict in dict, such as PARENT and CHILD coverage
						sample_role=qc_lines[0].split()[x][7:8]
						if sample_role == "C":
							if float(s.qc_criteria[splitline[0]]["CHILD"])> float(splitline[x+1]):
								warning+=[" ".join(["WARNING coverage for sample",qc[1].split()[x], 
                                                                          "is below threshold "+str(splitline[0])+ 
                                                                          " (<"+str(float(s.qc_criteria[splitline[0]]["CHILD"]))+")!\n"])]
						if sample_role == "P":
							if float(s.qc_criteria[splitline[0]]["PARENT"])> float(splitline[x+1]):
								warning+=[" ".join(["WARNING coverage for sample",qc[1].split()[x],
                                                                          "is below threshold "+str(splitline[0])+  
                                                                          " (<"+str(float(s.qc_criteria[splitline[0]]["PARENT"]))+")!\n"])]
					except:
						if float(s.qc_criteria[splitline[0]]) > float(splitline[x+1]):
							warning+=[" ".join(["WARNING coverage for sample",qc[1].split()[x],"is below threshold "
                                                                  +str(splitline[0])+" (<"+str(s.qc_criteria[splitline[0]])+")!\n"])]
				except:	
					pass
				x+=1

		except:
			pass
	qc= "".join(qc)

	"""GENDER CHECK"""
	gender_file = '{}/{}'.format(run_id, '/gender_check.out')
	gender_lines=open(str(gender_file),"r").readlines()
	gender=[]
	for line in gender_lines:
		splitline=line.split()
		gender+=["\t".join(splitline),"\n"]
		if "FID" not in line:
			if splitline[4] == "OK":
				pass
			else:
				warning+=["WARNING gender of sample "+str(splitline[1])+" is incorrect"]

	gender="".join(gender)

	"""KINSHIP """
	kinship_file = '{}/{}.kinship'.format(run_id, run_id.split("/")[-1])	
	kinship_lines=open(str(kinship_file),"r").readlines()
	kinship=[]
	for line in kinship_lines:
                splitline=line.split()
		kinship+=[splitline[1]+"\t"+splitline[3]+"\t"+splitline[4]+"\t"+splitline[5]+"\t"+splitline[6]+"\t"+splitline[7]+"\n"] 
		printline=[splitline[1],splitline[3],splitline[7]]
		warning+=check_kinship(printline)
	kinship="".join(kinship)

	"""
	# determine number of lanes?
	# open conversion mail of lanes, and parse info
	conversion= commands.getoutput("find "+str(run)+"/Data/Intensities/BaseCalls/Reports/html/*/all/all/all/lane.html")
	con_lines=open(str(conversion),"r").readlines()
	for line in con_lines:	
		print line
		pass
	"""
	return qc,gender,kinship,warning

def touch_link(run_id):
        """ Make softlink on Upload folder """
	os.system("ln -sd " + str(run_id) + " "+str(s.rsync_folder))
        # check transfer is with touched Transfer.done? (only accesible vai hpct's!)

def check_kinship(kinship):
	kin_warn=[]
	if "ID1" not in kinship[0]:
		sample1_fam=kinship[0][0:7]
		sample1_gen=kinship[0][7:8]
		sample2_fam=kinship[1][0:7]
		sample2_gen=kinship[1][7:8]
		"""Check relations with kinship values"""
		if sample1_fam == sample2_fam and sample1_gen=="C" and sample2_gen=="P" or sample1_gen=="P" and sample2_gen=="C":
			if float(kinship[2])>0.177 and float(kinship[2])< 0.354:
				pass # kinship is correct
			else:
				kin_warn+=["PARENT-CHILD relation is not correct "+" ".join(kinship)]
				
		elif sample1_fam == sample2_fam and  sample1_gen=="C" and sample2_gen=="C":
			if float(kinship[2])>0.177 and float(kinship[2])< 0.354:
        			pass # kinship is correct
			else:
				kin_warn+=["SIBLING-SIBLING relation is not correct "+" ".join(kinship)]

		elif sample1_fam == sample2_fam and sample1_gen=="P" and sample2_gen=="P":
			if float(kinship[2])>0.177 and float(kinship[2])< 0.354:
				kin_warn+=["PARENT-PARENT relation is not correct "+" ".join(kinship)]
			else:
				pass  # kinship is correct
				# note: kinship between 0.10 and 0.177 is not mentioned here!

		elif sample1_fam is not sample2_fam:
			if float(kinship[2])>0.177 and float(kinship[2])< 0.354:
				kin_warn+=["UNRELATED relation is not correct "+" ".join(kinship)]
			else:
                               	pass  # kinship is correct
				# note: kinship between 0.10 and 0.177 is not mentioned here!
	return kin_warn

def run_IAP(runs,IAP,INI):
	""" Check state of the run """
	state=[]
	for run in runs:
		run_id = runs[run][0]
		project= runs[run][1]
		running = '{}/{}'.format(run_id, 'processing.running')
                error = '{}/{}'.format(run_id, 'processing.error')
		error_mail = '{}/{}'.format(run_id, 'processing.error.mail')
		done ='{}/{}'.format(run_id, 'processing.done')
		if(isfile(running)) and not isfile(error_mail):
			""" check if run is completed, or failed (again) """	
			log_file='{}/{}'.format(run_id,'/logs/PipelineCheck.log')
			if(isfile(log_file)):
				lines=open(str(log_file),"r").readlines()
				if "completed" in lines[-1]:
					clean_up(run_id)
					qc=extract_qc(run_id,run)
					if not qc[3] or overrule == "force":
						"""run completed without warnings """
						touch_link(run_id)
						make_mail("completed",run_id,qc,None)
						state+=["completed",run_id]
        	                                os.system("touch "+ str(done))
                	                        os.system("rm "+str(running))
					else:
						"""run completed with warnings, do not transfer, chmod 775 """
						make_mail("warning",run_id,qc,None)
						state+=["completed_warnings",run_id]
			                        os.system("touch "+ str(done))
                                        	os.system("rm "+str(running))
                                        	os.system("chmod 775 -R "+ str(run_id))
					if(isfile(error)):
                                                os.system("rm "+str(error))

				elif "failed" in lines[-1]:
					os.system("touch "+ str(error))
					submitfiles=commands.getoutput("ls "+str(run_id)+"/submit*").split()
					if len(submitfiles) > 2:
						"""if IAP has been runned 3 times, fail run"""
						os.system("rm "+ str(running))
					else:
						"""re-submit IAP"""
						os.system("rm "+ str(log_file))
						submit_IAP(IAP,INI,run_id,run,s.email_from,project)	
						state+=["failed_resubmit",run_id]

		elif(isfile(error) and not isfile(running)) and not isfile(error_mail):
			""" Send email that run has failed, and only once """
			make_mail("failed",run_id)
			print "Run has failed 3 times !!! "
			os.system("touch "+str(error_mail))
			state+=["failed_third_time",run_id]
			os.system("chmod 775 -R "+ str(run_id))

		elif (not isfile(running) and not isfile(error) and not isfile(error_mail)):
			""" start IAP and touch .running file """
			submit_IAP(IAP,INI,run_id,run,s.email_from,project)
			os.system("touch "+str(running))
			state+=["iap_started",run_id]
	return state

#################
## Main script ##
#################
if __name__ == "__main__":
        parser = OptionParser();
        group = OptionGroup(parser, "Main options")
	group.add_option("-i", default="/hpc/cog_bioinf/diagnostiek/development/DEV_Dx_INI/UMCU_DX_EXOME.ini", dest="INI", type='string',\
                         help="path to INI file. Default = /hpc/cog_bioinf/diagnostiek/development/DEV_Dx_INI/UMCU_DX_EXOME.ini")
	group.add_option("-d", default="/hpc/local/CentOS7/cog_bioinf/IAP_Dx/", dest="IAP",\
                         help="path to IAP folder. Default = /hpc/local/CentOS7/cog_bioinf/IAP_Dx/]")
	#group.add_option("-r", default="/hpc/cog_bioinf/diagnostiek/raw_data/", dest="raw_data",\
        #                 help="path to RAW data folder. Default = /hpc/cog_bioinf/diagnostiek/raw_data/]")
	#group.add_option("-o", default="/hpc/cog_bioinf/diagnostiek/processed/Exome/", dest="output",\
        #                 help="output folder for processed data. Default = /hpc/cog_bioinf/diagnostiek/processed/Exome/]")
	group.add_option("-o", default="/hpc/cog_bioinf/diagnostiek/users/Martin/Test_automatic/processed/", dest="processed",\
                         help="output folder for processed data. Default = warning test]")
	group.add_option("-p", default="/hpc/cog_bioinf/diagnostiek/processed/PED_FILES/", dest="ped_folder",\
                         help="IAP PED folder. Default = /hpc/cog_bioinf/diagnostiek/processed/PED_FILES/]")
	group.add_option("-f", dest="overrule",choices=['force'], help="transfer data to bgarray even if FAILED. Default= None")
	group.add_option("-t", default="on", dest="time_check",choices=['on','off'], help="remove time_stamp check of 14 days. Default= on")
	parser.add_option_group(group)
        (opt, args) = parser.parse_args()
	INI=str(opt.INI)
	IAP=str(opt.IAP)
	processed=str(opt.processed)
	ped_folder=str(opt.ped_folder)
	time_check=str(opt.time_check)
	#RAW=str(opt.raw_data)
	if opt.overrule:
		overrule=str(opt.overrule)
	lists=check_raw_run(s.sequencer_dirs,processed,ped_folder)
	print run_IAP(lists[0], IAP, INI)


