#!/usr/bin/env python
# -*- coding: utf-8 -*-

__version__ = 'v2.2.5'
print "###########分析流程版本号:",__version__

import os
import sys
import gzip
import re
import time
import argparse
import getpass
import subprocess

from argparse import RawTextHelpFormatter

from Get_config import *
from Basic_function import *
from QualityControl import QC
from Mapping import Mapping,Assembly,Quantification
from RNAVariants import *
from Differential import *

#####
pipe_dir = os.path.dirname(os.path.abspath(__file__))
root_path = os.path.dirname(os.path.abspath(pipe_dir))

basic_analysis = {
	'1':'quality_control', 
	'2':'#mapping', 
	'2.1':'mapping_hisat2',
	'2.2':'mapping_tophat2',
	'2.3':'mapping_star',
	'3':'#quantification',
	'3.1':'quant_rsem',
	'3.2':'quant_htseq',
	'3.3':'quant_cuffdiff',
	'4':'#diff_enrichment',
	'4.1':'diff_edgeR',
	'4.2':'diff_deseq',
	'4.3':'diff_cuffdiff',
	'5':'lncRNA_target',
	'6':'#snp_indel',
	'6.1':'GATK',
	'6.2':'samtools',
	'6.3':'VarDict',
	'7':'#fusion_gene',
	'7.1':'Arriba',
	'7.2':'starfusion',
	'8':'Altersplice',
	'9':'TTMV_screen',
	'10':'#circ_novel', 
	'10.1':'find_circ',
	'10.2':'ciri',
	'10.3':'both',
	'11':'circ_target',
	'12':'circ_diff',
	'13':'circ_enrich',
	'14':'circ_coding'
}
parser = argparse.ArgumentParser(description="Human RNA Pipeline",formatter_class=RawTextHelpFormatter)
parser.add_argument('--infile',help="a file include patient sample and library information",required=True)
parser.add_argument('--analydir',help="a path to store analysis results",required=True)
parser.add_argument('--seqtype',help="product type, default: RNAseq",default='RNAseq',choices=['RNAseq','lncRNA','circRNA']) # teng
parser.add_argument('--libtype',help="library type, default: fr-firststrand",default=None,choices=['fr-firststrand','fr-secondstrand','fr-unstranded'])
parser.add_argument('--analy-array',help="which steps of analysis will do. Default (1)\n	" + \
	"\n\t".join(['%s %s'%(e,basic_analysis[e]) for e in \
	['1','2','2.1','2.2','2.3','3','3.1','3.2','3.3','4','4.1','4.2','5','6','6.1','6.2','6.3','7','7.1','7.2','8','9','10','11','12','13','14']])) # teng
parser.add_argument('--newjob',help="the new job file, default: datetime.job",default=time.strftime('%m.%d',time.localtime(time.time()))+'.job')
parser.add_argument('--compare-group',help="Compare groups for differential analysis", default=None)
parser.add_argument('--startpoint',help="give a start analysis point,if need!",default=None)
parser.add_argument('--assembly',help="If assembly transcript or not.",action='store_true')
parser.add_argument('--assembly-method',help="Assembly method, default: cufflinks.",choices=['cufflinks','stringtie'],default='cufflinks')
parser.add_argument('--WGCNA',help="Whether or not a WGCNA analysis is carried out,Select this parameter to carry out this analysis",action='store_true')
parser.add_argument('--somaticfusion',help="Whether somatic fusion gene analysis should be carried out,Select this parameter to carry out this analysis",action='store_true')
parser.add_argument('--genome',help="Genome for analysis, default human_B37",default='human_B37',choices=['human_B37','human_B38','human_B38_NCCL_spikein','mm10','rat'])
parser.add_argument('--read-length',help="Read length, default 150 bp.",default=150,type=int)
parser.add_argument('--fqcheck',help="If fqcheck for rawdata.",action='store_true')
parser.add_argument('--singlecell',help="single cell sequence data or not",default=None,choices=['vazyme','NEB'])
parser.add_argument('--TR',help="a bed file for sequence region",default="%s/bed/HMRNA-210729_all_targets.sort.extend.merge.bed" % root_path)
parser.add_argument('--NP',help="Sample clinical information file",required=True)
parser.add_argument('-v','--version',action='version',version=__version__)

argv = parser.parse_args()
###############################################################################
TR = argv.TR
if not os.path.isfile(TR):
	raise IOError('%s is not exsit' % TR)
###
infile = argv.infile
if not os.path.isfile(infile):
	raise IOError('%s is not exsit' % infile)
###
analydir = argv.analydir
if not (os.path.exists(analydir) and os.path.isdir(analydir)):
	raise IOError('%s is not exsit' % analydir)
###
if os.path.isfile(argv.NP):
	NP_tmp = argv.NP
	if not NP_tmp.startswith("/"):
		NP = os.path.join(analydir,NP_tmp)
	else:
		NP = NP_tmp
else:
	raise IOError('%s does not exist' % argv.NP)
####################################################################################################################################################################################
disease_dict = get_disease_dict()

dic_disease={}
DNARNA_combine={}  ## ID: 姓名，医院，核收日期，项目名
with open(NP,'r') as f:
	for line in f:
		Line=line.strip().split("\t")
		## solid
		if re.search("^RT",Line[0]) and Line[22] in disease_dict:
			dic_disease[Line[0]] = disease_dict[Line[22]]	
			if re.match('实体瘤785基因检测全景版\(DNA\+RNA\)',Line[22]):
				DNARNA_combine[Line[0]] = [Line[4],Line[8],Line[13],Line[22]]
		## blood
		elif Line[13] in disease_dict:
			dic_disease[Line[0]] = disease_dict[Line[13]]
			#### DNA + RNA project
			if re.match(r'血液肿瘤基因筛查全面版\(DNA\+RNA\)',Line[13]) or re.match(r'血液肿瘤全面版\(DNA\+RNA\)-RNA',Line[13]):
				DNARNA_combine[Line[0]] = [Line[1],Line[6],Line[11],Line[13]]
			elif re.match(r'Ph-like ALL/LBL基因筛查\(DNA\+RNA\)-RNA',Line[13]):
				DNARNA_combine[Line[0]] = [Line[1],Line[6],Line[11],Line[13]]
			elif 'ALL相关基因突变筛查(DNA+RNA)' in Line[13] or "ALL相关基因变异筛查\(DNA\+RNA\)-RNA" in Line[13]:
				DNARNA_combine[Line[0]] = [Line[1],Line[6],Line[11],Line[13]]
			elif re.match(r'ALL相关基因变异筛查\(DNA\+RNA\)-RNA',Line[13]) or re.match(r'AML相关变异筛查70基因\(DNA\+RNA\)-RNA',Line[13]) or re.match(r'MPN相关变异筛查56基因\(DNA\+RNA\)-RNA',Line[13]) or re.match(r'血液肿瘤相关变异筛查302基因\(DNA\+RNA\)-RNA',Line[13]) or re.match(r'血液肿瘤相关变异筛查302基因\(DNA\+RNA\)\(肿瘤\+对照\)-RNA',Line[13]):
				DNARNA_combine[Line[0]] = [Line[1],Line[6],Line[11],Line[13]]
		else:
			print "\033[1;41;43m %s 检测项目名称不标准 \033[0m" % (Line[0])
#print dic_disease
#print "DNA+RNA:",DNARNA_combine
####################################################################################################################################################################################
times = time.strftime('%m.%d',time.localtime(time.time()))

startpoints = []
if argv.startpoint:
	startpoints = [str(x) for x in argv.startpoint.strip().split(',')]
## cluster queue
user = getpass.getuser()
default_queues = ['all.q','novo.q','cancer1.q','cancer2.q']
user_queues = set([each.split('@')[0] for each in os.popen('qselect -U %s'%user).read().strip().split('\n')])
useful_queues = [x for x in default_queues if x in user_queues]
queue_list = ' '.join(['-q %s'%q for q in useful_queues])
## analysis array
analy_array = [x for x in argv.analy_array.strip().split(',')]
analysis = {}
for each in analy_array:
	if basic_analysis[each].startswith('#'):
		each = each[1:]
		each += '.1'
	x = each.split('.')[0]
	if x not in analysis:
		analysis[x] = each
includes = set(analysis.keys())

## software for mapping, mutation, fusiongene, circRNA
map_soft = 'hisat2'
if analysis.get('2','2.1') == '2.2':
	map_soft = 'tophat2'
elif analysis.get('2','2.1') == '2.3':
	map_soft = 'star'
exp_soft = 'rsem'
if analysis.get('3','3.1') == '3.2':
	exp_soft = 'htseq'
elif analysis.get('3','3.1') == '3.3':
	exp_soft = 'cuffdiff'
elif analysis.get('3','3.1') == '3.4':
	exp_soft = 'stringtie'
diff_soft = 'edgeR'
if analysis.get('4','4.1') == '4.2':
	diff_soft = 'DESeq2'
if analysis.get('4','4.1') == '4.3':
	diff_soft = 'cuffdiff'

mut_soft = 'gatk'
if analysis.get('6','6.1') == '6.2':
	mut_soft = 'samtools'
elif analysis.get('6','6.1') == '6.3':
	mut_soft = 'VarDict'

if analysis.get('7','7.1') == '7.2':
	fusion_soft = 'starfusion'

circ_soft ='find_circ' # teng
if analysis.get('10','10.1') == '10.3':
	circ_soft = 'both'
elif analysis.get('10','10.1') == '10.2':
	circ_soft = 'ciri'

if (exp_soft == 'cuffdiff' and diff_soft == 'edgeR') or (exp_soft == 'cuffdiff' and diff_soft == 'DESeq2') or (exp_soft == 'htseq' and diff_soft == 'cuffdiff') or (exp_soft == 'rsem' and diff_soft == 'cuffdiff'):
	print 'exp_soft and diff_soft are in conflict'
	sys.exit(1)

for i in analy_array:
	if i not in basic_analysis:
		print str(i) + ' is not in analysis'
		sys.ext(1)

#give a relationship of analysis steps
relationship = {
	'1':[],
	'2':['1'],
	'2.1':['1'],
	'2.2':['1'],
	'2.3':['1'],
	'3':['1','2'],
	'3.1':['1','2'],
	'3.2':['1','2'],
	'3.3':['1','2'],
	'4':['1','2','3'],
	'4.1':['1','2','3'],
	'4.2':['1','2','3'],
	'4.3':['1','2','3'],
	'5':['1','2','3'],
	'6':['1','2'],
	'6.1':['1','2'],
	'6.2':['1','2'],
	'6.3':['1','2'],
	'7':['1','2'],
	'7.1':['1','2'],
	'7.2':['1','2'],
	'8':['1','2'],
	'9':['1'],
	'10':['1'], # teng
	'10.1':['1'],
	'10.2':['1'],
	'10.3':['1'],
	'11':['10'],
	'12':['10'],
	'13':['12'],
	'14':['10']
}

will_quit = False
for i in analy_array:
	s = set(relationship[i])
	if not s.issubset(set(analysis.keys())):
		print 'you need do %s before %s' %(', '.join(relationship[i]), i)
		will_quit = True
	else:pass
if will_quit:
	sys.exit(0)

## update qclist
qclist_file = os.path.join(analydir,'qc_list')
update_qclist(infile,qclist_file)

## parse qc_list and sample config
list_in_sample = {}
list_in_qc = {}
patientInfo = {}

##get singlecell index
index={}
if argv.singlecell:
	for line in safe_open(qclist_file,'r'):
		if line.startswith('#'):
			continue
		array = line.strip().split('\t')
		if array[2] not in index.keys():
			index[array[2]] = str(array[4]).strip().split(";")[0]
		else:
			continue
#print index
for line in safe_open(qclist_file,'r'):
	if line.startswith('#'):
		continue
	array = line.strip().split('\t')
	if not array[2] in list_in_qc:
		list_in_qc[array[2]] = {}
	if array[3] not in list_in_qc[array[2]]:
		list_in_qc[array[2]][array[3]] = {}
	tmpid = array[-1].strip('"') ##path_LaneID
	list_in_qc[array[2]][array[3]][tmpid] = array[0]  ## sampleID,libID,path_LaneID end_of path_LaneID
	array[6] = array[6].strip() ## N or T
	if array[1] not in patientInfo:  # judge wherher patientInfo is in patientInfo list ,if not
		patientInfo[array[1]] = {'N':'','T':[]} # T is a list?
	if array[6] == 'N':
		patientInfo[array[1]]['N'] = array[2]
	elif array[6] == 'T':
		if array[2] not in patientInfo[array[1]]['T']:
			patientInfo[array[1]]['T'].append(array[2])
#print patientInfo
dic_sample_library={}
t_type={}
tumor_type={}
sample_type={}
sample_db = {}

for line in safe_open(infile,'r'):
	if line.startswith('#') or line.strip() == '':
		continue
	array = line.strip().split('\t')
	sampleID = array[2]
	dic_sample_library[array[2]]=array[3]
	assert re.search(u'(\d+)',array[0])
	laneid = re.search(u'(\d+)',array[0]).group(1)
	tmp_fq1 = os.path.join(array[5],array[3],'%s_L%s_1.fq.gz'%(array[3],laneid))
	tmp_fq2 = os.path.join(array[5],array[3],'%s_L%s_2.fq.gz'%(array[3],laneid))
	if ',' in array[5] and array[5].endswith('.gz'):
		tmp_fq1,tmp_fq2 = array[5].split(',')

	if array[2] not in list_in_sample:
		list_in_sample[array[2]] = {}
	if array[3] not in list_in_sample[array[2]]:
		list_in_sample[array[2]][array[3]] = {}
	tmpid = '_'.join([array[2],array[3],array[5],array[0]])
	lane = list_in_qc[array[2]][array[3]][tmpid]
	list_in_sample[array[2]][array[3]][lane] = [tmp_fq1, tmp_fq2]
	if array[6] == 'T' and  array[7] not in t_type.keys():
		t_type[array[7]] = [sampleID]
	elif array[6] == 'T' and array[7] in t_type.keys() and sampleID not in t_type[array[7]]:
		t_type[array[7]].append(sampleID)
	if array[6] == 'T' and  sampleID not in tumor_type.keys():
		tumor_type[sampleID] = array[7]
	elif array[6] == 'T' and  sampleID in tumor_type.keys():
		continue
	if sampleID not in sample_type.keys():
		if array[4] == 'amplicon' or array[4] == 'Amplicon':
			print '''\033[1;41;43m note : %s 是amplicon建库\033[0m '''%sampleID
		sample_type[sampleID] = array[4]
	else:
		continue
	sample_db[sampleID] = array[7]

#print list_in_sample.keys()
## cut NEG 
list_in_sample_noNEG = {k: v for k, v in list_in_sample.items() if 'NEG' not in k}
#print list_in_sample_noNEG.keys()


print "tumor_type=%s"%tumor_type
if argv.WGCNA:
	if len(list_in_sample.keys())<15:
		print 'Less than 15 samples, WGCNA can not be carried out'
		sys.exit(1)
		
#libtype
if argv.libtype:
	libtype=argv.libtype
else:
	libtype_dic={'RNAseq':'fr-unstranded','lncRNA':'fr-firststrand','other':'fr-secondstrand','circRNA':'fr-firststrand'} # teng sometime circRNA :fr-unstranded
	libtype=libtype_dic[argv.seqtype]


#print list_in_sample

repeat_flag=False
all_compare_groups = {}

if argv.compare_group:
	for line in open(argv.compare_group):
		if line.startswith('#'):
			continue
		array = line.strip().split('\t')
		if len(array)!=5:
			print 'For compare.txt,all columns must be filled in'
		grp_name = '%s_vs_%s' % (array[1],array[3])
		if grp_name != array[0]:
			print 'Please confirm your comparison'
			sys.exit(1)
		if ','  in array[2] or ','  in array[4]:
			repeat_flag=True
		all_compare_groups[grp_name] = {'T':{'name':array[1], 'samples':array[2]}, \
			'N':{'name':array[3],'samples':array[4]}}
else:
	for eachP in patientInfo:
		for each in patientInfo[eachP]['T']:
			if patientInfo[eachP]['N'] == "":
				continue
			grp_name = '%s_vs_%s' % (each,patientInfo[eachP]['N'])
			all_compare_groups[grp_name] = {'T':{'name':each, 'samples':each}, \
				'N':{'name':patientInfo[eachP]['N'],'samples':patientInfo[eachP]['N']}}

if all_compare_groups == {}:
	print 'The comapre is empty'
this_compare_groups=all_compare_groups
all_groups_tmp=[]
all_groupname_tmp=[]
for eachgrp in this_compare_groups:
	all_groups_tmp.append(this_compare_groups[eachgrp]['T']['samples'].replace(',',':'))
	all_groups_tmp.append(this_compare_groups[eachgrp]['N']['samples'].replace(',',':'))
	all_groupname_tmp.append(this_compare_groups[eachgrp]['T']['name'])
	all_groupname_tmp.append(this_compare_groups[eachgrp]['N']['name'])

all_groups=[]
all_groupname=[]
for each in list(set(all_groupname_tmp)):
	all_groupname.append(each)
	all_groups.append(all_groups_tmp[all_groupname_tmp.index(each)])
### jobstatus file
logdir = os.path.join(analydir,'log')
create_dir(logdir)
## dir
rawdir = os.path.join(analydir,'RawData')
qcdir = os.path.join(analydir,'QC')
mapdir = os.path.join(analydir,'Mapping')
assemdir = os.path.join(analydir,'Assembly')
expdir = os.path.join(analydir,'Quantification')
mutdir = os.path.join(analydir,'Mutation')
pindeldir = os.path.join(mutdir,'Pindel')
conpairdir = os.path.join(mutdir,'Conpair')
zhushidir = os.path.join(mutdir,'Annotate')
arribadir = os.path.join(analydir,'FusionGene','Arriba')
starfusiondir = os.path.join(analydir,'FusionGene','Starfusion')
diffdir = os.path.join(analydir,'DiffEnrichment')
as_diffdir= os.path.join(analydir,'Altersplice')
TTMVscreendir = os.path.join(analydir,'TTMV_screen')
targetdir = os.path.join(analydir,'LncTarget')
enrich_targetdir = os.path.join(analydir,'LncTargetEnrich')
circdir =  os.path.join(analydir,'CircNovel') # teng
circcoding = os.path.join(analydir,'CircCoding')
circdiff = os.path.join(analydir,'CircDiff')
circbind = os.path.join(analydir,'CircBind')
resultdir = os.path.join(analydir,'Result')
reportdir = os.path.join(analydir,'Report')

new_jobs = {}
orders = []

job_points = {} ## for startpoints

databases = genome_files[argv.genome]
## enrich_spe: hsa, mmu
enrich_spe = 'hsa'
rRNA_spe = 'hsa'
if argv.genome == 'mm10':
	enrich_spe = 'mmu'
	rRNA_spe = 'mm10'
elif argv.genome == 'rat':
	rRNA_spe = 'rat'
	enrich_spe = 'rno'
if set(['1']).issubset(includes):
	print "QC ... \n"
	create_dir(rawdir)
	job_points['qc'] = []
	job_points['md5_raw'] = []
	job_points['md5_clean'] = []
	job_points['rm_rRNA'] = []
	job_points['pollution'] = []
	for eachsample in list_in_sample:
		myrawdir = os.path.join(rawdir,eachsample)
		myqcdir = os.path.join(qcdir,eachsample)
		create_dir(myrawdir)
		create_dir(myqcdir)
		aQC = QC(eachsample, rawdir, myqcdir, list_in_sample[eachsample],softwares, databases)
		## rm rRNA
		rmr_cmd = aQC.rm_rRNA(rRNA_rate=10,species=rRNA_spe)
		new_jobs['_'.join(['rm_rRNA',eachsample])] = \
			{'name' : '_'.join(['rm_rRNA',eachsample]),
			'status' : 'waiting',
			'sched' : '-V -cwd %s' % queue_list,
			'memory' : compute_resource['rm_rRNA']['m'],
			'cmd' : 'sh %s' % os.path.join(myrawdir,'_'.join(['rm_rRNA',eachsample])+'.sh')}
		safe_open(os.path.join(myrawdir,'_'.join(['rm_rRNA',eachsample])+'.sh'),'w').write(rmr_cmd)
		job_points['rm_rRNA'].append('_'.join(['rm_rRNA',eachsample]))
		## QC wwt
		if argv.singlecell:
			qc_cmd,order1 = aQC.qc(fqcheck=argv.fqcheck,singlecell=argv.singlecell,index=index[eachsample])
			add_items(orders,order1)
		else:
			qc_cmd,order1 = aQC.qc(fqcheck=argv.fqcheck,singlecell=argv.singlecell,index="")
			add_items(orders,order1)
		new_jobs['_'.join(['qc',eachsample])] = \
			{'name' : '_'.join(['qc',eachsample]),
			'status' : 'waiting',
			'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['qc']['t']),
			'memory' : compute_resource['qc']['m'],
			'cmd' : 'sh %s' % os.path.join(myqcdir,'_'.join(['qc',eachsample])+'.sh')}
		safe_open(os.path.join(myqcdir,'_'.join(['qc',eachsample])+'.sh'),'w').write(qc_cmd)
		job_points['qc'].append('_'.join(['qc',eachsample]))

		job_points['qc_summary']=[]
		qc_summary_cmd,order1 = aQC.qc_summary(dic_sample_library[eachsample],analydir)
		add_items(orders,order1)
		new_jobs['_'.join(['qc_summary',eachsample])] = \
			{'name' : '_'.join(['qc_summary',eachsample]),
			'status' : 'waiting',
			'sched' : '-V -cwd %s' % queue_list,
			'memory' : compute_resource['qc_summary']['m'],
			'cmd' : 'sh %s' % os.path.join(myqcdir,'_'.join(['qc_summary',eachsample])+'.sh')}
		safe_open(os.path.join(myqcdir,'_'.join(['qc_summary',eachsample])+'.sh'),'w').write(qc_summary_cmd)
		job_points['qc_summary'].append('_'.join(['qc_summary',eachsample]))

		## md5sum raw
		md5_cmd,order1 = aQC.md5_raw()
		add_items(orders,order1)
		new_jobs['_'.join(['md5_raw',eachsample])] = \
			{'name' : '_'.join(['md5_raw',eachsample]),
			'status' : 'waiting',
			'sched' : '-V -cwd %s' % queue_list,
			'memory' : compute_resource['md5_raw']['m'],
			'cmd' : 'sh %s' % os.path.join(myrawdir,'_'.join(['md5_raw',eachsample])+'.sh')}
		safe_open(os.path.join(myrawdir,'_'.join(['md5_raw',eachsample])+'.sh'),'w').write(md5_cmd)
		job_points['md5_raw'].append('_'.join(['md5_raw',eachsample]))
				## md5 clean
		md5_cmd,order1 = aQC.md5_clean()
		add_items(orders,order1)
		new_jobs['_'.join(['md5_clean',eachsample])] = \
			{'name' : '_'.join(['md5_clean',eachsample]),
			'status' : 'waiting',
			'sched' : '-V -cwd %s' % queue_list,
			'memory' : compute_resource['md5_clean']['m'],
			'cmd' : 'sh %s' % os.path.join(myqcdir,'_'.join(['md5_clean',eachsample])+'.sh')}
		safe_open(os.path.join(myqcdir,'_'.join(['md5_clean',eachsample])+'.sh'),'w').write(md5_cmd)
		job_points['md5_clean'].append('_'.join(['md5_clean',eachsample]))
			

if set(['1','2']).issubset(includes):
	print "Mapping ...\n"
	create_dir(mapdir)
	job_points['mapping_rnaseq'] = []
	job_points['cal_readcount'] = []
	job_points['mapping_summary'] = []
	if map_soft == 'star':
		qsub_memory='25G'
		qsub_p='6'
	else:
		qsub_memory='5G'
		qsub_p='4'
	for eachsample in list_in_sample:
		myqcdir = os.path.join(qcdir,eachsample)
		mymapdir = os.path.join(mapdir,eachsample)
		create_dir(mymapdir)
		aMapping = Mapping(eachsample,myqcdir,mymapdir,softwares,databases,TR)
		map_cmd,order1 = aMapping.mapping_rnaseq(soft=map_soft,libtype=libtype)
		add_items(orders,order1)
		new_jobs['_'.join(['mapping_rnaseq',eachsample])] = \
			{'name' : '_'.join(['mapping_rnaseq',eachsample]),
			'status' : 'waiting',
			'sched' : '-V -cwd %s -l p=%s' % (queue_list,qsub_p),
			'memory' : qsub_memory,
			'cmd' : 'sh '+os.path.join(mymapdir,'_'.join(['mapping_rnaseq',eachsample])+'.sh')}
		safe_open(os.path.join(mymapdir,'_'.join(['mapping_rnaseq',eachsample])+'.sh'),'w').write(map_cmd)
		job_points['mapping_rnaseq'].append('_'.join(['mapping_rnaseq',eachsample]))

		mypolldir=os.path.join(mymapdir,'pollution')
		create_dir(mypolldir)
		pol_cmd,order1 = aMapping.pollution()
		safe_open(os.path.join(mypolldir,'_'.join(['pollution',eachsample])+'.sh'),'w').write(pol_cmd)
		add_items(orders,order1)
		new_jobs['_'.join(['pollution',eachsample])] = \
			{'name' : '_'.join(['pollution',eachsample]),
			'status' : 'waiting',
 			'sched' : '-V -cwd %s' % queue_list,
			'memory' : compute_resource['pollution']['m'],
			'cmd' : 'sh %s' % os.path.join(mypolldir,'_'.join(['pollution',eachsample])+'.sh')}
		job_points['pollution'].append('_'.join(['pollution',eachsample]))

	
		## map cal_readcount
		stat_cmd,order1 = aMapping.cal_readcount(libtype)
		add_items(orders,order1)
		new_jobs['_'.join(['cal_readcount',eachsample])] = \
			{'name' : '_'.join(['cal_readcount',eachsample]),
			'status' : 'waiting',
			'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['cal_readcount']['t']),
			'memory' : compute_resource['cal_readcount']['m'],
			'cmd' : 'sh '+os.path.join(mymapdir,'_'.join(['cal_readcount',eachsample])+'.sh')}
		safe_open(os.path.join(mymapdir,'_'.join(['cal_readcount',eachsample])+'.sh'),'w').write(stat_cmd)
		job_points['cal_readcount'].append('_'.join(['cal_readcount',eachsample]))
	
		## map stat
		stat_cmd,order1 = aMapping.mapping_summary()
		add_items(orders,order1)
		new_jobs['_'.join(['mapping_summary',eachsample])] = \
			{'name' : '_'.join(['mapping_summary',eachsample]),
			'status' : 'waiting',
			'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['mapping_summary']['t']),
			'memory' : compute_resource['mapping_summary']['m'],
			'cmd' : 'sh '+os.path.join(mymapdir,'_'.join(['mapping_summary',eachsample])+'.sh')}
		safe_open(os.path.join(mymapdir,'_'.join(['mapping_summary',eachsample])+'.sh'),'w').write(stat_cmd)
		job_points['mapping_summary'].append('_'.join(['mapping_summary',eachsample]))
	
if set(['1','2']).issubset(includes) and argv.assembly:
	print "Assembly ...\n"
	merge_assemdir = os.path.join(assemdir,'merge_gtf')
	filter_assemdir = os.path.join(assemdir,'lncRNA_filter')
	create_dir(assemdir)
	create_dir(merge_assemdir)
	create_dir(filter_assemdir)

	job_points['assembly'] = []
	job_points['merge_assembly'] = []
	if argv.seqtype == 'lncRNA':
		job_points['extract_fasta']= []
		job_points['split_fasta']= []
		job_points['lncRNA_cpc'] = []
		job_points['lncRNA_cnci'] = []
		job_points['lncRNA_pfam'] = []
		job_points['lncRNA_identification'] = []
		job_points['lncRNA_signatures'] = []
		job_points['merge_cpc'] = []
		job_points['merge_cnci'] = []
		job_points['merge_pfam'] = []
	for eachsample in list_in_sample:
		myqcdir = os.path.join(qcdir,eachsample)
		mymapdir = os.path.join(mapdir,eachsample)
		myassemdir = os.path.join(assemdir,eachsample)
		create_dir(myassemdir)
		aAssembly = Assembly(eachsample,mymapdir,myassemdir,assemdir,softwares,databases)

		assem_cmd,order1 = aAssembly.assembly(soft=argv.assembly_method,libtype=libtype)
		add_items(orders,order1)
		new_jobs['_'.join(['assembly',eachsample])] = \
			{'name' : '_'.join(['assembly',eachsample]),
			'status' : 'waiting',
			'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['assembly']['t']),
			'memory' : compute_resource['assembly']['m'],
			'cmd' : 'sh '+os.path.join(myassemdir,'_'.join(['assembly',eachsample])+'.sh')}
		safe_open(os.path.join(myassemdir,'_'.join(['assembly',eachsample])+'.sh'),'w').write(assem_cmd)
		job_points['assembly'].append('_'.join(['assembly',eachsample]))

	if argv.seqtype == 'lncRNA':
		merge_mode = 'transcript'
	else:
		merge_mode = 'gene'


	if len(list_in_sample.keys())<5:
		qsub_merge='2G'
	elif len(list_in_sample.keys())<10:
		qsub_merge='5G'
	elif len(list_in_sample.keys())<20:
		qsub_merge='10G'
	elif len(list_in_sample.keys())<30:
		qsub_merge='20G'
	elif len(list_in_sample.keys())<50:
		qsub_merge='30G'
	elif len(list_in_sample.keys())<100:
		qsub_merge='50G'
	else:
		qsub_merge='100G'
		
	merge_cmd,order1 = aAssembly.merge_assembly(list_in_sample.keys(),soft=argv.assembly_method,mode=merge_mode)
	add_items(orders,order1)
	new_jobs['merge_assembly'] = \
		{'name' : 'merge_assembly',
		'status' : 'waiting',
		'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['merge_assembly']['t']),
		'memory' : qsub_merge,
		'cmd' : 'sh '+os.path.join(merge_assemdir,'merge_assembly'+'.sh')}
	safe_open(os.path.join(merge_assemdir,'merge_assembly'+'.sh'),'w').write(merge_cmd)
	job_points['merge_assembly'].append('merge_assembly')
	
	if argv.seqtype == 'lncRNA':
		cmd,order1 = aAssembly.extract_fasta()
		add_items(orders,order1)
		new_jobs['extract_fasta'] = \
			{'name' : 'extract_fasta',
			'status' : 'waiting',
			'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['extract_fasta']['t']),
			'memory' : compute_resource['extract_fasta']['m'],
			'cmd' : 'sh '+os.path.join(filter_assemdir,'extract_fasta'+'.sh')}
		safe_open(os.path.join(filter_assemdir,'extract_fasta'+'.sh'),'w').write(cmd)
		job_points['extract_fasta'].append('extract_fasta')
		fa_split_dir=os.path.join(filter_assemdir,'split')
		create_dir(fa_split_dir)
		cmd,order1 = aAssembly.split_fasta(fa_split_num=20)
		add_items(orders,order1)
		new_jobs['split_fasta'] = \
			{'name' : 'split_fasta',
			'status' : 'waiting',
			'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['split_fasta']['t']),
			'memory' : compute_resource['split_fasta']['m'],
			'cmd' : 'sh '+os.path.join(filter_assemdir,'split_fasta'+'.sh')}
		safe_open(os.path.join(filter_assemdir,'split_fasta'+'.sh'),'w').write(cmd)
		job_points['split_fasta'].append('split_fasta')
		cpc_dir=os.path.join(filter_assemdir,'CPC')
		create_dir(cpc_dir)
		cnci_dir=os.path.join(filter_assemdir,'CNCI')
		create_dir(cnci_dir)
		pfam_dir=os.path.join(filter_assemdir,'PFAM')
		create_dir(pfam_dir)
		for each in range(1,21):
			each=str(each)
			infa=os.path.join(filter_assemdir,'split','novel.exp_filter_'+each+'.fasta')
			each_cpc_dir=os.path.join(filter_assemdir,'CPC','split_'+each)
			create_dir(each_cpc_dir)
			each_cnci_dir=os.path.join(filter_assemdir,'CNCI','split_'+each)
			create_dir(each_cnci_dir)
			each_pfam_dir=os.path.join(filter_assemdir,'PFAM','split_'+each)
			create_dir(each_pfam_dir)
			cmd,order1 = aAssembly.lncRNA_cpc(infa,each_cpc_dir,each)
			add_items(orders,order1)
			new_jobs['lncRNA_cpc_'+each] = \
				{'name' : 'lncRNA_cpc_'+each,
				'status' : 'waiting',
				'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['lncRNA_cpc']['t']),
				'memory' : compute_resource['lncRNA_cpc']['m'],
				'cmd' : 'sh '+os.path.join(each_cpc_dir,'lncRNA_cpc_'+each+'.sh')}
			safe_open(os.path.join(each_cpc_dir,'lncRNA_cpc_'+each+'.sh'),'w').write(cmd)
			job_points['lncRNA_cpc'].append(os.path.join('lncRNA_cpc_',each))
			cmd,order1 = aAssembly.lncRNA_cnci(infa,each_cnci_dir,each)
			add_items(orders,order1)
			new_jobs['lncRNA_cnci_'+each] = \
				{'name' : 'lncRNA_cnci_'+each,
				'status' : 'waiting',
				'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['lncRNA_cnci']['t']),
				'memory' : compute_resource['lncRNA_cnci']['m'],
				'cmd' : 'sh '+os.path.join(each_cnci_dir,'lncRNA_cnci_'+each+'.sh')}
			safe_open(os.path.join(each_cnci_dir,'lncRNA_cnci_'+each+'.sh'),'w').write(cmd)
			job_points['lncRNA_cnci'].append(os.path.join('lncRNA_cnci_',each))
			cmd,order1 = aAssembly.lncRNA_pfam(infa,each_pfam_dir,each)
			add_items(orders,order1)
			new_jobs['lncRNA_pfam_'+each] = \
				{'name' : 'lncRNA_pfam_'+each,
				'status' : 'waiting',
				'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['lncRNA_pfam']['t']),
				'memory' : compute_resource['lncRNA_pfam']['m'],
				'cmd' : 'sh '+os.path.join(each_pfam_dir,'lncRNA_pfam_'+each+'.sh')}
			safe_open(os.path.join(each_pfam_dir,'lncRNA_pfam_'+each+'.sh'),'w').write(cmd)
			job_points['lncRNA_pfam'].append(os.path.join('lncRNA_pfam_',each))
		cmd,order1 = aAssembly.merge_cpc()
		add_items(orders,order1)
		new_jobs['merge_cpc'] = \
			{'name' : 'merge_cpc',
			'status' : 'waiting',
			'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['merge_cpc']['t']),
			'memory' : compute_resource['merge_cpc']['m'],
			'cmd' : 'sh '+os.path.join(cpc_dir,'merge_cpc'+'.sh')}
		safe_open(os.path.join(cpc_dir,'merge_cpc'+'.sh'),'w').write(cmd)
		job_points['merge_cpc'] = ['merge_cpc']
		cmd,order1 = aAssembly.merge_cnci()
		add_items(orders,order1)
		new_jobs['merge_cnci'] = \
			{'name' : 'merge_cnci',
			'status' : 'waiting',
			'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['merge_cnci']['t']),
			'memory' : compute_resource['merge_cnci']['m'],
			'cmd' : 'sh '+os.path.join(cnci_dir,'merge_cnci'+'.sh')}
		safe_open(os.path.join(cnci_dir,'merge_cnci'+'.sh'),'w').write(cmd)
		job_points['merge_cnci'] = ['merge_cnci']
		cmd,order1 = aAssembly.merge_pfam()
		add_items(orders,order1)
		new_jobs['merge_pfam'] = \
			{'name' : 'merge_pfam',
			'status' : 'waiting',
			'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['merge_pfam']['t']),
			'memory' : compute_resource['merge_pfam']['m'],
			'cmd' : 'sh '+os.path.join(pfam_dir,'merge_pfam'+'.sh')}
		safe_open(os.path.join(pfam_dir,'merge_pfam'+'.sh'),'w').write(cmd)
		job_points['merge_pfam'] = ['merge_pfam']
		cmd,order1 = aAssembly.lncRNA_identification()
		add_items(orders,order1)
		new_jobs['lncRNA_identification'] = \
			{'name' : 'lncRNA_identification',
			'status' : 'waiting',
			'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['lncRNA_identification']['t']),
			'memory' : compute_resource['lncRNA_identification']['m'],
			'cmd' : 'sh '+os.path.join(filter_assemdir,'lncRNA_identification'+'.sh')}
		safe_open(os.path.join(filter_assemdir,'lncRNA_identification'+'.sh'),'w').write(cmd)
		job_points['lncRNA_identification'].append('lncRNA_identification')
		cmd,order1 = aAssembly.lncRNA_signatures()
		add_items(orders,order1)
		new_jobs['lncRNA_signatures'] = \
			{'name' : 'lncRNA_signatures',
			'status' : 'waiting',
			'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['lncRNA_signatures']['t']),
			'memory' : compute_resource['lncRNA_signatures']['m'],
			'cmd' : 'sh '+os.path.join(filter_assemdir,'lncRNA_signatures'+'.sh')}
		safe_open(os.path.join(filter_assemdir,'lncRNA_signatures'+'.sh'),'w').write(cmd)
		job_points['lncRNA_signatures'].append('lncRNA_signatures')

## Quantification
if set(['1','2','3']).issubset(includes):
	print "Quantification ...\n"
	create_dir(expdir)
		
	gtf_quant = os.path.join(expdir,'all_transcripts.gtf')
	create_dir(os.path.join(expdir,'Summary'))
	if argv.assembly and exp_soft == 'rsem':
		create_dir(os.path.join(expdir,'RSEM_INDEX'))
		job_points['prepare_rsem_index'] = []
	job_points["quantification"] = []
	job_points["quantification_FPKM_ratio"] = []
	job_points['quantification_summary'] = []
	if exp_soft=='cuffdiff' and repeat_flag:
		job_points['cuffquant']=[]
    
	for eachsample in list_in_sample_noNEG:
		myqcdir = os.path.join(qcdir,eachsample)
		mymapdir = os.path.join(mapdir,eachsample)
		myexpdir = os.path.join(expdir,eachsample)
		aQuant = Quantification(eachsample,myqcdir,mymapdir,myexpdir,softwares,databases,dic_disease)
		if exp_soft=='cuffdiff' and repeat_flag :
			create_dir(myexpdir)
			cmd,order1 = aQuant.cuffquant(libtype=libtype)
			add_items(orders,order1)
			new_jobs['_'.join(['cuffquant',eachsample])] = \
				{'name' : '_'.join(['cuffquant',eachsample]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['cuffquant_quantification']['t']),
				'memory' : compute_resource['cuffquant_quantification']['m'],
				'cmd' : 'sh '+os.path.join(myexpdir,'_'.join(['cuffquant',eachsample])+'.sh')}
			safe_open(os.path.join(myexpdir,'_'.join(['cuffquant',eachsample])+'.sh'),'w').write(cmd)
			job_points['cuffquant'].append('_'.join(['cuffquant',eachsample]))
		if exp_soft!='cuffdiff':
			create_dir(myexpdir)
			cmd,order1 = aQuant.quantification(gtf_quant,argv.assembly,soft=exp_soft,analysis_type=argv.seqtype,libtype=libtype,repeat_flag=repeat_flag)
			add_items(orders,order1)
			new_jobs['_'.join(['quantification',eachsample])] = \
				{'name' : '_'.join(['quantification',eachsample]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['quantification']['t']),
				'memory' : compute_resource['quantification']['m'],
				'cmd' : 'sh '+os.path.join(myexpdir,'_'.join(['quantification',eachsample])+'.sh')}
			safe_open(os.path.join(myexpdir,'_'.join(['quantification',eachsample])+'.sh'),'w').write(cmd)
			job_points['quantification'].append('_'.join(['quantification',eachsample]))
			cmd,order1 = aQuant.htseq_quantification_FPKM_ratio(libtype=libtype)
			add_items(orders,order1)
			new_jobs['_'.join(['quantification_FPKM_ratio',eachsample])] = \
				{'name' : '_'.join(['quantification_FPKM_ratio',eachsample]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['quantification_FPKM_ratio']['t']),
				'memory' : compute_resource['quantification_FPKM_ratio']['m'],
				'cmd' : 'sh '+os.path.join(myexpdir,'_'.join(['quantification_FPKM_ratio',eachsample])+'.sh')}
			safe_open(os.path.join(myexpdir,'_'.join(['quantification_FPKM_ratio',eachsample])+'.sh'),'w').write(cmd)
			job_points['quantification_FPKM_ratio'].append('_'.join(['quantification_FPKM_ratio',eachsample]))

	if argv.assembly:
		if argv.seqtype == 'lncRNA':
			cmd,order1 = aQuant.quant_prepare(assembly=True,assemdir=assemdir,analysis_type='lncRNA',sample_lis=list_in_sample.keys())
		elif argv.seqtype == 'RNAseq':
			cmd,order1 = aQuant.quant_prepare(assembly=True,assemdir=assemdir,analysis_type='RNAseq',sample_lis=list_in_sample.keys())
	else:
		cmd,order1 = aQuant.quant_prepare(assembly=False,assemdir=assemdir,sample_lis=list_in_sample.keys())
	add_items(orders,order1)
	job_points['quant_prepare']=[]
	new_jobs['quant_prepare'] = \
		{'name' : 'quant_prepare',
		'status' : 'waiting',
		'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['quant_prepare']['t']),
		'memory' : compute_resource['quant_prepare']['m'],
		'cmd' : 'sh '+os.path.join(expdir,'quant_prepare.sh')}
	safe_open(os.path.join(expdir,'quant_prepare.sh'),'w').write(cmd)
	job_points['quant_prepare'].append('quant_prepare')


	if argv.assembly and exp_soft == 'rsem':
		cmd,order1 = aQuant.prepare_rsem_index()
		add_items(orders,order1)
		new_jobs['prepare_rsem_index'] = \
			{'name' : 'prepare_rsem_index',
			'status' : 'waiting',
			'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['prepare_rsem_index']['t']),
			'memory' : compute_resource['prepare_rsem_index']['m'],
			'cmd' : 'sh '+os.path.join(expdir,'RSEM_INDEX','prepare_rsem_index.sh')}
		safe_open(os.path.join(expdir,'RSEM_INDEX','prepare_rsem_index.sh'),'w').write(cmd)
		job_points['prepare_rsem_index'].append('prepare_rsem_index')
	if exp_soft == 'cuffdiff':
		if len(list_in_sample.keys()) < 10:
			qsub_cuffdiff='5G'
		elif len(list_in_sample.keys()) < 20:
			qsub_cuffdiff='10G'
		elif len(list_in_sample.keys()) < 30:
			qsub_cuffdiff='25G'
		elif len(list_in_sample.keys()) < 50:
			qsub_cuffdiff='50G'
		else:
			qsub_cuffdiff='100G'
		cuffdiff_dir=os.path.join(expdir,'cuffdiff')
		create_dir(cuffdiff_dir)
		job_points['cuffdiff_quantification']=[]	
		if argv.assembly:
			if argv.seqtype == 'lncRNA':
				cmd,order1 = aQuant.cuffdiff_quantification(sample_lis=list_in_sample.keys(),groupname=','.join(all_groupname),groups=all_groups,libtype=libtype,assembly=True,analysis_type="lncRNA",repeat_flag=repeat_flag) 
			elif argv.seqtype == 'RNAseq':
				cmd,order1 = aQuant.cuffdiff_quantification(sample_lis=list_in_sample.keys(),groupname=','.join(all_groupname),groups=all_groups,libtype=libtype,assembly=True,analysis_type="RNAseq",repeat_flag=repeat_flag)
		else:
			cmd,order1 = aQuant.cuffdiff_quantification(sample_lis=list_in_sample.keys(),groupname=','.join(all_groupname),groups=all_groups,libtype=libtype,assembly=False,repeat_flag=repeat_flag)
		add_items(orders,order1)
		new_jobs['cuffdiff_quantification'] = \
			{'name' : 'cuffdiff_quantification',
			'status' : 'waiting',
			'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['cuffdiff_quantification']['t']),
			'memory' : qsub_cuffdiff,
			'cmd' : 'sh '+os.path.join(expdir,'cuffdiff_quantification.sh')}
		safe_open(os.path.join(expdir,'cuffdiff_quantification.sh'),'w').write(cmd)
		job_points['cuffdiff_quantification'].append('cuffdiff_quantification')

	if argv.seqtype=='lncRNA':
		mod_flag='both'
	else:
		mod_flag='gene'

	cmd,order1 = aQuant.quantification_summary(sample_lis=list_in_sample_noNEG.keys(),groups=','.join(all_groups),groupname=','.join(all_groupname),mod=mod_flag,seqtype=argv.seqtype,exp_soft=exp_soft,repeat_flag=repeat_flag)	
	add_items(orders,order1)
	new_jobs['quantification_summary'] = \
		{'name' : 'quantification_summary',
		'status' : 'waiting',
		'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['quantification_summary']['t']),
		'memory' : compute_resource['quantification_summary']['m'],
		'cmd' : 'sh '+os.path.join(expdir,'Summary','quantification_summary.sh')}
	safe_open(os.path.join(expdir,'Summary','quantification_summary.sh'),'w').write(cmd)
	job_points['quantification_summary'].append('quantification_summary')
		
	if argv.WGCNA:
		job_points['wgcna']=[]
		wgcna_dir=os.path.join(expdir,'WGCNA')
		create_dir(wgcna_dir)

		cmd,order1 = aQuant.wgcna()
		add_items(orders,order1)
		new_jobs['wgcna'] = \
			{'name' : 'wgcna',
			'status' : 'waiting',
			'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['wgcna']['t']),
			'memory' : compute_resource['wgcna']['m'],
			'cmd' : 'sh '+os.path.join(expdir,'WGCNA','wgcna.sh')}
		safe_open(os.path.join(expdir,'WGCNA','wgcna.sh'),'w').write(cmd)
		job_points['wgcna'].append('wgcna')
            

## RNAvariants
if set(['1','2','6']).issubset(includes):
	print "SNP/INDEL ...\n"
	create_dir(mutdir)
	create_dir(pindeldir)
	create_dir(conpairdir)
	create_dir(zhushidir)
	job_points['process_bam'] = []
	job_points['mutation_calling'] = []
	job_points['somatic_Pindel'] = []
	job_points['somatic_Pindel_annovar'] = []
	job_points['run_conpair'] = [] 
	#job_points['summary_conpair'] = []
	job_points['zhushi']=[]
	
	os.system('rm -rf  %s/filter_cmd.sh'%(zhushidir))
	
	for eachsample in list_in_sample:
		myqcdir = os.path.join(qcdir,eachsample)
		mymapdir = os.path.join(mapdir,eachsample)
		mymutdir = os.path.join(mutdir,eachsample)
		create_dir(mymutdir)
		aVariant = RNAVariants(eachsample,myqcdir,mymapdir,mymutdir,softwares,databases,tumor_type)

		cmd,order1 = aVariant.process_bam()
		add_items(orders,order1)
		new_jobs['_'.join(['process_bam',eachsample])] = \
			{'name' : '_'.join(['process_bam',eachsample]),
			'status' : 'waiting',
			'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['process_bam']['t']),
			'memory' : compute_resource['process_bam']['m'],
			'cmd' : 'sh '+os.path.join(mymutdir,'_'.join(['process_bam',eachsample])+'.sh')}
		safe_open(os.path.join(mymutdir,'_'.join(['process_bam',eachsample])+'.sh'),'w').write(cmd)
		job_points['process_bam'].append('_'.join(['process_bam',eachsample]))

		cmd,order1 = aVariant.mutation_calling(TR,soft=mut_soft)
		add_items(orders,order1)
		new_jobs['_'.join(['mutation_calling',eachsample])] = \
			{'name' : '_'.join(['mutation_calling',eachsample]),
			'status' : 'waiting',
			'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['mutation_calling']['t']),
			'memory' : compute_resource['mutation_calling']['m'],
			'cmd' : 'sh '+os.path.join(mymutdir,'_'.join(['mutation_calling',eachsample])+'.sh')}
		safe_open(os.path.join(mymutdir,'_'.join(['mutation_calling',eachsample])+'.sh'),'w').write(cmd)
		job_points['mutation_calling'].append('_'.join(['mutation_calling',eachsample]))

	##################  20250905 add Pindel
	for eachsample in list_in_sample_noNEG:
		mypindeldir = os.path.join(pindeldir,eachsample)
		create_dir(mypindeldir)

		cmd,order1 = aVariant.somatic_Pindel(mypindeldir,eachsample,mutdir)
		add_items(orders,order1)
		new_jobs['_'.join(['somatic_Pindel',eachsample])] = \
			{'name' : '_'.join(['somatic_Pindel',eachsample]),
			'status' : 'waiting',
			'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['somatic_Pindel']['t']),
			'memory' : compute_resource['somatic_Pindel']['m'],
			'cmd' : 'sh '+os.path.join(mypindeldir,'_'.join(['somatic_Pindel',eachsample])+'.sh')}
		safe_open(os.path.join(mypindeldir,'_'.join(['somatic_Pindel',eachsample])+'.sh'),'w').write(cmd)
		job_points['somatic_Pindel'].append('_'.join(['somatic_Pindel',eachsample]))
		
		cmd,order1 = aVariant.somatic_Pindel_annovar(mypindeldir,eachsample)
		add_items(orders,order1)
		new_jobs['_'.join(['somatic_Pindel_annovar',eachsample])] = \
			{'name' : '_'.join(['somatic_Pindel_annovar',eachsample]),
			'status' : 'waiting',
			'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['somatic_Pindel_annovar']['t']),
			'memory' : compute_resource['somatic_Pindel_annovar']['m'],
			'cmd' : 'sh '+os.path.join(mypindeldir,'_'.join(['somatic_Pindel_annovar',eachsample])+'.sh')}
		safe_open(os.path.join(mypindeldir,'_'.join(['somatic_Pindel_annovar',eachsample])+'.sh'),'w').write(cmd)
		job_points['somatic_Pindel_annovar'].append('_'.join(['somatic_Pindel_annovar',eachsample]))
 
 
	### 20250908 add conpair 
	for eachsample in DNARNA_combine.keys():
		myConpair = os.path.join(conpairdir,eachsample)
		create_dir(myConpair)
		
		cmd,order1 = aVariant.run_conpair(myConpair,NP,eachsample,mutdir)
		add_items(orders,order1)
		new_jobs['_'.join(['run_conpair',eachsample])] = \
		    {'name' : '_'.join(['run_conpair',eachsample]),
		    'status' : 'waiting',
		    'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['run_conpair']['t']),
		    'memory' : compute_resource['run_conpair']['m'],
		    'cmd' : 'sh '+os.path.join(myConpair,'_'.join(['run_conpair',eachsample])+'.sh')}
		safe_open(os.path.join(myConpair,'_'.join(['run_conpair',eachsample])+'.sh'),'w').write(cmd)
		job_points['run_conpair'].append('_'.join(['run_conpair',eachsample]))
	
		#cmd,order1 = aVariant.summary_conpair(conpairdir,myConpair,eachsample)
		#add_items(orders,order1)
		#new_jobs['summary_conpair'] = \
		#    {'name' : 'summary_conpair',
		#    'status' : 'waiting',
		#    'sched' : '-V -cwd %s' % (queue_list),
		#    'memory' : compute_resource['summary_conpair']['m'],
		#    'cmd' : 'sh ' + os.path.join(conpairdir,'summary_conpair.sh')}
		#safe_open(os.path.join(conpairdir,'summary_conpair.sh'),'w').write(cmd)
		#job_points['summary_conpair'] = ['summary_conpair']

		

	### 只加一次环境变量
	filter_cmd0 = softwares['perl_library'] + '\ncd ' + zhushidir + "\n"
	safe_open(os.path.join(zhushidir,'filter_cmd.sh'),'a').write(filter_cmd0)
	#print "==>",t_type
	## 遍历type---zhushi
	for i in t_type.keys():
		if i == "BRCA" or i == "CM" or i == "CORE" or i == "GIST" or i == "LICH" or i == "NSCLC" or i == "OST" or i == "OV" or i == "PAAD" or i == "STES"  or i == "THCA" or i == "UNKNOW":
			zhushi_cmd,order1 = aVariant.zhushi_solid(zhushidir,analydir,t_type,i,TR)
		else:
			zhushi_cmd,order1 = aVariant.zhushi(zhushidir,analydir,t_type,i,TR)
		add_items(orders,order1)
		new_jobs['_'.join(['zhushi',i])] = \
			{'name' : 'zhushi_%s'%i,
			'status' : 'waiting',
			'sched' : '-V -cwd %s' % queue_list,
			'memory' : compute_resource['zhushi']['m'],
			'cmd' : 'sh '+os.path.join(zhushidir,'_'.join(['zhushi',i])+'.sh')}
		safe_open(os.path.join(zhushidir,'_'.join(['zhushi',i])+'.sh'),'w').write(zhushi_cmd)
		job_points['zhushi'].append('_'.join(['zhushi',i]))
		job_points['filter_cmd'] = []
	
		### 遍历type里的每个sample-----filter_cmd
		for each in t_type[i]:
			if sample_type[each]=='cfDNA'  or sample_type[each]=='PM':
				type_tmp = 'PM'
			elif sample_type[each]=='UMI':
				type_tmp = 'UMI'
			else:
				type_tmp = 'FM'
			if re.search("NEG",each):
				filter_cmd1,order1 = aVariant.filter_cmd_NEG(zhushidir,i,each,dic_disease,analydir)
			else:
				if each.startswith("RT"):
					filter_cmd1,order1 = aVariant.filter_cmd_solid(zhushidir,i,each,dic_disease,analydir,type_tmp)
				else:
					filter_cmd1,order1 = aVariant.filter_cmd(zhushidir,i,each,dic_disease,analydir)
			add_items(orders,order1)
			new_jobs['filter_cmd'] = \
				{'name' : 'filter_cmd',
				'status' : 'waiting',
				'sched' : '-V -cwd %s' % queue_list,
				'memory' : compute_resource['filter_cmd']['m'],
				'cmd' : 'sh '+os.path.join(zhushidir,'filter_cmd.sh')}
			safe_open(os.path.join(zhushidir,'filter_cmd.sh'),'a').write(filter_cmd1)
			job_points['filter_cmd'] = ['filter_cmd']
## fusiongene
if set(['1','2','7']).issubset(includes):
	print "Fusion gene ...\n"
	create_dir(arribadir)
	create_dir(starfusiondir)
	if argv.genome=='human_B37':
		fusion_version='GRCh37'
	elif argv.genome=='human_B38':
		fusion_version='GRCh38'
	elif argv.genome=='mm10': 
		fusion_version='GRCm38'

    #### starfusion && Arriba
	job_points['fusion_step1'] = []
	job_points['fusion_step2'] = []
	job_points['arriba_step1'] = []
	job_points['arriba_step2'] = []
	
	for eachsample in list_in_sample:
		myqcdir = os.path.join(qcdir,eachsample)
		mymapdir = os.path.join(mapdir,eachsample)

		### starfusion
		myfusiondir1 = os.path.join(starfusiondir,eachsample)
		create_dir(myfusiondir1)
		aFusiongene = Fusiongene(eachsample,myqcdir,mymapdir,myfusiondir1,softwares,databases,tumor_type,analydir,NP,dic_disease)
		
		cmd,order1 = aFusiongene.fusion_step1()
		add_items(orders,order1)
		new_jobs['_'.join(['fusion_step1',eachsample])] = \
			{'name' : '_'.join(['fusion_step1',eachsample]),
			'status' : 'waiting',
			'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['fusion_step1']['t']),
			'memory' : compute_resource['fusion_step1']['m'],
			'cmd' : 'sh '+os.path.join(myfusiondir1,'_'.join(['fusion_step1',eachsample])+'.sh')}
		safe_open(os.path.join(myfusiondir1,'_'.join(['fusion_step1',eachsample])+'.sh'),'w').write(cmd)
		job_points['fusion_step1'].append('_'.join(['fusion_step1',eachsample]))

		cmd,order1 = aFusiongene.fusion_step2(version=fusion_version,sample_db=sample_db)
		add_items(orders,order1)
		new_jobs['_'.join(['fusion_step2',eachsample])] = \
			{'name' : '_'.join(['fusion_step2',eachsample]),
			'status' : 'waiting',
			'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['fusion_step2']['t']),
			'memory' : compute_resource['fusion_step2']['m'],
			'cmd' : 'sh '+os.path.join(myfusiondir1,'_'.join(['fusion_step2',eachsample])+'.sh')}
		safe_open(os.path.join(myfusiondir1,'_'.join(['fusion_step2',eachsample])+'.sh'),'w').write(cmd)
		job_points['fusion_step2'].append('_'.join(['fusion_step2',eachsample]))

		### Arriba
	for eachsample in list_in_sample_noNEG:
		myqcdir = os.path.join(qcdir,eachsample)
		myfusiondir3 = os.path.join(arribadir,eachsample)
		aFusiongene = Fusiongene(eachsample,myqcdir,mymapdir,myfusiondir3,softwares,databases,tumor_type,analydir,NP,dic_disease)
		create_dir(myfusiondir3)
		
		cmd,order1 = aFusiongene.arriba_step1(eachsample)
		add_items(orders,order1)
		new_jobs['_'.join(['arriba_step1',eachsample])] = \
			{'name' : '_'.join(['arriba_step1',eachsample]),
			'status' : 'waiting',
			'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['arriba_step1']['t']),
			'memory' : compute_resource['arriba_step1']['m'],
			'cmd' : 'sh '+os.path.join(myfusiondir3,'_'.join(['arriba_step1',eachsample])+'.sh')}
		safe_open(os.path.join(myfusiondir3,'_'.join(['arriba_step1',eachsample])+'.sh'),'w').write(cmd)
		job_points['arriba_step1'].append('_'.join(['arriba_step1',eachsample]))
		
		cmd,order1 = aFusiongene.arriba_step2(sample_db,eachsample)
		add_items(orders,order1)
		new_jobs['_'.join(['arriba_step2',eachsample])] = \
			{'name' : '_'.join(['arriba_step2',eachsample]),
			'status' : 'waiting',
			'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['arriba_step2']['t']),
			'memory' : compute_resource['arriba_step2']['m'],
			'cmd' : 'sh '+os.path.join(myfusiondir3,'_'.join(['arriba_step2',eachsample])+'.sh')}
		safe_open(os.path.join(myfusiondir3,'_'.join(['arriba_step2',eachsample])+'.sh'),'w').write(cmd)
		job_points['arriba_step2'].append('_'.join(['arriba_step2',eachsample]))



## TTMV_screen
if set(['1','2','9']).issubset(includes):
	print "TTMV::RARA病毒整合分析 ...\n"
	create_dir(TTMVscreendir)
	if argv.genome=='human_B37':
		fusion_version='GRCh37'
	elif argv.genome=='human_B38':
		fusion_version='GRCh38'
	elif argv.genome=='mm10':
		fusion_version='GRCm38'
	
	run_TTMVscreen0 = softwares['perl_library'] + '\ncd ' + TTMVscreendir + "\n"
	safe_open(os.path.join(TTMVscreendir,'run_TTMVscreen.sh'),'w').write(run_TTMVscreen0)

	job_points['run_TTMVscreen'] = []
	for eachsample in list_in_sample:
		myqcdir = os.path.join(qcdir,eachsample)
		mymapdir = os.path.join(mapdir,eachsample)

		aVariant = RNAVariants(eachsample,myqcdir,mymapdir,"",softwares,databases,tumor_type)	
		cmd,order1 = aVariant.TTMVscreen(eachsample,mymapdir,TTMVscreendir) 
		add_items(orders,order1)

		new_jobs['run_TTMVscreen'] = \
		    {'name' : 'run_TTMVscreen',
		    'status' : 'waiting',
		    'sched' : '-V -cwd %s' % queue_list,
		    'memory' : compute_resource['run_TTMVscreen']['m'],
		    'cmd' : 'sh '+os.path.join(TTMVscreendir,'run_TTMVscreen.sh')}
		safe_open(os.path.join(TTMVscreendir,'run_TTMVscreen.sh'),'a').write(cmd)
		job_points['run_TTMVscreen'] = ['run_TTMVscreen']



if set(['1','2','3','4']).issubset(includes):
	print "Different ...\n"
	create_dir(diffdir)
	job_points['diff_expression'] = []
	job_points['adjusting_parameter'] = []
	job_points['go_enrichment'] = []
	job_points['kegg_enrichment'] = []
	job_points['reactome_enrichment'] = []
	job_points['classification_of_RNA']=[]
	job_points[''] = []

	gene_info=os.path.join(diffdir,'geneInfo')
	condition_file=safe_open(os.path.join(diffdir,'condition.txt'),'w')
	condition_file.write('samples\tgroups1\ttype\n')
        condition_lis=[]
	tiaocan_list=[]
	class_list=[]

	if argv.seqtype=='lncRNA':
		lncrna_flag='yes'
		mod_flag='both'
	else:
		lncrna_flag='no'
		mod_flag='gene'
	for eachgrp in this_compare_groups:
		mydiffdir = os.path.join(diffdir,eachgrp)
		create_dir(mydiffdir)
		create_dir(os.path.join(mydiffdir,'GO'))
		create_dir(os.path.join(mydiffdir,'KEGG'))
		create_dir(os.path.join(mydiffdir,'KEGG','pathway'))
		create_dir(os.path.join(mydiffdir,'Reactome'))
		Ts = this_compare_groups[eachgrp]['T']['samples'].split(',')
		Ns = this_compare_groups[eachgrp]['N']['samples'].split(',')
		Tname = this_compare_groups[eachgrp]['T']['name']
		for each_T in Ts:
                        ele=each_T+'\t'+Tname
                        if ele not in condition_lis:
				condition_file.write(each_T+'\t'+Tname+'\n')
                        	condition_lis.append(ele)
		Nname = this_compare_groups[eachgrp]['N']['name']
		for each_N in Ns:
			ele=each_N+'\t'+Nname
			if ele not in condition_lis:
				condition_file.write(each_N+'\t'+Nname+'\n')
				condition_lis.append(ele)

		aDiff = Differential(Ts,Ns,Tname,Nname,eachgrp,mapdir,expdir,mydiffdir,softwares,databases)
		if diff_soft!='cuffdiff':
			cmd,order1 = aDiff.diff_expression(gene_info,soft=diff_soft,mod=mod_flag,lncRNA=lncrna_flag)
			cmd_tiaocan,order1_tiaocan = aDiff.adjusting_parameter(mod=mod_flag,FC='2',porq='pval',value='0.05')
			tiaocan_list.append(cmd_tiaocan)
			add_items(orders,order1)
			new_jobs['_'.join(['diff_expression',eachgrp])] = \
				{'name' : '_'.join(['diff_expression',eachgrp]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['diff_expression']['t']),
				'memory' : compute_resource['diff_expression']['m'],
				'cmd' : 'sh '+os.path.join(mydiffdir,'_'.join(['diff_expression',eachgrp])+'.sh')}
			safe_open(os.path.join(mydiffdir,'_'.join(['diff_expression',eachgrp])+'.sh'),'w').write(cmd)
			job_points['diff_expression'].append('_'.join(['diff_expression',eachgrp]))
		
		cmd,order1 = aDiff.classification_of_RNA(mod=mod_flag,soft=diff_soft,lncRNA=lncrna_flag)
		class_cmd,class_order1 = aDiff.classification_of_RNA(mod=mod_flag,soft=diff_soft,lncRNA=lncrna_flag,porq='pval')
		class_list.append(class_cmd)
		add_items(orders,order1)
		new_jobs['_'.join(['classification_of_RNA',eachgrp])] = \
			{'name' : '_'.join(['classification_of_RNA',eachgrp]),
			'status' : 'waiting',
			'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['classification_of_RNA']['t']),
			'memory' : compute_resource['classification_of_RNA']['m'],
			'cmd' : 'sh '+os.path.join(mydiffdir,'_'.join(['classification_of_RNA',eachgrp])+'.sh')}
		safe_open(os.path.join(mydiffdir,'_'.join(['classification_of_RNA',eachgrp])+'.sh'),'w').write(cmd)
		job_points['classification_of_RNA'].append('_'.join(['classification_of_RNA',eachgrp]))

		if enrich_spe in ['hsa','mmu']:
	
			cmd,order1 = aDiff.go_enrichment(enrich_spe,lncrna_flag)
			add_items(orders,order1)
			new_jobs['_'.join(['go_enrichment',eachgrp])] = \
				{'name' : '_'.join(['go_enrichment',eachgrp]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['go_enrichment']['t']),
				'memory' : compute_resource['go_enrichment']['m'],
				'cmd' : 'sh '+os.path.join(mydiffdir,'_'.join(['go_enrichment',eachgrp])+'.sh')}
			safe_open(os.path.join(mydiffdir,'_'.join(['go_enrichment',eachgrp])+'.sh'),'w').write(cmd)
			job_points['go_enrichment'].append('_'.join(['go_enrichment',eachgrp]))

			cmd,order1 = aDiff.kegg_enrichment(enrich_spe,lncrna_flag)
			add_items(orders,order1)
			new_jobs['_'.join(['kegg_enrichment',eachgrp])] = \
				{'name' : '_'.join(['kegg_enrichment',eachgrp]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['kegg_enrichment']['t']),
				'memory' : compute_resource['kegg_enrichment']['m'],
				'cmd' : 'sh '+os.path.join(mydiffdir,'_'.join(['kegg_enrichment',eachgrp])+'.sh')}
			safe_open(os.path.join(mydiffdir,'_'.join(['kegg_enrichment',eachgrp])+'.sh'),'w').write(cmd)
			job_points['kegg_enrichment'].append('_'.join(['kegg_enrichment',eachgrp]))

			cmd,order1 = aDiff.reactome_enrichment(enrich_spe,lncrna_flag)
			add_items(orders,order1)
			new_jobs['_'.join(['reactome_enrichment',eachgrp])] = \
				{'name' : '_'.join(['reactome_enrichment',eachgrp]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['reactome_enrichment']['t']),
				'memory' : compute_resource['reactome_enrichment']['m'],
				'cmd' : 'sh '+os.path.join(mydiffdir,'_'.join(['reactome_enrichment',eachgrp])+'.sh')}
			safe_open(os.path.join(mydiffdir,'_'.join(['reactome_enrichment',eachgrp])+'.sh'),'w').write(cmd)
			job_points['reactome_enrichment'].append('_'.join(['reactome_enrichment',eachgrp]))
		else:
			cmd,order1 = aDiff.go_enrichment2(enrich_spe,lncrna_flag)
			add_items(orders,order1)
			new_jobs['_'.join(['go_enrichment',eachgrp])] = \
				{'name' : '_'.join(['go_enrichment',eachgrp]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['go_enrichment']['t']),
				'memory' : compute_resource['go_enrichment']['m'],
				'cmd' : 'sh '+os.path.join(mydiffdir,'_'.join(['go_enrichment',eachgrp])+'.sh')}
			safe_open(os.path.join(mydiffdir,'_'.join(['go_enrichment',eachgrp])+'.sh'),'w').write(cmd)
			job_points['go_enrichment'].append('_'.join(['go_enrichment',eachgrp]))

			cmd,order1 = aDiff.kegg_enrichment2(enrich_spe,lncrna_flag)
			add_items(orders,order1)
			new_jobs['_'.join(['kegg_enrichment',eachgrp])] = \
				{'name' : '_'.join(['kegg_enrichment',eachgrp]),
 				'status' : 'waiting',
				'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['kegg_enrichment']['t']),
				'memory' : compute_resource['kegg_enrichment']['m'],
				'cmd' : 'sh '+os.path.join(mydiffdir,'_'.join(['kegg_enrichment',eachgrp])+'.sh')}
			safe_open(os.path.join(mydiffdir,'_'.join(['kegg_enrichment',eachgrp])+'.sh'),'w').write(cmd)
			job_points['kegg_enrichment'].append('_'.join(['kegg_enrichment',eachgrp]))

			
	if diff_soft != 'cuffdiff':
		job_points['diff_prepare']=[]
		cmd,order1 = aDiff.diff_prepare(samples_lis=list_in_sample.keys(),groups=','.join(all_groups),groupname=','.join(all_groupname),mod=mod_flag,seqtype=argv.seqtype)
		add_items(orders,order1)
		new_jobs['diff_prepare'] = \
			{'name' : 'diff_prepare',
			'status' : 'waiting',
			'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['diff_prepare']['t']),
			'memory' : compute_resource['diff_prepare']['m'],
			'cmd' : 'sh '+os.path.join(diffdir,'diff_prepare.sh')}
		safe_open(os.path.join(diffdir,'diff_prepare.sh'),'w').write(cmd)
		job_points['diff_prepare'].append('diff_prepare')

	if diff_soft == 'cuffdiff':
		job_points['extract_cuffdiff']=[]
		cmd,order1 = aDiff.extract_cuffdiff(samples_lis=list_in_sample.keys(),groups=','.join(all_groups),groupnames=','.     join(all_groupname),compare=','.join(this_compare_groups),repeat_flag=repeat_flag)
		cmd_tiaocan,order1_tiaocan = aDiff.extract_cuffdiff(samples_lis=list_in_sample.keys(),groups=','.join(all_groups),groupnames=','.join(all_groupname),compare=','.join(this_compare_groups),repeat_flag=repeat_flag,pvalue='0.05')
		tiaocan_list.append(cmd_tiaocan+' && \\')
		add_items(orders,order1)
		new_jobs['extract_cuffdiff'] = \
			{'name' : 'extract_cuffdiff',
			'status' : 'waiting',
			'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['extract_cuffdiff']['t']),
			'memory' : compute_resource['extract_cuffdiff']['m'],
			'cmd' : 'sh '+os.path.join(diffdir,'extract_cuffdiff.sh')}
		safe_open(os.path.join(diffdir,'extract_cuffdiff.sh'),'w').write(cmd)
		job_points['extract_cuffdiff'].append('extract_cuffdiff')

	cluster_dir=os.path.join(diffdir,'cluster')
	create_dir(cluster_dir)
	job_points['cluster']=[]
	cmd,order1 = aDiff.cluster(compare=','.join(this_compare_groups),mod=mod_flag,lncRNA=lncrna_flag)
	add_items(orders,order1)
	new_jobs['cluster'] = \
		{'name' : 'cluster',
		'status' : 'waiting',
		'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['cluster']['t']),
		'memory' : compute_resource['cluster']['m'],
		'cmd' : 'sh '+os.path.join(diffdir,'cluster','cluster.sh')}
	safe_open(os.path.join(diffdir,'cluster','cluster.sh'),'w').write(cmd)
	job_points['cluster'].append('cluster')

	pca_dir=os.path.join(diffdir,'PCA')
	create_dir(pca_dir)
	job_points['PCA']=[]
	cmd,order1 = aDiff.PCA(mod=mod_flag,lncRNA=lncrna_flag)
	add_items(orders,order1)
	new_jobs['PCA'] = \
		{'name' : 'PCA',
		'status' : 'waiting',
		'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['cluster']['t']),
		'memory' : compute_resource['cluster']['m'],
		'cmd' : 'sh '+os.path.join(diffdir,'PCA.sh')}
	safe_open(os.path.join(diffdir,'PCA.sh'),'w').write(cmd)
	job_points['PCA'].append('PCA')

	
	##tiaocan#####################################
	cmd_add1=' && \\\n'.join(tiaocan_list)
	cmd_add2='\n'.join(class_list)
 
	if lncrna_flag=='yes':
		prefix='.lncRNA.difftranscript.xls'
		if diff_soft=='cufdiff':
			xls_list=['.diffgene.xls','.difftranscript.xls','.gene_level.Differential_analysis_results.xls','.lncRNA.diffgene.xls','.lncRNA.difftranscript.xls','.lncRNA.gene_level.Differential_analysis_results.xls','.lncRNA.transcript_level.Differential_analysis_results.xls','.mRNA.diffgene.xls','.mRNA.difftranscript.xls','.mRNA.gene_level.Differential_analysis_results.xls','.mRNA.transcript_level.Differential_analysis_results.xls','.transcript_level.Differential_analysis_results.xls','.Unclassified.diffgene.xls','.Unclassified.difftranscript.xls','.Unclassified.gene_level.Differential_analysis_results.xls','.Unclassified.transcript_level.Differential_analysis_results.xls']
		else:
			xls_list=['.diffgene.xls','.difftranscript.xls','.lncRNA.diffgene.xls','.lncRNA.difftranscript.xls','.lncRNA.gene_level.Differential_analysis_results.xls','.lncRNA.transcript_level.Differential_analysis_results.xls','.mRNA.diffgene.xls','.mRNA.difftranscript.xls','.mRNA.gene_level.Differential_analysis_results.xls','.mRNA.transcript_level.Differential_analysis_results.xls','.Unclassified.diffgene.xls','.Unclassified.difftranscript.xls','.Unclassified.gene_level.Differential_analysis_results.xls','.Unclassified.transcript_level.Differential_analysis_results.xls','.diffgene.original_output']
		png_list=['.lncRNA.diffgene.Volcanoplot.png','.lncRNA.difftranscript.Volcanoplo.png','.mRNA.diffgene.Volcanoplot.png','.mRNA.difftranscript.Volcanoplot.png']
		pdf_list=['.lncRNA.diffgene.Volcanoplot.pdf','.lncRNA.difftranscript.Volcanoplo.pdf','.mRNA.diffgene.Volcanoplot.pdf','.mRNA.difftranscript.Volcanoplot.pdf']
	else:
		
		prefix='.diffgene.xls'
		xls_list=['.diffgene.xls','.diffgene.original_output']
		png_list=['.diffgene.Volcanoplot.png']
		pdf_list=['.diffgene.Volcanoplot.pdf']

	all_compare=' '.join(['"'+elem+'"' for elem in this_compare_groups])
	mv_list=xls_list+png_list+pdf_list
	mv_=' '.join(['"'+elem+'"' for elem in mv_list])
	#min_num_count=5

	cmd='''
compare=(%s)
mv_list=(%s)
diff_soft='%s'
prefix='%s'
for i in `ls %s/*/*${prefix}`;do wc -l $i>>%s/diffsum.txt;done 
min_num=$(sort -nk1 %s/diffsum.txt|head -1|awk -F ' ' '{print $1}')
min_num=$(echo $min_num/1|bc)
rm %s/diffsum.txt

if [[ $min_num -lt 5 ]]; then
    mkdir %s/qvalue
    for i in  ${compare[@]}
    do
        mkdir %s/qvalue/$i
        mv %s/$i/classification_of_RNA_${i}.sh %s/qvalue/$i
        for j in ${mv_list[@]}
        do
            mv %s/$i/${i}${j} %s/qvalue/$i;
        done

        if  [[ $diff_soft != 'cuffdiff' ]]; then
            mv %s/$i/diff_expression_$i.sh %s/qvalue/$i;
        fi
    done

    if  [[ $diff_soft = 'cuffdiff' ]]; then
        mv %s/extract_cuffdiff.sh %s/qvalue
    fi

    %s
    %s

fi
	''' % (all_compare,mv_,diff_soft,prefix,diffdir,diffdir,diffdir,diffdir,diffdir,diffdir,diffdir,diffdir,diffdir,diffdir,diffdir,diffdir,diffdir,diffdir,cmd_add1,cmd_add2)
	order1 = ['order adjusting_parameter after classification_of_RNA_%s' % each  for each in this_compare_groups]
	add_items(orders,order1)
	new_jobs['adjusting_parameter'] = \
		{'name' : 'adjusting_parameter',
		'status' : 'waiting',
		'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['adjusting_parameter']['t']),
		'memory' : compute_resource['adjusting_parameter']['m'],
		'cmd' : 'sh '+os.path.join(diffdir,'adjusting_parameter.sh')}
	safe_open(os.path.join(diffdir,'adjusting_parameter.sh'),'w').write(cmd)
	job_points['adjusting_parameter'].append('adjusting_parameter')
	##############################################
	

## rMTAS, as different analysis
if set(['1','2','8']).issubset(includes):
	print "Altersplice ...\n"
	create_dir(as_diffdir)
	job_points['rnaseq_mats'] = []
	job_points['rnaseq_mats_anno'] = []

	for eachgrp in this_compare_groups:
		myas_diffdir = os.path.join(as_diffdir,eachgrp)
		create_dir(myas_diffdir)

		Ts = this_compare_groups[eachgrp]['T']['samples'].split(',')
		Ns = this_compare_groups[eachgrp]['N']['samples'].split(',')
		Tname = this_compare_groups[eachgrp]['T']['name']
		Nname = this_compare_groups[eachgrp]['N']['name']

		if len(Ts) >= 3 and len(Ts)==len(Ns):
			analysis='P'
		else:
			analysis='U'
		aASdiff = ASdifferential(Ts,Ns,Tname,Nname,eachgrp,mapdir,myas_diffdir,softwares,databases)

		cmd,order1 = aASdiff.rnaseq_mats(analysis=analysis,len=argv.read_length,libtype=libtype)
		add_items(orders,order1)
		new_jobs['_'.join(['rnaseq_mats',eachgrp])] = \
			{'name' : '_'.join(['rnaseq_mats',eachgrp]),
			'status' : 'waiting',
			'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['rnaseq_mats']['t']),
			'memory' : compute_resource['rnaseq_mats']['m'],
			'cmd' : 'sh '+os.path.join(myas_diffdir,'_'.join(['rnaseq_mats',eachgrp])+'.sh')}
		safe_open(os.path.join(myas_diffdir,'_'.join(['rnaseq_mats',eachgrp])+'.sh'),'w').write(cmd)
		job_points['rnaseq_mats'].append('_'.join(['rnaseq_mats',eachgrp]))
		cmd2,order2 = aASdiff.rnaseq_mats_anno()
		add_items(orders,order2)
		new_jobs['_'.join(['rnaseq_mats_anno',eachgrp])] = \
			{'name' : '_'.join(['rnaseq_mats_anno',eachgrp]),
			'status' : 'waiting',
			'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['rnaseq_mats_anno']['t']),
			'memory' : compute_resource['rnaseq_mats_anno']['m'],
			'cmd' : 'sh '+os.path.join(myas_diffdir,'_'.join(['rnaseq_mats_anno',eachgrp])+'.sh')}
		safe_open(os.path.join(myas_diffdir,'_'.join(['rnaseq_mats_anno',eachgrp])+'.sh'),'w').write(cmd2)
		job_points['rnaseq_mats_anno'].append('_'.join(['rnaseq_mats_anno',eachgrp]))
			
## Lnc target
if set(['1','2','3','4','5']).issubset(includes) and argv.seqtype == 'lncRNA':
	print "LncRNA taregt  ...\n"
	create_dir(targetdir)
	create_dir(os.path.join(targetdir,'co_location'))
	create_dir(enrich_targetdir)
	job_points['co_location'] = []
	job_points['lnctarget_enrich_colocation_prepare'] = []
	job_points['colocation_go_enrichment']=[]
	job_points['colocation_kegg_enrichment']=[]
	job_points['colocation_reactome_enrichment']=[]
	if 5<=len(list_in_sample.keys())<=55:
		create_dir(os.path.join(targetdir,'co_expression'))
		job_points['co_expression'] = []
		job_points['lnctarget_enrich_coexpression_prepare'] = []
		job_points['coexpression_go_enrichment']=[]
		job_points['coexpression_kegg_enrichment']=[]
		job_points['coexpression_reactome_enrichment']=[]
	for eachgrp in this_compare_groups:
		mytargetdir = os.path.join(enrich_targetdir,eachgrp)
		create_dir(mytargetdir)
		create_dir(os.path.join(mytargetdir,'co_location'))
		Tname = this_compare_groups[eachgrp]['T']['name']
		Nname = this_compare_groups[eachgrp]['N']['name']
		aTarget = LncTarget(Tname,Nname,eachgrp,targetdir,mytargetdir,assemdir,expdir,diffdir,softwares,databases)

		cmd,order1 = aTarget.lnctarget_enrich_colocation_prepare(species=enrich_spe,diff_soft=diff_soft)
		add_items(orders,order1)
		new_jobs['_'.join(['lnctarget_enrich_colocation_prepare',eachgrp])] = \
			{'name' : '_'.join(['lnctarget_enrich_colocation_prepare',eachgrp]),
			'status' : 'waiting',
			'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['lnctarget_enrich_colocation_prepare']['t']),
			'memory' : compute_resource['lnctarget_enrich_colocation_prepare']['m'],
			'cmd' : 'sh '+os.path.join(mytargetdir,'_'.join(['lnctarget_enrich_colocation_prepare',eachgrp])+'.sh')}
		safe_open(os.path.join(mytargetdir,'_'.join(['lnctarget_enrich_colocation_prepare',eachgrp])+'.sh'),'w').write(cmd)
		job_points['lnctarget_enrich_colocation_prepare'].append('_'.join(['lnctarget_enrich_colocation_prepare',eachgrp]))

		co_location_go_dir=os.path.join(mytargetdir,'co_location','GO')
		co_location_kegg_dir=os.path.join(mytargetdir,'co_location','KEGG')
		co_location_reactome_dir=os.path.join(mytargetdir,'co_location','Reactome')
		create_dir(co_location_go_dir)
		create_dir(co_location_kegg_dir)
		create_dir(os.path.join(co_location_kegg_dir,'pathway'))
		create_dir(co_location_reactome_dir)
	
		if enrich_spe in ['hsa','mmu']:
		
			cmd = aTarget.go_enrichment(species=enrich_spe,label='co_location')
			order1 = 'order colocation_go_enrichment_%s after lnctarget_enrich_colocation_prepare_%s' % (eachgrp,eachgrp)
			add_items(orders,order1)
			new_jobs['_'.join(['colocation_go_enrichment',eachgrp])] = \
				{'name' : '_'.join(['colocation_go_enrichment',eachgrp]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['go_enrichment']['t']),
				'memory' : compute_resource['go_enrichment']['m'],
				'cmd' : 'sh '+os.path.join(mytargetdir,'co_location','_'.join(['colocation_go_enrichment',eachgrp])+'.sh')}
			safe_open(os.path.join(mytargetdir,'co_location','_'.join(['colocation_go_enrichment',eachgrp])+'.sh'),'w').write(cmd)
			job_points['colocation_go_enrichment'].append('_'.join(['colocation_go_enrichment',eachgrp]))			
		
			cmd = aTarget.kegg_enrichment(species=enrich_spe,label='co_location')
			order1 = 'order colocation_kegg_enrichment_%s after lnctarget_enrich_colocation_prepare_%s' % (eachgrp,eachgrp)
			add_items(orders,order1)
			new_jobs['_'.join(['colocation_kegg_enrichment',eachgrp])] = \
				{'name' : '_'.join(['colocation_kegg_enrichment',eachgrp]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['kegg_enrichment']['t']),
				'memory' : compute_resource['kegg_enrichment']['m'],
				'cmd' : 'sh '+os.path.join(mytargetdir,'co_location','_'.join(['colocation_kegg_enrichment',eachgrp])+'.sh')}
			safe_open(os.path.join(mytargetdir,'co_location','_'.join(['colocation_kegg_enrichment',eachgrp])+'.sh'),'w').write(cmd)
			job_points['colocation_kegg_enrichment'].append('_'.join(['colocation_kegg_enrichment',eachgrp]))
		
			cmd = aTarget.reactome_enrichment(species=enrich_spe,label='co_location')
			order1 = 'order colocation_reactome_enrichment_%s after lnctarget_enrich_colocation_prepare_%s' % (eachgrp,eachgrp)
			add_items(orders,order1)
			new_jobs['_'.join(['colocation_reactome_enrichment',eachgrp])] = \
				{'name' : '_'.join(['colocation_reactome_enrichment',eachgrp]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['reactome_enrichment']['t']),
				'memory' : compute_resource['reactome_enrichment']['m'],
				'cmd' : 'sh '+os.path.join(mytargetdir,'co_location','_'.join(['colocation_reactome_enrichment',eachgrp])+'.sh')}
			safe_open(os.path.join(mytargetdir,'co_location','_'.join(['colocation_reactome_enrichment',eachgrp])+'.sh'),'w').write(cmd)
			job_points['colocation_go_enrichment'].append('_'.join(['colocation_reactome_enrichment',eachgrp]))

		if 5<=len(list_in_sample.keys())<=55:
			create_dir(os.path.join(mytargetdir,'co_expression'))
			cmd,order1 = aTarget.lnctarget_enrich_coexpression_prepare(species=enrich_spe,diff_soft=diff_soft)
			add_items(orders,order1)
			new_jobs['_'.join(['lnctarget_enrich_coexpression_prepare',eachgrp])] = \
				{'name' : '_'.join(['lnctarget_enrich_coexpression_prepare',eachgrp]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['lnctarget_enrich_coexpression_prepare']['t']),
				'memory' : compute_resource['lnctarget_enrich_coexpression_prepare']['m'],
				'cmd' : 'sh '+os.path.join(mytargetdir,'_'.join(['lnctarget_enrich_coexpression_prepare',eachgrp])+'.sh')}
			safe_open(os.path.join(mytargetdir,'_'.join(['lnctarget_enrich_coexpression_prepare',eachgrp])+'.sh'),'w').write(cmd)
			job_points['lnctarget_enrich_coexpression_prepare'].append('_'.join(['lnctarget_enrich_coexpression_prepare',eachgrp]))
			co_expression_go_dir=os.path.join(mytargetdir,'co_expression','GO')
			co_expression_kegg_dir=os.path.join(mytargetdir,'co_expression','KEGG')
			co_expression_reactome_dir=os.path.join(mytargetdir,'co_expression','Reactome')
			create_dir(co_expression_go_dir)
			create_dir(co_expression_kegg_dir)
			create_dir(os.path.join(co_expression_kegg_dir,'pathway'))
			create_dir(co_expression_reactome_dir)

			cmd = aTarget.go_enrichment(species=enrich_spe,label='co_expression')
			order1 = 'order coexpression_go_enrichment_%s after lnctarget_enrich_coexpression_prepare_%s' % (eachgrp,eachgrp)

			if enrich_spe in ['hsa','mmu']:

				add_items(orders,order1)
				new_jobs['_'.join(['coexpression_go_enrichment',eachgrp])] = \
					{'name' : '_'.join(['coexpression_go_enrichment',eachgrp]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['reactome_enrichment']['t']),
					'memory' : compute_resource['reactome_enrichment']['m'],
					'cmd' : 'sh '+os.path.join(mytargetdir,'co_expression','_'.join(['coexpression_go_enrichment',eachgrp])+'.sh')}
				safe_open(os.path.join(mytargetdir,'co_expression','_'.join(['coexpression_go_enrichment',eachgrp])+'.sh'),'w').write(cmd)
				job_points['coexpression_go_enrichment'].append('_'.join(['coexpression_go_enrichment',eachgrp]))

				cmd = aTarget.kegg_enrichment(species=enrich_spe,label='co_expression')
				order1 = 'order coexpression_kegg_enrichment_%s after lnctarget_enrich_coexpression_prepare_%s' % (eachgrp,eachgrp)
				add_items(orders,order1)
				new_jobs['_'.join(['coexpression_kegg_enrichment',eachgrp])] = \
					{'name' : '_'.join(['coexpression_kegg_enrichment',eachgrp]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['kegg_enrichment']['t']),
					'memory' : compute_resource['kegg_enrichment']['m'],
					'cmd' : 'sh '+os.path.join(mytargetdir,'co_expression','_'.join(['coexpression_kegg_enrichment',eachgrp])+'.sh')}
				safe_open(os.path.join(mytargetdir,'co_expression','_'.join(['coexpression_kegg_enrichment',eachgrp])+'.sh'),'w').write(cmd)
				job_points['coexpression_go_enrichment'].append('_'.join(['coexpression_go_enrichment',eachgrp]))

				cmd = aTarget.reactome_enrichment(species=enrich_spe,label='co_expression')
				order1 = 'order coexpression_reactome_enrichment_%s after lnctarget_enrich_coexpression_prepare_%s' % (eachgrp,eachgrp)
				add_items(orders,order1)
				new_jobs['_'.join(['coexpression_reactome_enrichment',eachgrp])] = \
					{'name' : '_'.join(['coexpression_reactome_enrichment',eachgrp]),
					'status' : 'waiting',
					'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['reactome_enrichment']['t']),
					'memory' : compute_resource['reactome_enrichment']['m'],
					'cmd' : 'sh '+os.path.join(mytargetdir,'co_expression','_'.join(['coexpression_reactome_enrichment',eachgrp])+'.sh')}
				safe_open(os.path.join(mytargetdir,'co_expression','_'.join(['coexpression_reactome_enrichment',eachgrp])+'.sh'),'w').write(cmd)
				job_points['coexpression_go_enrichment'].append('_'.join(['coexpression_reactome_enrichment',eachgrp]))	

	cmd,order1 = aTarget.co_location()
	add_items(orders,order1)
	new_jobs['co_location'] = \
		{'name' : 'co_location',
		'status' : 'waiting',
		'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['co_location']['t']),
		'memory' : compute_resource['co_location']['m'],
		'cmd' : 'sh '+os.path.join(targetdir,'co_location','co_location.sh')}
	safe_open(os.path.join(targetdir,'co_location','co_location.sh'),'w').write(cmd)
	job_points['co_location'].append('co_location')

	if 5<=len(list_in_sample.keys())<=55:
		co_expression_split_dir=os.path.join(targetdir,'co_expression','split')
		create_dir(co_expression_split_dir)
		job_points['co_expression_split']=[]
		cmd,order1 = aTarget.co_expression_split()
		add_items(orders,order1)
		new_jobs['co_expression_split'] = \
			{'name' : 'co_expression_split',
			'status' : 'waiting',
			'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['co_expression_split']['t']),
			'memory' : compute_resource['co_expression_split']['m'],
			'cmd' : 'sh '+os.path.join(targetdir,'co_expression','co_expression_split.sh')}
		safe_open(os.path.join(targetdir,'co_expression','co_expression_split.sh'),'w').write(cmd)
		job_points['co_expression_split'].append('co_expression_split')
	
		for each in range(1,21):
			cmd,order1 = aTarget.co_expression(each)
			add_items(orders,order1)
			new_jobs['co_expression_'+str(each)] = \
				{'name' : 'co_expression_'+str(each),
				'status' : 'waiting',
				'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['co_expression']['t']),
				'memory' : compute_resource['co_expression']['m'],
				'cmd' : 'sh '+os.path.join(targetdir,'co_expression','split','co_expression_'+str(each)+'.sh')}
			safe_open(os.path.join(targetdir,'co_expression','split','co_expression_'+str(each)+'.sh'),'w').write(cmd)
			job_points['co_expression'].append('co_expression')

		job_points['co_expression_merge']=[]
		cmd,order1 = aTarget.co_expression_merge()
		add_items(orders,order1)
		new_jobs['co_expression_merge'] = \
			{'name' : 'co_expression_merge',
			'status' : 'waiting',
			'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['co_expression_merge']['t']),
			'memory' : compute_resource['co_expression_merge']['m'],
			'cmd' : 'sh '+os.path.join(targetdir,'co_expression','co_expression_merge.sh')}
		safe_open(os.path.join(targetdir,'co_expression','co_expression_merge.sh'),'w').write(cmd)
		job_points['co_expression_split'].append('co_expression_merge')

if argv.assembly:
	ass_add='--assembly'
else:
	ass_add=''

if set(['1','10']).issubset(includes):
	print "Circ novel ...\n"
	create_dir(circdir)
	circ_sum_dir  = os.path.join(circdir,'Summary')
	create_dir(circ_sum_dir)
	sampleList = list_in_sample.keys()

    #find_circ
	if circ_soft=='find_circ' or circ_soft=='both':
		job_points['prepare_find_circ'] = []
		job_points['find_circ'] = []
		job_points['find_circ_sum'] = []
		CircNovel = circ_novel(sampleList,'tem',qcdir,circdir,circ_soft,argv.genome)
		cmd,order = CircNovel.prepare_find_circ()
		add_items(orders,order)
		new_jobs['prepare_find_circ'] = \
			{'name':'prepare_find_circ',
			'status' : 'waiting',
			'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['prepare_find_circ']['t']),
			'memory' : compute_resource['prepare_find_circ']['m'],
			'cmd' : 'sh '+os.path.join(circdir,'prepare_find_circ'+'.sh')}
		safe_open(os.path.join(circdir,'prepare_find_circ'+'.sh'),'w').write(cmd)
		job_points['prepare_find_circ'].append('prepare_find_circ')


		for eachsample in sampleList:
			mycircdir = os.path.join(circdir,eachsample,'find_circ')
			os.system("mkdir -p %s" % mycircdir)
			CircNovel = circ_novel(sampleList,eachsample,qcdir,circdir,circ_soft,argv.genome)
			cmd,order1,order2 = CircNovel.find_circ()
			add_items(orders,order1)
			add_items(orders,order2)
			new_jobs['_'.join(['find_circ',eachsample])] = \
				{'name' : '_'.join(['find_circ',eachsample]),
				'status' : 'waiting',
				'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['find_circ']['t']),
				'memory' : compute_resource['find_circ']['m'],
				'cmd' : 'sh '+os.path.join(mycircdir,'_'.join(['find_circ',eachsample])+'.sh')}
			safe_open(os.path.join(mycircdir,'_'.join(['find_circ',eachsample])+'.sh'),'w').write(cmd)
			job_points['find_circ'].append('_'.join(['find_circ',eachsample]))

        #find_circ_sum
		CircNovel = circ_novel(sampleList,'tem',qcdir,circdir,circ_soft,argv.genome)
		cmd = CircNovel.find_circ_summary()		
		new_jobs['find_circ_sum'] = \
			{'name':'find_circ_sum',
			'status' : 'waiting',
			'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['find_circ_sum']['t']),
			'memory' : compute_resource['find_circ_sum']['m'],
			'cmd' : 'sh '+os.path.join(circ_sum_dir,'find_circ_summary.sh')}
		safe_open(os.path.join(circ_sum_dir,'find_circ_summary.sh'),'w').write(cmd)
		job_points['find_circ_sum'].append('find_circ_sum')

    #ciri
	if circ_soft=='ciri' or circ_soft=='both':
		job_points['ciri_bwa'] = []
		job_points['ciri'] = []
		job_points['ciri_sum'] = []

		for eachsample in sampleList:
			mycircdir = os.path.join(circdir,eachsample,'ciri')
			create_dir(mycircdir)
			CircNovel = circ_novel(sampleList,eachsample,qcdir,mycircdir,circ_soft,argv.genome)


            #ciri_bwa
			cmd,order = CircNovel.ciri_bwa()
			add_items(orders,order)
			new_jobs['_'.join(['ciri_bwa',eachsample])] = \
                        	{'name':'_'.join(['ciri_bwa',eachsample]),
                        	'status' : 'waiting',
                        	'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['ciri_bwa']['t']),
                        	'memory' : compute_resource['ciri_bwa']['m'],
                        	'cmd' : 'sh '+os.path.join(mycircdir,'ciri_bwa_'+eachsample+'.sh')}
			safe_open(os.path.join(mycircdir,'ciri_bwa_'+eachsample+'.sh'),'w').write(cmd)
			job_points['ciri_bwa'].append('_'.join(['ciri_bwa_',eachsample]))

           #ciri
			cmd,order1,order2 = CircNovel.ciri()
			add_items(orders,order1)
			add_items(orders,order2)
			new_jobs['_'.join(['ciri',eachsample])] = \
                                {'name' : '_'.join(['ciri',eachsample]),
                                'status' : 'waiting',
                                'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['ciri']['t']),
                                'memory' : compute_resource['ciri']['m'],
                                'cmd' : 'sh '+os.path.join(mycircdir,'_'.join(['ciri',eachsample])+'.sh')}
                        safe_open(os.path.join(mycircdir,'_'.join(['ciri',eachsample])+'.sh'),'w').write(cmd)
                        job_points['ciri'].append('_'.join(['find_circ',eachsample]))

       #ciri_sum
		CircNovel = circ_novel(sampleList,'tem',qcdir,circdir,circ_soft,argv.genome)
		cmd = CircNovel.ciri_summary()
		new_jobs['ciri_sum'] = \
			{'name' : 'ciri_sum',
			'status' : 'waiting',
			'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['ciri_sum']['t']),
			'memory' : compute_resource['ciri_sum']['m'],
			'cmd' : 'sh '+os.path.join(circ_sum_dir,'ciri_summary.sh')}
		safe_open(os.path.join(circ_sum_dir,'ciri_summary.sh'),'w').write(cmd)
		job_points['ciri_sum'].append('ciri_sum')


    #intersection
	if circ_soft == 'both':
		job_points['intersection'] = []
		CircNovel = circ_novel(sampleList,'tem',qcdir,circ_sum_dir,circ_soft,argv.genome)
		order1,order2,cmd = CircNovel.intersection() 
		add_items(orders,order1)
		add_items(orders,order2)
		new_jobs['circ_intersection'] = \
			{'name' : 'circ_intersection',
			'status' : 'waiting',
			'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['intersection']['t']),
			'memory' : compute_resource['intersection']['m'],
			'cmd' : 'sh '+os.path.join(circ_sum_dir,'intersection.sh')}
		safe_open(os.path.join(circ_sum_dir,'intersection.sh'),'w').write(cmd)
		job_points['intersection'].append('intersection')




#annotation
	CircNovel = circ_novel(sampleList,'tem',qcdir,circ_sum_dir,circ_soft,argv.genome)
	cmd,order = CircNovel.circ_annotation()
	job_points['circ_annotation'] = []
	add_items(orders,order)
	new_jobs['circ_annotation'] = \
	{'name':'circ_annotation',
	'status' : 'waiting',
	'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['circ_annotation']['t']),
	'memory' : compute_resource['circ_annotation']['m'],
	'cmd' : 'sh '+os.path.join(circ_sum_dir,'circ_annotation'+'.sh')}
	safe_open(os.path.join(circ_sum_dir,'circ_annotation'+'.sh'),'w').write(cmd)
	job_points['circ_annotation'].append('circ_annotation')

#getseq

	job_points['get_seq'] = []
	circ_seq = os.path.join(circdir,"Sequence")
	create_dir(circ_seq)	
	CircNovel = circ_novel(sampleList,'tem',qcdir,circdir,circ_soft,argv.genome)
	cmd,order = CircNovel.get_seq()
	add_items(orders,order)
	new_jobs['circ_get_seq'] = \
			{'name':'circ_get_seq',
			'status' : 'waiting',
			'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['get_seq']['t']),
			'memory' : compute_resource['get_seq']['m'],
			'cmd' : 'sh '+os.path.join(circ_seq,'get_seq'+'.sh')}
	safe_open(os.path.join(circ_seq,'get_seq.sh'),'w').write(cmd)
	job_points['get_seq'].append('circ_get_seq')

#length

	job_points['length'] = []
	circ_feature_dir = os.path.join(circdir,"Feature")
	circ_F = circ_novel(sampleList,'tem',qcdir,circdir,circ_soft,argv.genome)
	scmd,rcmd,order = circ_F.length()
	out_dir = os.path.join(circdir,"Feature","length")
	os.system('mkdir -p %s'% out_dir)
	safe_open(os.path.join(out_dir,'length.sh'),'w').write(scmd)
	safe_open(os.path.join(out_dir,'length.R'),'w').write(rcmd)
	new_jobs['circ_length'] = \
			{'name':'circ_length',
			'status' : 'waiting',
			'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['length']['t']),
			'memory' : compute_resource['length']['m'],
			'cmd' : 'sh '+os.path.join(out_dir,'length'+'.sh')}
	job_points['length'].append('length')
	add_items(orders,order)

#feature
	job_points['feature'] = ['circ_feature']
	CircNovel = circ_novel(sampleList,'tem',qcdir,circ_feature_dir,circ_soft,argv.genome)
	gfcode,order = CircNovel.genomic_feature()
	out_dir = os.path.join(circdir,"Feature/genomic_feature")
	os.system('mkdir -p %s'% out_dir)
	safe_open(os.path.join(out_dir,'genomic_feature.sh'),'w').write(gfcode)
	new_jobs['circ_feature'] = \
			{'name':'circ_feature',
			'status' : 'waiting',
			'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['feature']['t']),
			'memory' : compute_resource['feature']['m'],               
			'cmd' : 'sh '+os.path.join(out_dir,'genomic_feature.sh')}
	add_items(orders,order)

#circos
	job_points['circos'] = []
	circos_dir = os.path.join(circdir,"Feature/circos")
	create_dir(circos_dir)
	for  eachsample in sampleList:
		sample_circos_dir = os.path.join(circos_dir,eachsample)
		create_dir(sample_circos_dir)
		CircNovel = circ_novel(sampleList,eachsample,qcdir,circdir,circ_soft,argv.genome)
		shcode,order = CircNovel.circos()
		safe_open(os.path.join(sample_circos_dir,'circos.sh'),'w').write(shcode)
		new_jobs[eachsample+'_circ_circos'] = \
			{'name':eachsample+'_circ_circos',
			'status' : 'waiting',
			'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['circos']['t']),
			'memory' : compute_resource['circos']['m'],
			'cmd' : 'sh '+os.path.join(sample_circos_dir,'circos'+'.sh')}
		add_items(orders,order)
		job_points['circos'].append('_'.join(eachsample+'_circos'))

##circ_Target

if set(['1','10','11']).issubset(includes):
	print "CircRna Target ...\n"
	Circ_target = os.path.join(analydir,"CircTarget")
	Circ_target_cf = os.path.join(Circ_target,"ChaiFen")
	create_dir(Circ_target)
	create_dir(Circ_target_cf)
	miRNA_Binding = miRNA_binding(analydir,Circ_target,argv.genome)
	shcode,order,chifen_shell,plot_shell = miRNA_Binding.generate_binding()
	safe_open(os.path.join(Circ_target,'generate_binding.sh'),'w').write(shcode)
	safe_open(os.path.join(Circ_target_cf,'Split_Runtarget.sh'),'w').write(chifen_shell)
	safe_open(os.path.join(Circ_target,'plot_runtarget.sh'),'w').write(plot_shell)
	add_items(orders,order)
	order1 = "order circ_runtarget after circ_generate_binding"
	add_items(orders,order1)
	order2 = 'order plot_runtarget after circ_runtarget'
	add_items(orders,order2)
	new_jobs['circ_generate_binding'] = \
			{'name':'circ_generate_binding',
			'status' : 'waiting',
			'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['generate_binding']['t']),
			'memory' : compute_resource['generate_binding']['m'],
			'cmd' : 'sh '+os.path.join(Circ_target,'generate_binding.sh')}
	job_points['generate_binding'] = ['generate_binding']
	job_points['runtarget'] = ['runtarget']
	new_jobs['circ_runtarget'] = \
			{'name':'circ_runtarget',
			'status' : 'waiting',
			'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['runtarget']['t']),
			'memory' : compute_resource['runtarget']['m'],
			'cmd' : 'sh '+os.path.join(Circ_target,'runtarget.sh')}
	new_jobs['plot_runtarget'] = \
			{'name':'plot_runtarget',
                        'status' : 'waiting',
                        'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['generate_binding']['t']),
                        'memory' : compute_resource['generate_binding']['m'],
                        'cmd' : 'sh '+os.path.join(Circ_target,'plot_runtarget.sh')}

##circ_diff includ 3 steps :circ_exp , generate_diff and DE

if set(['1','10','11','12']).issubset(includes):
	dict_group = {} ;compare_list = []; com =[]
	for line in open(argv.compare_group):
		if line.startswith('#'):
			continue
		array = line.strip().split('\t')
		if array[1] not in dict_group.keys() :
			dict_group[array[1]] = ':'.join(array[2].strip().split(','))
		if array[3] not in dict_group.keys() :
			dict_group[array[3]] = ':'.join(array[4].strip().split(','))
		compare_list.append(array[0])
	groupname = ','.join(dict_group.keys())
	group = ','.join(dict_group.values())
	for i in compare_list :
		I = i.split('_vs_')
		b = dict_group.keys()
		c = str(b.index(I[0])+1) +':'+str(b.index(I[1])+1)
		com.append(c)
	compare = ','.join(com)
	enrich_com = ','.join(compare_list)
		
	print "CircRna diffrence ...\n"
	Circ_dif = os.path.join(analydir,"CircDiff")
	Circ_novel = os.path.join(analydir,"CircNovel")
	create_dir(Circ_dif)
	Circ_enrich = os.path.join(analydir,"CircEnrich")
	create_dir(Circ_enrich)
	Circ_diff = circ_diff(sampleList,compare,group,groupname,Circ_dif,argv.genome,Circ_novel,Circ_enrich)
	shcode_exp,order1 = Circ_diff.circexp() 
	shcode_diff,order2 = Circ_diff.gene_diff()
	order3 = "order circ_DE after circ_generate_diff"
	safe_open(os.path.join(Circ_dif,'circ_exp.sh'),'w').write(shcode_exp)
	safe_open(os.path.join(Circ_dif,'generate_diff.sh'),'w').write(shcode_diff)
	add_items(orders,order1)
	add_items(orders,order2)
	add_items(orders,order3)

	job_points['circ_exp'] = ['circ_exp']
	job_points['circ_generate_diff'] = ['circ_generate_diff']
	job_points['circ_DE'] = ['circ_DE']
	
	new_jobs['circ_exp'] = \
                        {'name':'circ_exp',
                        'status' : 'waiting',
                        'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['circ_exp']['t']),
                        'memory' : compute_resource['circ_exp']['m'],
                        'cmd' : 'sh '+os.path.join(Circ_dif,'circ_exp.sh')}

	new_jobs['circ_generate_diff'] = \
                        {'name':'circ_generate_diff',
                        'status' : 'waiting',
                        'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['generate_diff']['t']),
                        'memory' : compute_resource['generate_diff']['m'],
                        'cmd' : 'sh '+os.path.join(Circ_dif,'generate_diff.sh')}

	new_jobs['circ_DE'] = \
                        {'name':'circ_DE',
                        'status' : 'waiting',
                        'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['DE']['t']),
                        'memory' : compute_resource['DE']['m'],
                        'cmd' : 'sh '+os.path.join(Circ_dif,'DE.sh')}

#Enrichment step3 generate_enrich + run_enrich + eachsample GO and KEGG
	print "CircRna Enrichment ...\n"
	Circ_diff = circ_diff(sampleList,enrich_com,group,groupname,Circ_dif,argv.genome,Circ_novel,Circ_enrich)
	shcode,order1 = Circ_diff.gene_enrich()
	add_items(orders,order1)
	safe_open(os.path.join(Circ_enrich,'circ_generate_enrich.sh'),'w').write(shcode)
	new_jobs['circ_generate_enrich'] = \
                        {'name':'circ_generate_enrich',
                        'status' : 'waiting',
                        'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['generate_enrich']['t']),
                        'memory' : compute_resource['generate_enrich']['m'],
                        'cmd' : 'sh '+os.path.join(Circ_enrich,'circ_generate_enrich.sh')}
	
	order2 = "order circ_run_enrich after circ_generate_enrich"
	add_items(orders,order2)
	new_jobs['circ_run_enrich'] = \
                        {'name':'circ_run_enrich',
                        'status' : 'waiting',
                        'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['run_enrich']['t']),
                        'memory' : compute_resource['run_enrich']['m'],                        'cmd' : 'sh '+os.path.join(Circ_enrich,'run_enrich.sh')}

	for each in compare_list:
		order_GO = "order circ_%s_GO_runenrich after circ_run_enrich" % each
		order_kegg = "order circ_%s_KEGG_runenrich after circ_run_enrich" % each
		add_items(orders,order_GO);add_items(orders,order_kegg)
		GOjob = "circ_"+each+'_GO_runenrich' 
		KEGGjob = "circ_"+each+'_KEGG_runenrich'
		new_jobs[GOjob] = \
                        {'name':GOjob,
                        'status' : 'waiting',
                        'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['GO_runenrich']['t']),
                        'memory' : compute_resource['GO_runenrich']['m'],                        'cmd' : 'sh '+os.path.join(Circ_enrich,each+'_GO_runenrich.sh')}

		new_jobs[KEGGjob] = \
                        {'name':KEGGjob,
                        'status' : 'waiting',
                        'sched' : '-V -cwd %s -l p=%s' % (queue_list,compute_resource['KEGG_runenrich']['t']),
                        'memory' : compute_resource['KEGG_runenrich']['m'],                        'cmd' : 'sh '+os.path.join(Circ_enrich,each+'_KEGG_runenrich.sh')}

##circRNA potential prediction
if set(['10','14']).issubset(includes):
	print "Circ Coding ...\n"
	create_dir(circcoding)
	circpro_dir = circcoding

	IRES=circpro_dir+"/IRES"
	CNCI=circpro_dir+"/ORF/CNCI"
	CPC=circpro_dir+"/ORF/CPC"
	PFAM=circpro_dir+"/ORF/PFAM"
	Stat=circpro_dir+"/Stat"
	os.system("mkdir -p %s %s %s %s %s" %(IRES,CNCI,CPC,PFAM,Stat))

	order1,order2,order3,order4,order5,order6,order7,ires_code,cnci_code,pfam_code,cpc_code,stat_code = circ_Pro(circdir+"/Sequence/novel_circRNA_seq_pro.fa",circpro_dir)
	safe_open(os.path.join(IRES,'circ_ires.sh'),'w').write(ires_code)
	safe_open(os.path.join(CPC,'circ_cpc.sh'),'w').write(cpc_code)
	safe_open(os.path.join(PFAM,'circ_pfam.sh'),'w').write(pfam_code)
	safe_open(os.path.join(CNCI,'circ_cnci.sh'),'w').write(cnci_code)
	safe_open(os.path.join(Stat,'Stat.sh'),'w').write(stat_code)
	map(lambda x: add_items(orders,x), [order1,order2,order3,order4,order5,order6,order7])
	new_jobs['circ_ires'] = {'name' : 'circ_ires',
		'status' : 'waiting',
		'sched' : '-V -cwd %s' % queue_list,
		'memory' : compute_resource['circ_ires']['m'],
		'cmd' : 'sh '+os.path.join(IRES,'circ_ires.sh')}
	new_jobs['circ_cnci'] = {'name' : 'circ_cnci',
		'status' : 'waiting',
                'sched' : '-V -cwd %s' % queue_list,
                'memory' : compute_resource['circ_cnci']['m'],
                'cmd' : 'sh '+os.path.join(CNCI,'circ_cnci.sh')}
	new_jobs['circ_cpc'] = {'name' : 'circ_cpc',
                'status' : 'waiting',
                'sched' : '-V -cwd %s' % queue_list,
                'memory' : compute_resource['circ_cpc']['m'],
                'cmd' : 'sh '+os.path.join(CPC,'circ_cpc.sh')}
	new_jobs['circ_pfam'] = {'name' : 'circ_pfam',
                'status' : 'waiting',
                'sched' : '-V -cwd %s' % queue_list,
                'memory' : compute_resource['circ_pfam']['m'],
                'cmd' : 'sh '+os.path.join(PFAM,'circ_pfam.sh')}
	new_jobs['circ_stat'] = {'name' : 'circ_stat',
                'status' : 'waiting',
                'sched' : '-V -cwd %s' % queue_list,
                'memory' : compute_resource['circ_stat']['m'],
                'cmd' : 'sh '+os.path.join(Stat,'Stat.sh')}

## final Result
create_dir(resultdir)
new_jobs['result'] = {'name' : 'result',
	'status' : 'waiting',
	'sched' : '-V -cwd %s' % queue_list,
	'memory' : compute_resource['result']['m'],
	'cmd' : 'sh '+os.path.join(resultdir,'result.sh')}
shell = safe_open(os.path.join(resultdir,'result.sh'),'w')
create_dir(os.path.join(resultdir,'Results'))
shell.writelines('echo done')
shell.close()

## final Result order
exclude_jobs = ['result']+job_points['pollution']
for eachjob in new_jobs:
	if eachjob in exclude_jobs:
		continue
	order1 = 'order result after %s' % eachjob
	add_items(orders,order1)

## final Report  刘家航修改
create_dir(reportdir)
print "Report and ...\n"

################ 找日期最新的DNA路径-----根据医院和核收日期找同一个人的DNA信息
def grep_with_output(pattern,hospital,heshou_data,projectname,file_pattern):
	DNA_projectname = re.sub(r'-RNA$', '-DNA', projectname)
	try:
		cmd = "grep -w %s %s"%(pattern, file_pattern)
		result = subprocess.check_output(cmd, shell=True, universal_newlines=True, stderr=subprocess.STDOUT)
		### grep结果里挑日期最新的路径		
		if result:
			### 命令执行成功，有输出
			NP_line,DNA_NPpath = '',''
			for outs in result.strip().split('\n'):
				xingming = outs.split("\t")[1]
				yiyuan = outs.split("\t")[6]
				tmp_data = outs.split("\t")[11]
				if xingming == pattern and yiyuan == hospital and tmp_data == heshou_data and DNA_projectname == outs.split("\t")[13]:
					NP_line = outs.split(':')[1]
					DNA_NPpath = outs.split(':')[0]
					
		return NP_line,DNA_NPpath

	except subprocess.CalledProcessError as e:
		### 命令执行失败; 错误输出:e.output   // grep未搜索到相关信息
		#print "DNA+RNA:",pattern,"未找到匹配DNA样本"
		NP_line,DNA_NPpath = '',''
		return NP_line,DNA_NPpath

	except Exception as e:
		return "Error:%s" % str(e)

def grep_with_output_solid(pattern,hospital,heshou_data,projectname,file_pattern):	
	try:
		cmd = "grep -w %s %s"%(pattern, file_pattern)
		result = subprocess.check_output(cmd, shell=True, universal_newlines=True, stderr=subprocess.STDOUT)
		if result:
			NP_line,DNA_NPpath = '',''
			for outs in result.strip().split('\n'):
				xingming = outs.split("\t")[4]
				yiyuan = outs.split("\t")[8]
				tmp_data = outs.split("\t")[13]
				if xingming == pattern and yiyuan == hospital and tmp_data == heshou_data and projectname == outs.split("\t")[22]:
					NP_line = outs.split(':')[1]
					DNA_NPpath = outs.split(':')[0]
		return NP_line,DNA_NPpath

	except subprocess.CalledProcessError as e:
		NP_line,DNA_NPpath = '',''
		return NP_line,DNA_NPpath
	except Exception as e:
		return "Error:%s" % str(e)


#####################################################################
dict_id_path={}
for line in DNARNA_combine.keys():
	DNA_RNA_txt = open('%s/DNA2RNA.%s.txt' % (reportdir,line), 'w')
	DNA_RNA_NP_txt = open('%s/%s.%s.DNA.RNA.NP.txt' % (reportdir,NP.strip().split('/')[-1].split('.NP.txt')[0],line), 'w')
	Line = line.strip()
	print "DNA+RNA:",Line,DNARNA_combine[Line][0]
	##############
	xingming = DNARNA_combine[Line][0]
	hospital = DNARNA_combine[Line][1]
	heshou_data = DNARNA_combine[Line][2]
	projectname = DNARNA_combine[Line][3]
	DNA_sample_name = ''
	if re.search("^RT",Line):
		NP_line,DNA_NPpath = grep_with_output_solid(xingming,hospital,heshou_data,projectname,"/public/project/autohm/workdir/Solid_tumor/DNA/*/*NP.txt")
		if NP_line != "":
			DNA_sample_name = NP_line.split('\t')[2]
	else:
		NP_line,DNA_NPpath = grep_with_output(xingming,hospital,heshou_data,projectname,"/public/project/autohm/workdir/Blood_tumor/DNA/*/*NP.txt")
		DNA_sample_name = NP_line.split('\t')[0]

	dict_id_path[line]=DNA_NPpath
	DNA_RNA_txt.write('%s\t%s\n' % (DNA_sample_name, Line))
	DNA_RNA_NP_txt.write('%s\n' % NP_line)
	DNA_RNA_txt.close()
	DNA_RNA_NP_txt.close()


new_jobs['analysis_report'] = {'name' : 'analysis_report',
	'status' : 'waiting',
	'sched' : '-V -cwd %s' % queue_list,
	'memory' : compute_resource['analysis_report']['m'],
	'cmd' : 'sh '+os.path.join(reportdir,'analysis_report.sh')}
shell = safe_open(os.path.join(reportdir,'analysis_report.sh'),'w')


shell.writelines('export PATH=/DATA01/Software/bin:$PATH\n')
shell.writelines('python %s/script/finalcheck/qc_NEG_sexCheck_RNA.py %s && \\\n' % (root_path, analydir))
###blood
if re.search("Blood_tumor/RNA",analydir):
	shell.writelines('##blood\npython %s/script/Report/BloodPanelReportDNARNA.py -R -p %s -n %s && \\\n' % (root_path, analydir,NP))
	shell.writelines('python %s/script/Report/FPKM/Report_FPKM.py -p %s -n %s && \\\n##DNA+RNA\n' % (root_path, analydir,NP))

	if DNARNA_combine:
		for line in DNARNA_combine.keys():
			DNA_NPpath=dict_id_path[line]
			### 如果DNA路径下能找到对应患者的样本信息，直接执行整合报告脚本信息
			if os.path.exists('/'.join(DNA_NPpath.strip().split('/')[0:-1])):
				shell.writelines('python %s/script/Report/BloodPanelReportDNARNA.py -DR -p %s -n %s/%s.%s.DNA.RNA.NP.txt -r %s -c %s/DNA2RNA.%s.txt && \\\n' % (root_path, '/'.join(DNA_NPpath.strip().split('/')[0:-1]),reportdir,NP.strip().split('/')[-1].split('.NP.txt')[0],line,analydir,reportdir,line))
			### 如果DNA路径下没找到/暂时不存在 对应患者的样本信息，先不执行整合报告脚本信息
			else:
				shell.writelines('#python %s/script/Report/BloodPanelReportDNARNA.py -DR -p %s -n %s/%s.%s.DNA.RNA.NP.txt -r %s -c %s/DNA2RNA.%s.txt && \\\n' % (root_path, '/'.join(DNA_NPpath.strip().split('/')[0:-1]),reportdir,NP.strip().split('/')[-1].split('.NP.txt')[0],line,analydir,reportdir,line))

###solid
if re.search("Solid_tumor/RNA",analydir):
	shell.writelines('##solid\npython %s/script/Report/Solid_Tumor_report.py %s %s && \\\n' % (root_path, analydir,NP))

	if DNARNA_combine:
		shell.writelines('##DNA+RNA solid\n')
		for line in DNARNA_combine.keys():
			DNA_NPpath=dict_id_path[line]
			if os.path.exists('/'.join(DNA_NPpath.strip().split('/')[0:-1])):
				shell.writelines('python %s/script/Report/Solid_Tumor_report.py %s %s/%s.%s.DNA.RNA.NP.txt --dna_path %s --DNA2RNA %s/DNA2RNA.%s.txt && \\\n' % (root_path, analydir, reportdir,NP.strip().split('/')[-1].split('.NP.txt')[0],line, '/'.join(DNA_NPpath.strip().split('/')[0:-1]),reportdir,line ))
			else:
				shell.writelines('#python %s/script/Report/Solid_Tumor_report.py %s %s/%s.%s.DNA.RNA.NP.txt --dna_path %s --DNA2RNA %s/DNA2RNA.%s.txt && \\\n' % (root_path, analydir, reportdir,NP.strip().split('/')[-1].split('.NP.txt')[0],line, '/'.join(DNA_NPpath.strip().split('/')[0:-1]),reportdir,line ))

if set(['8']).issubset(includes):
	shell.writelines('\ncd %s/Results\npython %s/script/Result/result.py --projdir %s --analy_array %s --seqtype %s --qclist %s --NP %s --compare_group %s && \\\n' % (resultdir,root_path,analydir,argv.analy_array,argv.seqtype,os.path.join(analydir,'qc_list'),NP,','.join(this_compare_groups)))
else:
	shell.writelines('\ncd %s/Results\npython %s/script/Result/result.py --projdir %s --analy_array %s --seqtype %s --qclist %s --NP %s && \\\n' % (resultdir,root_path,analydir,argv.analy_array,argv.seqtype,os.path.join(analydir,'qc_list'),NP))
#### send Email
shell.writelines('python %s/script/Email/Email.zip.py -p %s -n %s \n' % (root_path,analydir,NP))

shell.close()

'''
if set(['10']).issubset(includes): #teng
        con,proname,pro,pici = open("%s/pn.txt" % analydir).readline().strip().split("\t")
        shell.write('\n\n'+'\\\n    '.join(['python /PUBLIC/software/CANCER/Pipeline/CircRNA/circ_module/pipeline/circ_report.py',
			'--project_dir %s' % analydir,
			'--project %s' % pro,
			'--project_name %s' % proname,
			'--compare %s' % enrich_com,
			'--sample %s' % ','.join(sampleList),
			'--project_number %s' % pro+"-B"+pici+"-6",
			'--prog %s' % circ_soft]))
'''
shell.close()


## final Report order
exclude_jobs = ['analysis_report']+job_points['pollution']
for eachjob in new_jobs:
	if eachjob in exclude_jobs:
		continue
	order1 = 'order analysis_report after %s' % eachjob
	add_items(orders,order1)

## teng report check
shell_re = safe_open(os.path.join(reportdir,'Report_check.sh'),'w')
if argv.seqtype == "RNAseq":
	datasize = "6"
else:
	datasize = "12"

#### 好像没用
shell_re.write('python /PUBLIC/software/CANCER/Pipeline/Cancer_RNAseq/V3/bin/Report_check.py --analydir %s --array %s --reportdir %s --datesize %s' % (analydir,argv.analy_array,reportdir,datasize))

## qc
if set(['1']).issubset(includes):
	qcreportdir = os.path.join(reportdir,'qc')
	create_dir(qcreportdir)
	for eachjob in new_jobs:
		if eachjob.startswith('qc_'):
			order1 = 'order qc_report after %s' % eachjob
			add_items(orders,order1)
	new_jobs['qc_report'] = {'name' : 'qc_report',
		'status' : 'waiting',
		'sched' : '-V -cwd %s' % queue_list,
		'memory' : compute_resource['qc_report']['m'],
		'cmd' : 'sh '+os.path.join(qcreportdir,'qc_report.sh')}
	shell = safe_open(os.path.join(qcreportdir,'qc_report.sh'),'w')
	job_points['qc_report']=['qc_report'] #teng
	shell.write('echo done')
	shell.close()

## mapping
job_points['mapping_report']=[]
if set(['1','2']).issubset(includes):
	mappingreportdir = os.path.join(reportdir,'mapping')
	create_dir(mappingreportdir)
	for eachjob in new_jobs:
		if eachjob.startswith('mapping_summary_'):
			order1 = 'order mapping_report after %s' % eachjob
			add_items(orders,order1)
		if eachjob.startswith('cal_readcount_'):
			order1 = 'order mapping_report after %s' % eachjob
			add_items(orders,order1)
	job_points['mapping_report'].append('mapping_report')
	new_jobs['mapping_report'] = {'name' : 'mapping_report',
		'status' : 'waiting',
		'sched' : '-V -cwd %s' % queue_list,
		'memory' : compute_resource['mapping_report']['m'],
		'cmd' : 'sh '+os.path.join(mappingreportdir,'mapping_report.sh')}
	shell = safe_open(os.path.join(mappingreportdir,'mapping_report.sh'),'w')
	shell.write('python %s/script/Result/qc_mapping_result.py  --infile %s --analydir %s --outfile %s\n\nrm -rf %s/Mapping/RH001 %s/Mapping/RA101  %s/Mapping/RT102 %s/Mapping/R21071502 && \n\nln -sf %s/RH001 %s/Mapping\nln -sf %s/RA101 %s/Mapping\nln -sf %s/RT102 %s/Mapping\nln -sf %s/R21071502 %s/Mapping\n\n' % (root_path,os.path.join(analydir,'qc_list'),analydir,os.path.join(mappingreportdir,'qc_mapping_summary.xls'),analydir,analydir,analydir,analydir,softwares['Control_bam'],analydir,softwares['Control_bam'],analydir,softwares['Control_bam'],analydir,softwares['Control_bam'],analydir))
	
	shell.close()
	## result after mapping_report
	add_items(orders,'order result after mapping_report')


###  for startpoint
if argv.startpoint:
	print "Startpoint ..."
	print '   ... from %s\n' % ', '.join(startpoints)
	order_relation = {}  
	##  relationship; jobA before jobB  ==>  order_relation[jobA] = jobB;
	for one in orders:
		## order ln_EC1001_N_DHE00213_C4D7FACXX_L6 before qc_EC1001_N_DHE00213_C4D7FACXX_L6
		flag,jobA,relation,jobB = one.strip().split( )
		if relation == 'before':
			if jobA not in order_relation:
				order_relation[jobA] = []
			order_relation[jobA].append(jobB)
		else:
			if jobB not in order_relation:
				order_relation[jobB] = []
			order_relation[jobB].append(jobA)

	## get the succeed jobs in startpoints list
	point_succeeds = set()
	for eachstart in startpoints:
		if eachstart not in job_points:
			print 'POINT '+eachstart+' not in you analysis, Please check your script!\nAvailable POINTS are:\n'+'  \n'.join(job_points.keys())+'\n'
			sys.exit(1)
		for each in job_points[eachstart]:
			tmp_succeeds = get_succeeds(each,order_relation)
			point_succeeds |= set(tmp_succeeds) # 类似append
	for eachjob in new_jobs:
		if eachjob not in point_succeeds:
			new_jobs[eachjob]['status'] = 'done'
	
## jobfile
if True:
#	jobfile = safe_open(os.path.join(logdir,argv.newjob),'w')
	jobfile = safe_open(os.path.join(analydir,argv.newjob),'w')
	jobfile.write('log_dir %s\n' % logdir)
	for eachjob in new_jobs:
		sjmcmd = '\n  '.join(['  name %s' % new_jobs[eachjob]['name'],
#			'memory %s' % new_jobs[eachjob]['memory'],
			'status %s' % new_jobs[eachjob]['status'],
			'sched_options %s -l vf=%s' % (new_jobs[eachjob]['sched'],new_jobs[eachjob]['memory']),
			'cmd_begin',
			'  %s' % new_jobs[eachjob]['cmd'],
			'cmd_end'])
		jobfile.write('job_begin\n'+sjmcmd+'\njob_end\n')
	## orders
	jobfile.write('\n'.join(orders)+'\n')
	jobfile.close()


## Help
help = '''
\033[1;41;43m
1. Please carefully check the paramters and sample_list file. This is VERY important.
2. Submit the job file with "sjm" command: sjm myjobfile.
3. Prepare the "pn.txt" file, with project no., project name and contract no. seperated by tab.
4. If your job is not correctly finished, Please carefully check the log files.
5. After project finished, execute the shell "backup_project.sh" to backup your project information.
6. Donot forget to execute the shell "final_project.sh".
\033[0m
'''
###############
##add the stat info for the project
if os.path.exists(os.path.join(analydir,'pn.txt')):
	cmd_stat = 'python %s --pwd %s --qclist_info %s --pn %s --odir /PROJ/HUMAN/share/PROJ_STAT/Sample_list_stat/' %(os.path.join(softwares['bin'],'proj.stat.py'),analydir, os.path.join(analydir,'qc_list'), os.path.join(analydir,'pn.txt'))
	os.system(cmd_stat)
###############
print "DONE!"
print help
known = os.path.join(analydir,'KNOWN_RULES')
if not os.path.exists(known):
	open(known,'w').write(help)
