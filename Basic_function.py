#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os,gzip,re


#########################################################################################
Reference = "/public/database/Reference"


#########################################################################################
pipe_dir = os.path.dirname(os.path.abspath(__file__))
root_path = os.path.dirname(os.path.abspath(pipe_dir))

#################################
def safe_open(file_name,mode='r'):
	try:
		if not file_name.endswith('.gz'):
			return open(file_name,mode)
		else:
			return gzip.open(file_name,mode)
	except IOError:
		print file_name + ' do not exist!'


def add_items(a,b):
	if type(b) == str:
		a.append(b)
	else:
		for mm in b:
			a.append(mm)

def create_dir (dir):
	if not os.path.exists(dir):
		assert not os.system('mkdir -p %s' % dir)

def mv_0(x):
	if int(x) == x:
		return int(x)
	else:
		return x

def update_qclist(samplelist,qclist_file):
	if not os.path.exists(qclist_file):
		assert not os.system("touch "+qclist_file)
	unique_list = [each.split('\t')[-1] for each in open(qclist_file).read().strip().split('\n')]
	assert not os.system('chmod 640 %s'%qclist_file)
	qclist_fh = open(qclist_file,'a')
	for line in open(samplelist):
		if line.startswith('#'):
			continue
		array = line.strip().split('\t')
		assert re.search(u'(\d+)',array[0])
		laneid = re.search(u'(\d+)',array[0]).group(1)
		project_id = os.path.basename(array[5].rstrip('/')).strip('/').split('.')[0]
		fcid = '%s_%s_L' % (project_id, os.path.basename(array[5].rstrip('/')).strip('/').split('_')[-1])
		if ',' in array[5] and array[5].endswith('.gz'):
			fcid = 'L'
		uniq_id = '%s_%s_%s_%s' % (array[2],array[3],array[5],array[0])
		if uniq_id in unique_list:
			continue
		qclist_fh.write('%s%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % \
			(fcid,laneid,array[1],array[2],array[3],array[4],array[5],array[6],uniq_id))
	qclist_fh.close()

## get recursively downstream jobs
def get_succeeds (job,order_relation):
	recursiveSucceeds = [job]
	succeeds = []
	if order_relation.has_key(job):
		succeeds = order_relation[job]
	if len(succeeds) >0:
		for each in succeeds:
			recursiveSucceeds.extend(get_succeeds(each,order_relation))
	return recursiveSucceeds

genome_files = {
	'human_B37':{
		'build':'b37',
		'gtf':'%s/GRCh37.exon.gtf'% Reference,
		'mRNA_gtf':'/PUBLIC/software/CANCER/Database/Genome/human/b37/RNA_reference_data_20171103/mRNA.gtf',
		'lncRNA_gtf':'/PUBLIC/software/CANCER/Database/Genome/human/b37/RNA_reference_data_20171103/lncRNA.gtf',
		'star_index':'/PUBLIC/software/CANCER/Database/Genome/human/b37/RNA_reference_data_20171103/star_index_150',
		'fasta':'%s/hg19.fasta' % Reference,
		'Annotated_lncRNA_information.xls':'/PUBLIC/software/CANCER/Database/Genome/human/b37/RNA_reference_data_20171103/Annotated_lncRNA_information.xls',
		'Annotated_lncRNA_transcript_info.json':'/PUBLIC/software/CANCER/Database/Genome/human/b37/RNA_reference_data_20171103/Annotated_lncRNA_transcript_info.json',
		'Annotated_mRNA_information.xls':'/PUBLIC/software/CANCER/Database/Genome/human/b37/RNA_reference_data_20171103/Annotated_mRNA_information.xls',
		'Annotated_mRNA_transcript_info.json':'/PUBLIC/software/CANCER/Database/Genome/human/b37/RNA_reference_data_20171103/Annotated_mRNA_transcript_info.json',
		'annovar_version':'hg19',
		'annovar_db':'/PUBLIC/software/CANCER/Database/ANNOVAR/humandb',
		'rsem_index_bowtie2':'/PUBLIC/software/CANCER/Database/Genome/human/b37/RNA_reference_data_20171103/rsem-prepare-reference/bowtie2/all_transcripts',
		'rsem_index_star':'/PUBLIC/software/CANCER/Database/Genome/human/b37/RNA_reference_data_20171103/rsem-prepare-reference/star/all_transcripts',
		'starfusion':'%s/GRCh37_gencode_v19_CTAT_lib_Apr032020.plug-n-play/ctat_genome_lib_build_dir'% Reference,
		'pfam_hmm':'/PUBLIC/database/Common/PFAM',
		'rRNA_hsa':'%s/rRNA_hsa/hsa'% Reference,
		'rRNA_E':'/PUBLIC/software/CANCER/Database/RNA_seq/rRNA_other/Eukaryote_rRNA_NR',
		'chrom_bed':'/PUBLIC/software/CANCER/Database/Genome/human/b37/b37.chr25Region.bed',
		'intron_gtf':'%s/GRCh37.intron.gtf'% Reference,
		'soapfuse':'/data2/weiwenting/software/SOAPfuse/Database/'
	},
	'human_B38':{
		'build':'b38',
		'gtf':'%s/hg38/Homo_sapiens.GRCh38.107.gtf'% Reference,
		'mRNA_gtf':'/PUBLIC/software/CANCER/Database/Genome/human/b38/RNA_reference_data_20171103/mRNA.gtf',
		'lncRNA_gtf':'/PUBLIC/software/CANCER/Database/Genome/human/b38/RNA_reference_data_20171103/lncRNA.gtf',
		'star_index':'/PUBLIC/software/CANCER/Database/Genome/human/b38/RNA_reference_data_20171103/star_index_150',
		'fasta':'%s/hg38/Homo_sapiens.GRCh38.dna.toplevel.fa'% Reference,
		'Annotated_lncRNA_information.xls':'/PUBLIC/software/CANCER/Database/Genome/human/b38/RNA_reference_data_20171103/Annotated_lncRNA_information.xls',
		'Annotated_lncRNA_transcript_info.json':'/PUBLIC/software/CANCER/Database/Genome/human/b38/RNA_reference_data_20171103/Annotated_lncRNA_transcript_info.json',
		'Annotated_mRNA_information.xls':'/PUBLIC/software/CANCER/Database/Genome/human/b38/RNA_reference_data_20171103/Annotated_mRNA_information.xls',
		'Annotated_mRNA_transcript_info.json':'/PUBLIC/software/CANCER/Database/Genome/human/b38/RNA_reference_data_20171103/Annotated_mRNA_transcript_info.json',
		'annovar_version':'hg38',
		'annovar_db':'/PUBLIC/software/CANCER/Database/ANNOVAR/humandb/humandb_B38',
		'rsem_index_bowtie2':'/PUBLIC/software/CANCER/Database/Genome/human/b38/RNA_reference_data_20171103/rsem-prepare-reference/bowtie2/all_transcripts',
		'rsem_index_star':'/PUBLIC/software/CANCER/Database/Genome/human/b38/RNA_reference_data_20171103/rsem-prepare-reference/star/all_transcripts',
		'starfusion':'/PUBLIC/software/CANCER/Database/RNA_seq/Fusion/STAR/GRCh38',
		'pfam_hmm':'/PUBLIC/database/Common/PFAM',
		'rRNA_hsa':'%s/rRNA_hsa/hsa' % Reference,
		'rRNA_E':'/PUBLIC/software/CANCER/Database/RNA_seq/rRNA_other/Eukaryote_rRNA_NR',
		'chrom_bed':'/PUBLIC/software/CANCER/Database/Genome/human/b38/B38.chr25Region.bed',
		'intron_gtf':'/PUBLIC/software/CANCER/Database/Genome/human/b38/RNA_reference_data_20171103/GRCh38.intron.gtf',
		'soapfuse':'/PUBLIC/software/CANCER/Database/RNA_seq/Fusion'
	},
	## no 1000indel, dbsnp, cosmic.mutect
	'human_B38_NCCL_spikein':{
		'build':'b38',
		'gtf':'%s/hg38/NCCL_spike_in/Homo_sapiens.GRCh38.107.spikein.gtf' % Reference,
		'mRNA_gtf':'/PUBLIC/software/CANCER/Database/Genome/human/b38/RNA_reference_data_20171103/mRNA.gtf',
		'lncRNA_gtf':'/PUBLIC/software/CANCER/Database/Genome/human/b38/RNA_reference_data_20171103/lncRNA.gtf',
		'star_index':'/PUBLIC/software/CANCER/Database/Genome/human/b38/RNA_reference_data_20171103/star_index_150',
		'fasta':'%s/hg38/NCCL_spike_in/Homo_sapiens.GRCh38.dna.toplevel.spikein.fa' % Reference,
		'Annotated_lncRNA_information.xls':'/PUBLIC/software/CANCER/Database/Genome/human/b38/RNA_reference_data_20171103/Annotated_lncRNA_information.xls',
		'Annotated_lncRNA_transcript_info.json':'/PUBLIC/software/CANCER/Database/Genome/human/b38/RNA_reference_data_20171103/Annotated_lncRNA_transcript_info.json',
		'Annotated_mRNA_information.xls':'/PUBLIC/software/CANCER/Database/Genome/human/b38/RNA_reference_data_20171103/Annotated_mRNA_information.xls',
		'Annotated_mRNA_transcript_info.json':'/PUBLIC/software/CANCER/Database/Genome/human/b38/RNA_reference_data_20171103/Annotated_mRNA_transcript_info.json',
		'annovar_version':'hg38',
		'annovar_db':'/PUBLIC/software/CANCER/Database/ANNOVAR/humandb/humandb_B38',
		'rsem_index_bowtie2':'/PUBLIC/software/CANCER/Database/Genome/human/b38/RNA_reference_data_20171103/rsem-prepare-reference/bowtie2/all_transcripts',
		'rsem_index_star':'/PUBLIC/software/CANCER/Database/Genome/human/b38/RNA_reference_data_20171103/rsem-prepare-reference/star/all_transcripts',
		'starfusion':'/PUBLIC/software/CANCER/Database/RNA_seq/Fusion/STAR/GRCh38',
		'pfam_hmm':'/PUBLIC/database/Common/PFAM',
		'rRNA_hsa':'%s/rRNA_hsa/hsa' % Reference,
		'rRNA_E':'/PUBLIC/software/CANCER/Database/RNA_seq/rRNA_other/Eukaryote_rRNA_NR',
		'chrom_bed':'/PUBLIC/software/CANCER/Database/Genome/human/b38/B38.chr25Region.bed',
		'intron_gtf':'/PUBLIC/software/CANCER/Database/Genome/human/b38/RNA_reference_data_20171103/GRCh38.intron.gtf',
		'soapfuse':'/PUBLIC/software/CANCER/Database/RNA_seq/Fusion'
	},
	'mm10':{
		'build':'mm10',
		'gtf':'/PUBLIC/software/CANCER/Database/Genome/mouse/mm10_RNA_reference_data_20171103/mm10.exon.gtf',
		'mRNA_gtf':'/PUBLIC/software/CANCER/Database/Genome/mouse/mm10_RNA_reference_data_20171103/mRNA.gtf',
		'lncRNA_gtf':'/PUBLIC/software/CANCER/Database/Genome/mouse/mm10_RNA_reference_data_20171103/lncRNA.gtf',
		'star_index':'/PUBLIC/software/CANCER/Database/Genome/mouse/mm10_RNA_reference_data_20171103/star_index_150',
		'fasta':'/PUBLIC/software/CANCER/Database/Genome/mouse/mm10_RNA_reference_data_20171103/mm10.fa',
		'Annotated_lncRNA_information.xls':'/PUBLIC/software/CANCER/Database/Genome/mouse/mm10_RNA_reference_data_20171103/Annotated_lncRNA_information.xls',
		'Annotated_lncRNA_transcript_info.json':'/PUBLIC/software/CANCER/Database/Genome/mouse/mm10_RNA_reference_data_20171103/Annotated_lncRNA_transcript_info.json',
		'Annotated_mRNA_information.xls':'/PUBLIC/software/CANCER/Database/Genome/mouse/mm10_RNA_reference_data_20171103/Annotated_mRNA_information.xls',
		'Annotated_mRNA_transcript_info.json':'/PUBLIC/software/CANCER/Database/Genome/mouse/mm10_RNA_reference_data_20171103/Annotated_mRNA_transcript_info.json',
		'annovar_version':'mm10',
		'annovar_db':'/PUBLIC/software/CANCER/Database/ANNOVAR/mousedb',
		'rsem_index_bowtie2':'/PUBLIC/software/CANCER/Database/Genome/mouse/mm10_RNA_reference_data_20171103/rsem-prepare-reference/bowtie2/all_transcripts',
		'rsem_index_star':'/PUBLIC/software/CANCER/Database/Genome/mouse/mm10_RNA_reference_data_20171103/rsem-prepare-reference/star/all_transcripts',
		'starfusion':'/PUBLIC/software/CANCER/Database/RNA_seq/Fusion/STAR/GRCm38',
		'chrom_bed':'/PUBLIC/software/CANCER/Database/Genome/mouse/mm10_RNA_reference_data_20171103/mm10.chrom.bed',
		'pfam_hmm':'/PUBLIC/database/Common/PFAM',
		'rRNA_hsa':'/PUBLIC/software/CANCER/Database/RNA_seq/rRNA_hsa/hsa',
		'rRNA_E':'/PUBLIC/software/CANCER/Database/RNA_seq/rRNA_other/Eukaryote_rRNA_NR',
		'intron_gtf':'/PUBLIC/software/CANCER/Database/Genome/mouse/mm10_RNA_reference_data_20171103/mm10.intron.gtf',
	},
	'rat':{
		'build':'rn6',
		'gtf':'/PUBLIC/software/CANCER/Database/Genome/rat/rn6_RNA_reference_data_20171206/rn6.exon.gtf',
		'mRNA_gtf':'/PUBLIC/software/CANCER/Database/Genome/rat/rn6_RNA_reference_data_20171206/mRNA.gtf',
		'lncRNA_gtf':'/PUBLIC/software/CANCER/Database/Genome/rat/rn6_RNA_reference_data_20171206/lncRNA.gtf',
		'star_index':'/PUBLIC/software/CANCER/Database/Genome/rat/rn6_RNA_reference_data_20171206/star_index_150',
		'fasta':'/PUBLIC/software/CANCER/Database/Genome/rat/rn6_RNA_reference_data_20171206/rn6.fa',
		'Annotated_lncRNA_information.xs':'/PUBLIC/software/CANCER/Database/Genome/rat/rn6_RNA_reference_data_20171206/Annotated_lncRNA_information.xls',
		'Annotated_lncRNA_transcript_info.json':'/PUBLIC/software/CANCER/Database/Genome/rat/rn6_RNA_reference_data_20171206/Annotated_lncRNA_transcript_info.json',
		'Annotated_mRNA_information.xls':'/PUBLIC/software/CANCER/Database/Genome/rat/rn6_RNA_reference_data_20171206/Annotated_mRNA_information.xls',
		'Annotated_mRNA_transcript_info.json':'/PUBLIC/software/CANCER/Database/Genome/rat/rn6_RNA_reference_data_20171206/Annotated_mRNA_transcript_info.json',
		'annovar_version':'mm10',
		'annovar_db':'/PUBLIC/software/CANCER/Database/ANNOVAR/mousedb',
		'rsem_index_bowtie2':'/PUBLIC/software/CANCER/Database/Genome/rat/rn6_RNA_reference_data_20171206/rsem-prepare-reference/bowtie2/all_transcripts',
		'rsem_index_star':'/PUBLIC/software/CANCER/Database/Genome/rat/rn6_RNA_reference_data_20171206/rsem-prepare-reference/star/all_transcripts',
		'chrom_bed':'/PUBLIC/software/CANCER/Database/Genome/rat/rn6_RNA_reference_data_20171206/rn6.chrom.bed',
		'pfam_hmm':'/PUBLIC/database/Common/PFAM',
		'rRNA_hsa':'/PUBLIC/software/CANCER/Database/RNA_seq/rRNA_hsa/hsa',
		'rRNA_E':'/PUBLIC/software/CANCER/Database/RNA_seq/rRNA_other/Eukaryote_rRNA_NR',
		'intron_gtf':'/PUBLIC/software/CANCER/Database/Genome/rat/rn6_RNA_reference_data_20171206/rn6.intron.gtf',
		'go':'/PUBLIC/software/CANCER/Database/Genome/rat/rn6_RNA_reference_data_20171206/go.txt',
		'kobas_blast_database':'/PUBLIC/software/CANCER/Database/Genome/rat/rn6_RNA_reference_data_20171206/rno.pep.fasta',
		'kobas_blast_xml':'/PUBLIC/software/CANCER/Database/Genome/rat/rn6_RNA_reference_data_20171206/kobas_blast.xml',
        }

}

softwares = {
	#### Control_bam
	'Control_bam':'/public/pipeline/RNA/Control_bam/',
	'fastp':'/public/software/anaconda3/envs/fastp/bin/fastp',
	'qc_stat':'%s/script/QC/get.qc.stat.py' % root_path,
	'script':'%s/script/' % root_path,
	'BED':'%s/bed' % root_path,
	'table_annovar':'%s/script/Mutation/table_annovar.pl' % root_path,
	'R4.0.3':'/public/software/anaconda3/envs/R_403/bin/Rscript',
	### mapping_summary
	'mapping_summary_lib':'PATH=/public/software/anaconda3/envs/R_403/bin:/public/software/anaconda3_new/envs/py27/bin:$PATH',
	'hisat2':'export PATH=/public/software/anaconda3_new/envs/perl_5.26.2/bin:$PATH\n/public/software/anaconda3/envs/hisat2/bin/hisat2',
	'bamdst':'/public/software/bamdst-master/bamdst',
	'java':'/public/software/java/bin/java',
	'picard':'/public/software/tools/picard/CollectInsertSizeMetrics.jar',
	'starfusion':'/public/software/tools/STAR-Fusion/STAR-Fusion',
	'STAR':'/public/software/anaconda3_new/envs/STAR-2.7.11/bin/STAR',
	'STAR-index':'/public/software/anaconda3_new/envs/STAR-2.7.11/STAR-index/',
	'Arriba':'/public/software/anaconda3_new/envs/Arriba/bin/arriba',
	'iHMALLdb':'/public/database/BloodDB/iHMALLdb/',
	'fusioncatcherdb':'/public/software/anaconda3/envs/fusioncatcher/share/fusioncatcher-1.20/db/current/', 
	#### Altersplice 
	'rmats':'%s/script/Altersplice/rmats/RNASeq-MATS.py' % root_path,
	'rMATSresult':'%s/script/Altersplice/rMATSresult.py' % root_path,
	'altersplice_anno':'%s/script/Altersplice/altersplice_anno.py'% root_path,
	'altersplice.filter.blood':'%s/script/Altersplice/altersplice.filter.blood.py' % root_path,
	'altersplice.filter.solid':'%s/script/Altersplice/altersplice.filter.solid.py' % root_path,
	'excel.combine.altersplice':'%s/script/Altersplice/excel.combine.altersplice.pl'  % root_path,
	'rmats2sashimiplot':'%s/script/Altersplice/rmats2sashimiplot' % root_path, 
	###
	'featurecount':'/public/software/tools/subread/subread-2.0.1-Linux-x86_64/bin/featureCounts',
	'perl_library':'export PATH=/public/software/anaconda3_new/envs/perl_5.26.2/bin:$PATH',
	## 
	'ttmv_rara_SR':'/public/software/tools/TTMV_RARA/software/ttmv_rara_SR-main/ttmv_rara_SR.pl',
	'ttmv_rara_results':'%s/script/TTMV_screen/result.py' % root_path,
	'ttmv_rara_excel':'%s/script/TTMV_screen/excel.combine.pl'% root_path,
	'samtools':'/public/software/tools/samtools',
	'bcftools':'/PUBLIC/software/HUMAN/bin/bcftools',
	'bowtie2':'\nexport PATH=/public/software/anaconda3_new/envs/bowtie2/bin:$PATH\nbowtie2',
	## fusion scripts
	'fusion_anno_blood':'%s/script/RNAfusion/Starfusion/fusion_anno.blood.py'% root_path,
	'fusion_filter':'%s/script/RNAfusion/Starfusion/fusion_filter.py' % root_path,
	'filter.result.sort':'%s/script/RNAfusion/filter.result.sort.pl' % root_path,
	'fusioncatcher.py':'/public/software/anaconda3/envs/fusioncatcher/bin/fusioncatcher.py',
	'fusion_filter.fusioncatcher':'%s/script/RNAfusion/FusionCatcher/fusion_filter.fusioncatcher.py' % root_path,
	'fusion_anno.blood.fusioncatcher':'%s/script/RNAfusion/FusionCatcher/fusion_anno.blood.fusioncatcher.py' % root_path,
	'fusion_hotspot':'%s/hotspot/Fusion/Fusion_iHMRNA.Lib.txt'  % root_path,
	'fusion_hotspot_phlike':'%s/hotspot/Fusion/phlike.Lib.txt'  % root_path,
	'fusion_hotspot_solid':'/public/pipeline/Solid_tumor/script_new/hotspot/Fusion_hotspot/hotspot/',
	'excel.combine.fusion':'%s/script/RNAfusion/excel.combine.fusion.pl' % root_path,
	'Fusion_Casevar_path':'%s/database/RNAFusionCasvar' % root_path,
	## fusion FalsePositive
	'False_positive_path':'%s/script/RNAfusion/filter_FalsePositive/' % root_path,
	'false_positive_filter.py':'%s/script/RNAfusion/filter_FalsePositive/filter.false.postive.py' % root_path,
	# arriba 
	'retrieve.fusions':'%s/script/RNAfusion/Arriba/retrieve.fusions.py' % root_path,
	'arriba_hotanno':'%s/script/RNAfusion/Arriba/fusion.hot.anno.py' % root_path,
	'arriba_hotanno_sort':'%s/script/RNAfusion/Arriba/sort_by_reads.py' % root_path,
	'arriba_filter':'%s/script/RNAfusion/Arriba/fusion.filter.py' % root_path,
	'arriba_fusion':'%s/script/RNAfusion/Arriba/filter.fusiongene.py' % root_path,
	'arriba_rearrangement':'%s/script/RNAfusion/Arriba/filter.rearrangement.py' % root_path,
	#'arriba_comb':'%s/script/RNAfusion/Arriba/excel.combine.py' % root_path,
	'arriba_comb':'%s/script/RNAfusion/Arriba/excel.combine.fusion.arriba.pl' % root_path,
	'result_sort':'%s/script/RNAfusion/Arriba/filterResult.sort.pl' % root_path,
	#
	'rearrangement_genelist':'/public/pipeline/Blood_tumor/script_new/Fusion/Lumpy/bin/rearrangement_genelist',
	'solid_rearrangement_genelist':'/public/pipeline/Solid_tumor/script_new/script/Fusion/Lumpy/solid_rearrangement_genelist',
	## annotation 
	#'vcfmerge.pl':'%s/script/Mutation/vcfmerge.test.20191023.pl' % root_path,
	#'Analyses.2017.test2.pl':'%s/script/Mutation/Analyses.2017.test2.pl' % root_path,
	## hot database
	'RNA.selfdb':'%s/hotspot/RNA.selfdb.txt' % root_path,
	'Somatic_Lib_iHMAP':'%s/hotspot/Somatic/iHMAP/iHMAP.Lib.txt' % root_path,
	'Somatic_Lib_iHMLYMP':'%s/hotspot/Somatic/iHMLYMP/iHMLYMP.Lib.txt' % root_path,
	'Somatic_Lib_phlike':'%s/hotspot/Somatic/phlike/phlike.Lib.txt' % root_path,
	'AA_position_iHMAP':'%s/hotspot/Somatic/iHMAP/AA_position.txt' % root_path,
	'AA_position_iHMLYMP':'%s/hotspot/Somatic/iHMLYMP/AA_position.txt' % root_path,
	'AA_position_phlike':'%s/hotspot/Somatic/phlike/AA_position.txt' % root_path, 
	##
	'excel.combine.panel.pl':'%s/script/Mutation/excel.combine.panel.pl'  % root_path,
	'filter_path':'/public/pipeline/Blood_tumor/script/',
	### mutation
	'vardict-java':'/public/software/java/bin/vardict-java',
	'var2vcf_valid.pl':'/public/software/java/bin/var2vcf_valid.pl',
	'teststrandbias.R':'/public/software/java/bin/teststrandbias.R',
	'filter_VarDict_single':'%s/script/Mutation/filter.vcf.VarDict_iHMRNA.py' % root_path,
	### pindel
	'pindel':'/public/software/anaconda3/bin/pindel',
	'pindel2vcf':'/public/software/anaconda3/bin/pindel2vcf',
	'FLT3.KMT2A.bed':'%s/bed/FLT3.KMT2A.bed' % root_path,
	'get_FLT3.v2.AD.py':'%s/script/Pindel/get_FLT3.v2.AD.py' % root_path,
	'pindel_format.py':'%s/script/Pindel/pindel_format.py' % root_path,
	'combine.Pindel':'%s/script/Pindel/excel.combine.Pindel.pl'% root_path,
	### Conpair
	'prepare_conpair':'%s/script/Conpair/prepare.Conpair.py' % root_path,
	'run_conpair':'%s/script/Conpair/Conpair.py'% root_path,
	'conpair_result':'%s/script/Conpair/result.py'% root_path,
	################################################### 以下脚本直接调用 DNA 路径 ############################################################
	## anno v2.1.1
	'vcfmerge.pl':'/public/pipeline/Blood_tumor/script_new/Annotate/vcfmerge.pl',
	'Analyses.2024.pl':'/public/pipeline/Blood_tumor/script_new/Annotate/Analyses.2024.pl',
	## filter v2.1.0
	'iHMAP_XHHM':'/public/pipeline/Blood_tumor/script_new/filter/filter_script/iHMAP_XHHM.filter.pl',
	'iHMLYMP_XHHM':'/public/pipeline/Blood_tumor/script_new/filter/filter_script/iHMLYMP_XHHM.filter.pl',
	'filter.result.sort.iHM303':'/public/pipeline/Blood_tumor/script_new/filter/filter_script/filter.result.sort.iHM303.pl',
	'filter.genelist.false.postive':'/public/pipeline/Blood_tumor/script_new/filter/filter_script/filter.genelist.false.postive_addM.py',
	'anno.transvar':'/public/pipeline/Blood_tumor/script_new/filter/filter_script/anno.transvar.py',
	'result.filter':'%s/script/Mutation/result.filter.py'% root_path,
	## pindel 
	'get_cvt.py':'/public/pipeline/Blood_tumor/script_new/filter/filter_script/get_cvt.py',
	################################################### 以下脚本直接调用 solid 路径 ############################################################
	#'solid_vcfmerge.pl':'/public/pipeline/Solid_tumor/script_new/script/Annotation/vcfmerge.pl',
	#'Analyses.2017.test2.XP784.pl':'/public/pipeline/Solid_tumor/script_new/script/Annotation/Analyses.2017.test2.XP784.pl',
	'XP_784_Somatic_Lib':'/public/pipeline/Solid_tumor/script_new/hotspot/Lib/Somatic',
	'solid_selfdb.txt':'/public/pipeline/Solid_tumor/script_new/hotspot/Lib/selfdb.txt',
	'refData':'/public/database/genome/hg19/hg19.fasta',
	## filter_cmd.sh
	#'tumor_XHHM.XP.filter':'/public/pipeline/Solid_tumor/script_new/script/Annotation/tumor_XHHM.filter.XP.py',
	'tumor_XHHM.XP.filter':'%s/script/Mutation/tumor_XHHM.filter.XP.py'% root_path,
	'anno_transvar_solid':'/public/pipeline/Solid_tumor/script_new/script/Annotation/anno.transvar.py',
	'filter.genelist_solid':'/public/pipeline/Solid_tumor/script_new/script/Annotation/filter.genelist.py',
	'filter.result.sort_solid':'/public/pipeline/Solid_tumor/script_new/script/Annotation/filter.result.sort.pl',
	'solid_excel.combine.panel.pl':'/public/pipeline/Solid_tumor/script_new/script/Annotation/excel.combine.panel.pl',
}
## Other script or path used

compute_resource = {
  'prepare_find_circ' : {'m':'1G' , 't':'1'}, # teng
  'find_circ' : {'m':'8G' , 't':'4'},
  'ciri_bwa' :{'m':'5G' , 't':'4'},
  'ciri' : {'m':'25G' , 't':'2'},
  'find_circ_sum' : {'m':'1G' , 't':'1'},
  'ciri_sum' : {'m':'1G' , 't':'1'},
  'intersection' : {'m':'1G' , 't':'1'},
  'circ_annotation' : {'m':'1G' , 't':'1'},
  'novel_circ_stat' : {'m':'1G' , 't':'1'},
  'get_seq' : {'m':'4G' , 't':'2'},
  'length' : {'m':'1G' , 't':'1'},
  'feature' : {'m':'1G' , 't':'1'},
  'circos' : {'m':'1G' , 't':'1'},
  'generate_binding' : {'m':'2G' , 't':'1'},
  'runtarget' : {'m':'4G' , 't':'2'},
  'circ_exp' : {'m':'1G' , 't':'2'},
  'circ_dplot' : {'m':'2G' , 't':'2'},
  'adj_parameter' : {'m':'2G' , 't':'2'},
  'generate_diff' : {'m':'1G' , 't':'1'},
  'DE' : {'m':'1G' , 't':'1'},
  'generate_enrich' : {'m':'100M' , 't':'1'},
  'run_enrich' : {'m':'100M' , 't':'1'},
  'GO_runenrich' : {'m':'4G' , 't':'1'},
  'KEGG_runenrich' : {'m':'4G' , 't':'1'}, # teng
  'rm_rRNA' : {'m':'1G', 't':'1'},
  'md5_raw' : {'m':'1G', 't':'1'},
  'md5_clean' : {'m':'1G', 't':'1'},
  'qc' : {'m':'1G', 't':'2'},
  'qc_summary' : {'m':'1G', 't':'1'},
  'pollution' : {'m':'1G', 't':'1'},
  'mapping_rnaseq' : {'m':'5G', 't':'1'},
  'mapping_summary' : {'m':'1G', 't':'1'},
  'cal_readcount' : {'m':'1G', 't':'1'},
  'assembly' : {'m':'2G', 't':'1'},
  'merge_assembly' : {'m':'2G', 't':'1'},
  'cuffquant_quantification' : {'m':'3G', 't':'1'},
  'expression_filter' : {'m':'1G', 't':'1'},
  'lncRNA_identification' : {'m':'1G', 't':'1'},
  'lncRNA_signatures' : {'m':'1G', 't':'1'},
  'quant_prepare' : {'m':'1G','t':'1'},
  'rsem_quantification' : {'m':'5G', 't':'1'},
  'quantification' : {'m':'5G', 't':'1'},
  'quantification_FPKM_ratio' : {'m':'1G', 't':'1'},
  'cuffdiff_quantification' : {'m':'5G','t':'1'},
  'prepare_rsem_index' : {'m':'5G', 't':'1'},
  'quantification_summary' : {'m':'1G', 't':'1'},
  'wgcna' :{'m':'10G','t':'2'},
  'extract_fasta' : {'m':'1G', 't':'1'},
  'split_fasta' : {'m':'1G', 't':'1'},
  'lncRNA_cpc' : {'m':'2G', 't':'1'},
  'lncRNA_cnci' : {'m':'1G', 't':'1'},
  'lncRNA_pfam' : {'m':'2G', 't':'1'},
  'merge_cpc' : {'m':'1G','t':'1'},
  'merge_cnci' : {'m':'1G','t':'1'},
  'merge_pfam' : {'m':'1G','t':'1'}, 
  'process_bam' : {'m':'40G', 't':'4'},
  'mutation_calling' : {'m':'40G', 't':'3'},
  'somatic_Pindel' : {'m':'2G', 't':'1'},
  'somatic_Pindel_annovar' : {'m':'2G', 't':'1'},
  'run_conpair' : {'m':'2G', 't':'1'},
  'summary_conpair' : {'m':'1G', 't':'1'},
  'filter_annotation' : {'m':'1G', 't':'1'},
  'fusion_step1' : {'m':'35G', 't':'8'},
  'fusion_step2' : {'m':'1G', 't':'1'},
  'arriba_step1' : {'m':'40G', 't':'5'},
  'arriba_step2' : {'m':'1G', 't':'1'},
  'run_TTMVscreen' : {'m':'1G', 't':'1'},
  'soapfuse_prepare' : {'m':'1G','t':'1'},
  'soapfuse_step1' : {'m':'10G','t':'8'},
  'soapfuse_step2' : {'m':'2G','t':'2'},
  'FusionCatcher_step1':{'m':'1G', 't':'0'},
  'FusionCatcher_step2':{'m':'1G', 't':'1'},
  'diff_prepare' : {'m':'1G', 't':'1'},
  'extract_cuffdiff' : {'m':'1G', 't':'1'},
  'diff_expression' : {'m':'1G', 't':'1'},
  'adjusting_parameter' : {'m':'1G', 't':'1'},
  'classification_of_RNA' : {'m':'1G', 't':'1'},
  'cluster' : {'m':'1G', 't':'1'},
  'go_enrichment' : {'m':'1G', 't':'1'},
  'kegg_enrichment' : {'m':'1G', 't':'1'},
  'reactome_enrichment' : {'m':'1G', 't':'1'},
  'rnaseq_mats' : {'m':'5G', 't':'1'},
  'rnaseq_mats_anno' : {'m':'1G', 't':'1'},
  'co_location' : {'m':'1G', 't':'1'},
  'co_expression_split' : {'m':'1G', 't':'1'},
  'co_expression' : {'m':'1G', 't':'1'},
  'co_expression_merge' : {'m':'1G', 't':'1'},
  'lnctarget_enrich_colocation_prepare' : {'m':'1G', 't':'1'},
  'lnctarget_enrich_coexpression_prepare' : {'m':'1G', 't':'1'},
  'cuffdiff' : {'m':'1G', 't':'1'},
  'circ_ires' : {'m':'1G', 't':'1'},
  'circ_cnci' : {'m':'2G', 't':'1'},
  'circ_cpc' : {'m':'2G', 't':'1'},
  'circ_pfam' : {'m':'2G', 't':'1'},
  'circ_stat' : {'m':'1G', 't':'1'},
  'result' : {'m':'1G', 't':'1'},
  'qc_report' : {'m':'1G', 't':'1'},
  'mapping_report' : {'m':'1G', 't':'1'},
  'analysis_report' : {'m':'1G', 't':'1'},
  'zhushi':{'m':'4G', 't':'1'},
  'filter_cmd':{'m':'1G', 't':'1'},
}
