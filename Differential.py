import os,sys
import string

class Differential:

	def __init__(self,groupTs,groupNs,groupTname,groupNname,grp_name,alignDir,expDir,diffDir,softwares,databases):
		self.groupTs = groupTs
		self.groupNs = groupNs
		self.groupTname = groupTname
		self.groupNname = groupNname
		self.grp_name = grp_name
		self.alignDir = alignDir
		self.diffRoot = os.path.dirname(diffDir)
		self.Root = os.path.dirname(os.path.dirname(diffDir))
		self.diffDir = diffDir
		self.expDir = expDir
		self.diffName = '%s_vs_%s' % (groupTname,groupNname)
		self.softwares = softwares
		self.databases = databases

	def adjusting_parameter(self,mod='both',FC='2',porq='pval',value='0.05'):
		order = 'order adjusting_parameter after classification_of_RNA_%s' % self.diffName
		if porq not in ['pval','padj']:
			print 'Please enter the correct parameters,pval or padj'
			sys.exit(1)
		if mod=='both':
			step_lis=[]
			step1 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['bin'],'adjusting_parameter.py'),
				'--infile %s' % os.path.join(self.diffDir,self.groupTname+'_vs_'+self.groupNname+'.transcript_level.Differential_analysis_results.xls'),
				'--type %s' % porq,
				'--value %s' % value,
				'--outfile %s' % os.path.join(self.diffDir,self.groupTname+'_vs_'+self.groupNname+'.difftranscript.xls')])
			step_lis.append(step1)
			step2 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['bin'],'adjusting_parameter.py'),
				'--infile %s' % os.path.join(self.diffDir,self.groupTname+'_vs_'+self.groupNname+'.gene_level.Differential_analysis_results.xls'),
				'--type %s' % porq,
				'--value %s' % value,
				'--outfile %s' % os.path.join(self.diffDir,self.groupTname+'_vs_'+self.groupNname+'.diffgene.xls')])
			step_lis.append(step2)
			return ' && \\\n'.join(step_lis),order
		else:
			step1 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['bin'],'adjusting_parameter.py'),
				'--infile %s' % os.path.join(self.diffDir,self.groupTname+'_vs_'+self.groupNname+'.Differential_analysis_results.xls'),
				'--type %s' % porq,
				'--value %s' % value,
				'--outfile %s' % os.path.join(self.diffDir,self.groupTname+'_vs_'+self.groupNname+'.diffgene.xls')])
			return step1,order
		

	def diff_prepare(self,samples_lis,groups,groupname,mod='both',seqtype='RNAseq'):
		order = 'order diff_prepare after quantification_summary'
		cmd_genes_readcount = ' '.join([os.path.join(self.expDir,eachsample,eachsample+'.genes.readcount') for eachsample in samples_lis])
		step_lis=[]
		step1 = 'perl %s %s > %s' % (os.path.join(self.softwares['bin'],'merge_geneID.pl'),cmd_genes_readcount,os.path.join(self.diffRoot,'geneIDList'))
		step_lis.append(step1)
		step2 = '\\\n\t'.join(['perl %s' % os.path.join(self.softwares['bin'],'extractInfo_v2.pl'),
			'-i %s' % os.path.join(self.expDir,'all_transcripts.gtf'),
			'-g %s' % os.path.join(self.diffRoot,'geneIDList'),
			'-o %s' % os.path.join(self.diffRoot,'geneInfo'),
			'-l %s' % os.path.join(self.diffRoot,'genelength')])
		step_lis.append(step2)
		#step2_1 = 'python %s %s %s %s' % (os.path.join(self.softwares['bin'],'gene_description.py'),os.path.join(self.expDir,'all_transcripts.gtf'),os.path.join(self.diffRoot,'trans.description'),os.path.join(self.diffRoot,'gene.description'))
		#step_lis.append(step2_1)
		if mod=='gene' or mod=='both':
			step3 = '\\\n\t'.join(['perl %s' % os.path.join(self.softwares['bin'],'calrowmeans_v2.2.pl'),
				'-rpkm %s' % os.path.join(self.expDir,'Summary','genes.FPKM.xls'),
				'-group %s' % groups,
				'-groupname %s' % groupname,
				'-out-rowmeans %s' % os.path.join(self.expDir,'Summary','genes_rowmeans.FPKM.xls')])
			step_lis.append(step3)
		if mod=='isoform' or mod=='both':
			step4  = '\\\n\t'.join(['perl %s' % os.path.join(self.softwares['bin'],'calrowmeans_v2.2.pl'),
				'-rpkm %s' % os.path.join(self.expDir,'Summary','transcripts.FPKM.xls'),
				'-group %s' % groups,
				'-groupname %s' % groupname,
				'-out-rowmeans %s' % os.path.join(self.expDir,'Summary','transcripts_rowmeans.FPKM.xls')])
			step_lis.append(step4)
			
		return ' && \\\n'.join(step_lis),order

	def edgeR(self,mod='both',lncRNA='no',porq='padj',value='0.05'):
		if porq=='padj':
			type='--padj'
		else:
			type='--pvalue'
		order = 'order diff_expression_%s after diff_prepare' % self.diffName
		step_lis=[]
		if lncRNA!='no':
			if mod=='gene' or mod=='both':
				step1 = '\\\n\t'.join([os.path.join(self.softwares['bin'],'edgeR.R'),
					'--rawcount %s' % os.path.join(self.expDir,'Summary','genes.readcount.xls'),
					'--condition %s' % os.path.join(self.diffRoot,'condition.txt'),
					'--comparename %s' % self.groupTname+'_vs_'+self.groupNname,
					'--mod gene --design normal --foldchange 1 %s %s' % (type,value),
					'--outdir %s ' % self.diffDir])
				step_lis.append(step1)
				step1_1 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['bin'],'anno_diffresult.py'),
					'--infile %s' % os.path.join(self.diffDir,self.groupTname+'_vs_'+self.groupNname+'.diffgene.original_output'),
					'--fpkm %s' % os.path.join(self.expDir,'Summary','genes_rowmeans.FPKM.xls'),
					'--genedes %s' % os.path.join(self.expDir,'gene.description')])
				step_lis.append(step1_1)
				step1_2 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['bin'],'anno_diffresult.py'),
					'--infile %s' % os.path.join(self.diffDir,self.groupTname+'_vs_'+self.groupNname+'.gene_level.Differential_analysis_results.original_output'),
					'--fpkm %s' % os.path.join(self.expDir,'Summary','genes_rowmeans.FPKM.xls'),
					'--genedes %s' % os.path.join(self.expDir,'gene.description')])
				step_lis.append(step1_2)
			if mod=='isoform' or mod=='both':
				step2 = '\\\n\t'.join([os.path.join(self.softwares['bin'],'edgeR.R'),
					'--rawcount %s' % os.path.join(self.expDir,'Summary','transcripts.readcount.xls'),
					'--condition %s' % os.path.join(self.diffRoot,'condition.txt'),
					'--comparename %s' % self.groupTname+'_vs_'+self.groupNname,
					'--mod transcript --design normal --foldchange 1 %s %s' % (type,value),
					'--outdir %s ' % self.diffDir])
				step_lis.append(step2)
				step2_1 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['bin'],'anno_diffresult.py'),
					'--infile %s' % os.path.join(self.diffDir,self.groupTname+'_vs_'+self.groupNname+'.difftranscript.original_output'),
					'--fpkm %s' % os.path.join(self.expDir,'Summary','transcripts_rowmeans.FPKM.xls'),
					'--genedes %s' % os.path.join(self.expDir,'trans.description'),
					'--tr'])
				step_lis.append(step2_1)
				step2_2 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['bin'],'anno_diffresult.py'),
					'--infile %s' % os.path.join(self.diffDir,self.groupTname+'_vs_'+self.groupNname+'.transcript_level.Differential_analysis_results.original_output'),
					'--fpkm %s' % os.path.join(self.expDir,'Summary','transcripts_rowmeans.FPKM.xls'),
					'--genedes %s' % os.path.join(self.expDir,'trans.description'),
					'--tr'])
				step_lis.append(step2_2)
				
		
		else:
			step1 = '\\\n\t'.join([os.path.join(self.softwares['bin'],'edgeR.R'),
				'--rawcount %s' % os.path.join(self.expDir,'Summary','genes.readcount.xls'),
				'--condition %s' % os.path.join(self.diffRoot,'condition.txt'),
				'--comparename %s' % self.groupTname+'_vs_'+self.groupNname,
				'--design normal --foldchange 1 %s %s' % (type,value),
				'--outdir %s ' % self.diffDir])
			step_lis.append(step1)
			step1_1 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['bin'],'anno_diffresult.py'),
				'--infile %s' % os.path.join(self.diffDir,self.groupTname+'_vs_'+self.groupNname+'.diffgene.original_output'),
				'--fpkm %s' % os.path.join(self.expDir,'Summary','genes_rowmeans.FPKM.xls'),
				'--genedes %s' % os.path.join(self.expDir,'gene.description')])
			step_lis.append(step1_1)
			step1_2 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['bin'],'anno_diffresult.py'),
				'--infile %s' % os.path.join(self.diffDir,self.groupTname+'_vs_'+self.groupNname+'.Differential_analysis_results.original_output'),
				'--fpkm %s' % os.path.join(self.expDir,'Summary','genes_rowmeans.FPKM.xls'),
				'--genedes %s' % os.path.join(self.expDir,'gene.description')])
			step_lis.append(step1_2)
		return ' && \\\n'.join(step_lis),order

	#def DESeq(self,geneinfo):
	#	order = 'order diff_expression_%s after diff_prepare' % self.diffName
	#	step1 = '\\\n\t'.join(['perl %s '% os.path.join(self.softwares['bin'],'DESeq.changeFC_v2.2.pl'),
	#		'-i %s' % os.path.join(self.expDir,'Summary','readcount.xls'),
	#		'-a %s -b %s' % (':'.join(self.groupTs),':'.join(self.groupNs)),
	#		'-c %s' % geneinfo,
	#		'-n1 %s -n2 %s -f 1 -p 0.05' % (self.groupTname,self.groupNname),
	#		'-op %s' % self.diffDir])
	#	return step1,order

	def DESeq2(self,geneinfo,porq='padj',value='0.05'):
		order = 'order diff_expression_%s after diff_prepare' % self.diffName
		if porq=='padj':
			type='qval'
		else:
			type='pval'

		step1 = '\\\n\t'.join(['perl %s '% os.path.join(self.softwares['bin'],'DESeq2.pl'),
			'-i %s' % os.path.join(self.expDir,'Summary','readcount.xls'),
			'-a %s -b %s' % (':'.join(self.groupTs),':'.join(self.groupNs)),
			'-n1 %s -n2 %s -f 1 -p %s -ty %s ' % (self.groupTname,self.groupNname,value,type),
			'-op %s' % self.diffDir])
		return step1,order

	def DEGseq(self,geneinfo,porq='padj',value='0.05'):
		order = 'order diff_expression_%s after diff_prepare' % self.diffName
		if porq=='padj':
			type='padj'
		else:
			type='pval'

		step1 = '\\\n\t'.join(['perl %s '% os.path.join(self.softwares['bin'],'DEGseq.pl'),
			'-r %s' % os.path.join(self.expDir,'Summary','readcount.xls'),
			'-a %s -b %s' % (self.groupTs[0],self.groupNs[0]),
			'-s %s,%s' % (self.groupTs[0],self.groupNs[0]),
			'-c %s' % geneinfo,
			'-ty %s' % type,
			'-p %s' % value,
			'-o %s' % self.diffDir])
		return step1,order

	def extract_cuffdiff(self,samples_lis,groups,groupnames,compare,repeat_flag,pvalue='1'):
		if pvalue=='1':
			add_p=''
		else:
			add_p='--pvalue '+pvalue
		if repeat_flag:
			replication='replication'
		else:
			replication='non_replication'
		order = 'order extract_cuffdiff after cuffdiff_quantification'
		step_lis=[]
		step1 = '\\\n\t'.join(['python %s '% os.path.join(self.softwares['bin'],'cuffdiff_extract_'+replication+'.py'),
			'--out-dir %s' % self.diffRoot,
			'--Quantification_dir %s' % os.path.join(self.expDir,'cuffdiff'),
			'--labels %s' % ','.join(samples_lis),
			'--groups %s' % groups,
			'--groupnames %s' % groupnames,
			'--compare %s' % compare,
			'--foldchange 0 --fdr 0.05',
			add_p])
		step_lis.append(step1)
		for each_compare in compare.split(','):
			step_01='mv %s %s.original_output' % (os.path.join(self.diffRoot,each_compare,each_compare+'.diffgene.xls'),os.path.join(self.diffRoot,each_compare,each_compare+'.diffgene'))
			step_lis.append(step_01)
			step_02 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['bin'],'anno_cuffdiff.py'),
				'--infile %s' % os.path.join(os.path.join(self.diffRoot,each_compare,each_compare+'.diffgene.original_output')),
				'--genedes %s' % os.path.join(self.expDir,'gene.description'),
				'--outfile %s' % os.path.join(self.diffRoot,each_compare,each_compare+'.diffgene.xls')])
			step_lis.append(step_02)
			step_03='mv %s %s.original_output' % (os.path.join(self.diffRoot,each_compare,each_compare+'.difftranscript.xls'),os.path.join(self.diffRoot,each_compare,each_compare+'.difftranscript'))
			step_lis.append(step_03)
			step_04 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['bin'],'anno_cuffdiff.py'),
				'--infile %s' % os.path.join(os.path.join(self.diffRoot,each_compare,each_compare+'.difftranscript.original_output')),
				'--genedes %s' % os.path.join(self.expDir,'trans.description'),
				'--outfile %s' % os.path.join(self.diffRoot,each_compare,each_compare+'.difftranscript.xls'),
				'--tr'])
			step_lis.append(step_04)
		return ' && \\\n'.join(step_lis),order

	def classification_of_RNA(self,mod='both',soft='cuffdiff',lncRNA='yes',porq='padj'):
		if soft=='cuffdiff':
			order = 'order classification_of_RNA_%s after extract_cuffdiff' % self.diffName
		else:
			order = 'order classification_of_RNA_%s after diff_expression_%s' % (self.diffName,self.diffName)
		if porq=='padj':
			Volcanoplot=os.path.join(self.softwares['bin'],'Volcanoplot_v2.R')
		else:
			Volcanoplot=os.path.join(self.softwares['bin'],'Volcanoplot_pvalue_v2.R')
		step_list=[]
		if lncRNA=='yes':
			if mod=='gene' or mod=='both':
				step1 = '\\\n\t'.join(['perl %s '% os.path.join(self.softwares['bin'],'extract_info_from_list_v2.pl'),
					'-list %s' % os.path.join(self.expDir,'lncRNA_gene.id'),
					'-table %s' % os.path.join(self.diffDir,self.diffName+'.gene_level.Differential_analysis_results.xls'),
					'-output %s' % os.path.join(self.diffDir,self.diffName+'.lncRNA.gene_level.Differential_analysis_results.xls')])
				step_list.append(step1)
				step2 = '\\\n\t'.join(['perl %s '% os.path.join(self.softwares['bin'],'extract_info_from_list_v2.pl'),
					'-list %s' % os.path.join(self.expDir,'lncRNA_gene.id'),
					'-table %s' % os.path.join(self.diffDir,self.diffName+'.diffgene.xls'),
					'-output %s' % os.path.join(self.diffDir,self.diffName+'.lncRNA.diffgene.xls')])
				step_list.append(step2)
				step3 = '\\\n\t'.join(['perl %s '% os.path.join(self.softwares['bin'],'extract_info_from_list_v2.pl'),
					'-list %s' % os.path.join(self.expDir,'mRNA_gene.id'),
					'-table %s' % os.path.join(self.diffDir,self.diffName+'.gene_level.Differential_analysis_results.xls'),
					'-output %s' % os.path.join(self.diffDir,self.diffName+'.mRNA.gene_level.Differential_analysis_results.xls')])
				step_list.append(step3)
				step4 = '\\\n\t'.join(['perl %s '% os.path.join(self.softwares['bin'],'extract_info_from_list_v2.pl'),
					'-list %s' % os.path.join(self.expDir,'mRNA_gene.id'),
					'-table %s' % os.path.join(self.diffDir,self.diffName+'.diffgene.xls'),
					'-output %s' % os.path.join(self.diffDir,self.diffName+'.mRNA.diffgene.xls')])
				step_list.append(step4)
				step5 = '\\\n\t'.join(['perl %s '% os.path.join(self.softwares['bin'],'extract_info_from_list_v2.pl'),
					'-list %s' % os.path.join(self.expDir,'Unclassified_gene.id'),
					'-table %s' % os.path.join(self.diffDir,self.diffName+'.gene_level.Differential_analysis_results.xls'),
					'-output %s' % os.path.join(self.diffDir,self.diffName+'.Unclassified.gene_level.Differential_analysis_results.xls')])
				step_list.append(step5)
				step6 = '\\\n\t'.join(['perl %s '% os.path.join(self.softwares['bin'],'extract_info_from_list_v2.pl'),
					'-list %s' % os.path.join(self.expDir,'Unclassified_gene.id'),
					'-table %s' % os.path.join(self.diffDir,self.diffName+'.diffgene.xls'),
					'-output %s' % os.path.join(self.diffDir,self.diffName+'.Unclassified.diffgene.xls')])
				step_list.append(step6)
				step_plot1 = '\\\n\t'.join(['/home/weiwenting/anaconda3/envs/R3.3.2/bin/Rscript %s' % Volcanoplot,
					'%s' % os.path.join(self.diffDir,self.diffName+'.lncRNA.gene_level.Differential_analysis_results.xls'),
					'%s' % self.diffDir,
					'%s.lncRNA.diffgene 0.05' % self.diffName])
				step_list.append(step_plot1)
				step_plot2 = '\\\n\t'.join(['/home/weiwenting/anaconda3/envs/R3.3.2/bin/Rscript %s' % Volcanoplot,
					'%s' % os.path.join(self.diffDir,self.diffName+'.mRNA.gene_level.Differential_analysis_results.xls'),
					'%s' % self.diffDir,
					'%s.mRNA.diffgene 0.05' % self.diffName])
				step_list.append(step_plot2)
			if mod=='isoform' or mod=='both':
				step7 = '\\\n\t'.join(['perl %s '% os.path.join(self.softwares['bin'],'extract_info_from_list_v2.pl'),
					'-list %s' % os.path.join(self.expDir,'lncRNA_transcript.id'),
					'-table %s' % os.path.join(self.diffDir,self.diffName+'.transcript_level.Differential_analysis_results.xls'),
					'-output %s' % os.path.join(self.diffDir,self.diffName+'.lncRNA.transcript_level.Differential_analysis_results.xls')])
				step_list.append(step7)
				step8 = '\\\n\t'.join(['perl %s '% os.path.join(self.softwares['bin'],'extract_info_from_list_v2.pl'),
					'-list %s' % os.path.join(self.expDir,'lncRNA_transcript.id'),
					'-table %s' % os.path.join(self.diffDir,self.diffName+'.difftranscript.xls'),
					'-output %s' % os.path.join(self.diffDir,self.diffName+'.lncRNA.difftranscript.xls')])
				step_list.append(step8)
				step9 = '\\\n\t'.join(['perl %s '% os.path.join(self.softwares['bin'],'extract_info_from_list_v2.pl'),
					'-list %s' % os.path.join(self.expDir,'mRNA_transcript.id'),
					'-table %s' % os.path.join(self.diffDir,self.diffName+'.transcript_level.Differential_analysis_results.xls'),
					'-output %s' % os.path.join(self.diffDir,self.diffName+'.mRNA.transcript_level.Differential_analysis_results.xls')])
				step_list.append(step9)
				step10 = '\\\n\t'.join(['perl %s '% os.path.join(self.softwares['bin'],'extract_info_from_list_v2.pl'),
					'-list %s' % os.path.join(self.expDir,'mRNA_transcript.id'),
					'-table %s' % os.path.join(self.diffDir,self.diffName+'.difftranscript.xls'),
					'-output %s' % os.path.join(self.diffDir,self.diffName+'.mRNA.difftranscript.xls')])
				step_list.append(step10)
				step11 = '\\\n\t'.join(['perl %s '% os.path.join(self.softwares['bin'],'extract_info_from_list_v2.pl'),
					'-list %s' % os.path.join(self.expDir,'Unclassified_transcript.id'),
					'-table %s' % os.path.join(self.diffDir,self.diffName+'.transcript_level.Differential_analysis_results.xls'),
					'-output %s' % os.path.join(self.diffDir,self.diffName+'.Unclassified.transcript_level.Differential_analysis_results.xls')])
				step_list.append(step11)
				step12 = '\\\n\t'.join(['perl %s '% os.path.join(self.softwares['bin'],'extract_info_from_list_v2.pl'),
					'-list %s' % os.path.join(self.expDir,'Unclassified_transcript.id'),
					'-table %s' % os.path.join(self.diffDir,self.diffName+'.difftranscript.xls'),
					'-output %s' % os.path.join(self.diffDir,self.diffName+'.Unclassified.difftranscript.xls')])
				step_list.append(step12)
				step_plot3 = '\\\n\t'.join(['/home/weiwenting/anaconda3/envs/R3.3.2/bin/Rscript %s' % Volcanoplot,
					'%s' % os.path.join(self.diffDir,self.diffName+'.lncRNA.transcript_level.Differential_analysis_results.xls'),
					'%s' % self.diffDir,
					'%s.lncRNA.difftranscript 0.05' % self.diffName])
				step_list.append(step_plot3)
				step_plot4 = '\\\n\t'.join(['/home/weiwenting/anaconda3/envs/R3.3.2/bin/Rscript %s' % Volcanoplot,
					'%s' % os.path.join(self.diffDir,self.diffName+'.mRNA.transcript_level.Differential_analysis_results.xls'),
					'%s' % self.diffDir,
					'%s.mRNA.difftranscript 0.05' % self.diffName])
				step_list.append(step_plot4)
			return ' && \\\n'.join(step_list),order
		else:
			step_plot = '\\\n\t'.join(['/home/weiwenting/anaconda3/envs/R3.3.2/bin/Rscript %s' % Volcanoplot,
				'%s' % os.path.join(self.diffDir,self.diffName+'.Differential_analysis_results.xls'),
				'%s' % self.diffDir,
				'%s.diffgene 0.05' % self.diffName])
			return step_plot,order



	def diff_expression(self,geneinfo,soft='edgeR',mod='both',lncRNA='no',porq='padj',value='0.05'):
		#if soft == 'DESeq':
		#	return self.DESeq(geneinfo)
		if soft == 'DESeq2':
			return self.DESeq2(geneinfo)
		elif soft == 'DEGseq':
			return self.DEGseq(geneinfo)
		elif soft == 'edgeR':
			return self.edgeR(mod,lncRNA)

	def go_enrichment(self,species='hsa',lncRNA='no'):
		order = 'order go_enrichment_%s after adjusting_parameter' % (self.diffName)
		if lncRNA == 'yes':
			mod1='.gene_level'
			mod2='.mRNA'
			#order = 'order go_enrichment_%s after adjusting_parameter' % (self.diffName)
		elif lncRNA == 'no':
			mod1=''
			mod2=''
			#order = 'order go_enrichment_%s after diff_expression_%s' % (self.diffName,self.diffName)
		step1 = '\\\n\t'.join(['/home/weiwenting/anaconda3/envs/R4.1/bin/Rscript %s ' % os.path.join(self.softwares['bin'],'GO.R'),
			'--diffgene %s' % os.path.join(self.diffDir,self.diffName+mod2+'.diffgene.xls'),
			'--diffresult %s' % os.path.join(self.diffDir,self.diffName+mod2+mod1+'.Differential_analysis_results.xls'),
			'--species %s --prefix %s' % (species,os.path.join(self.diffDir,'GO',self.diffName))])
		step2 = '\\\n\t'.join(['/home/weiwenting/anaconda3/envs/R4.1/bin/Rscript %s ' % os.path.join(self.softwares['bin'],'EnrichBar.R'),
			'--enrich %s' % os.path.join(self.diffDir,'GO',self.diffName+'.GOenrich.xls'),
			'--nbar 20 --type GO --title %s --prefix %s' % (self.diffName,os.path.join(self.diffDir,'GO',self.diffName))])
		step3 = '\\\n\t'.join(['/home/weiwenting/anaconda3/envs/R4.1/bin/Rscript %s ' % os.path.join(self.softwares['bin'],'EnrichDot.R'),
			'--enrich %s' % os.path.join(self.diffDir,'GO',self.diffName+'.GOenrich.xls'),
			'--ndot 20 --title %s --prefix %s' % (self.diffName,os.path.join(self.diffDir,'GO',self.diffName))])
		return ' && \\\n'.join([step1,step2,step3]),order
	def go_enrichment2(self,species='hsa',lncRNA='no'):
		order = 'order go_enrichment_%s after adjusting_parameter' % (self.diffName)
		if lncRNA == 'yes':
			mod1='.gene_level'
			mod2='.mRNA'
			step1_1='perl %s %s %s %s' %s (self.softares['extractcDNAfromFA.pl'],os.path.join(self.Root,'Assembly','lncRNA_filter','Novel_mRNA.gtf'),self.databases['fasta'],os.path.join(self.diffRoot,'Novel_mRNA_gene.fa'))
			step1_2='perl %s -fa %s -n 1 -out %s' %s (self.softwares['hmm_pfam_go.pl'],os.path.join(self.diffRoot,'Novel_mRNA_gene.fa'),self.diffRoot)
			step1_3='sed -i -e \'s/\t$//g\' %s/Novelgene.hmm_go.txt ' % self.diffRoot
			step1_4='cat %s %s/Novelgene.hmm_go.txt > %s/go.txt' % (self.databases['go'],self.diffRoot,self.diffRoot)
			step1=step1_1+'\n'+step1_2+'\n'+step1_3+'\n'+step1_4
		elif lncRNA == 'no':
			mod1=''
			mod2=''
			#step1='ln -s %s %s/go.txt' % (self.databases['go'],self.diffRoot)
			step1='if [ ! -f %s ];then ln -s %s %s;fi' % (os.path.join(self.diffRoot,'go.txt'),self.databases['go'],os.path.join(self.diffRoot,'go.txt'))
		step2='cut -f1 %s|sed \'1d\' > %s' % (os.path.join(self.diffDir,self.diffName+mod2+'.diffgene.xls'),os.path.join(self.diffDir,'GO',self.diffName+mod2+'.diffgene.id'))
		step3='perl %s -i %s -goann %s -gtf %s -o %s -p %s' % (os.path.join(self.softwares['bin'],'goseq_graph.pl'),os.path.join(self.diffDir,self.diffName+mod2+'.diffgene.id'),os.path.join(self.diffRoot,'go.txt'),self.databases['gtf'],os.path.join(self.diffDir,'GO'),self.diffName) 
		return ' && \\\n'.join([step1,step2,step3]),order	
	def kegg_enrichment(self,species='hsa',lncRNA='no'):
		order = 'order kegg_enrichment_%s after adjusting_parameter' % (self.diffName)
		if lncRNA == 'yes':
			mod1='.gene_level'
			mod2='.mRNA'
			#order = 'order kegg_enrichment_%s after adjusting_parameter' % (self.diffName)
		elif lncRNA == 'no':
			mod1=''
			mod2=''
			#order = 'order kegg_enrichment_%s after diff_expression_%s' % (self.diffName,self.diffName)
		step0 = 'export LD_LIBRARY_PATH=/WORK/software/src/gcc-4.92/gcc/lib64:/WORK/software/src/gcc-4.92/gmp/lib:/WORK/software/src/gcc-4.92/mpc/lib:/WORK/software/src/gcc-4.92/mpfr/lib:LD_LIBRARY_PATH\n'
		step1 = '\\\n\t'.join([os.path.join(self.softwares['bin'],'KEGG.R'),
			'--diffgene %s' % os.path.join(self.diffDir,self.diffName+mod2+'.diffgene.xls'),
			'--diffresult %s' % os.path.join(self.diffDir,self.diffName+mod2+mod1+'.Differential_analysis_results.xls'),
			'--pathway %s' % os.path.join(self.diffDir,'KEGG','pathway'),
			'--species %s --prefix %s' % (species,os.path.join(self.diffDir,'KEGG',self.diffName))])
		step2 = '\\\n\t'.join([os.path.join(self.softwares['bin'],'EnrichBar.R'),
			'--enrich %s' % os.path.join(self.diffDir,'KEGG',self.diffName+'.KEGGenrich.xls'),
			'--nbar 20 --type KEGG --title %s --prefix %s' % (self.diffName,os.path.join(self.diffDir,'KEGG',self.diffName))])
		step3 = '\\\n\t'.join([os.path.join(self.softwares['bin'],'EnrichDot.R'),
			'--enrich %s' % os.path.join(self.diffDir,'KEGG',self.diffName+'.KEGGenrich.xls'),
			'--ndot 20 --title %s --prefix %s' % (self.diffName,os.path.join(self.diffDir,'KEGG',self.diffName))])
		return step0+' && \\\n'.join([step1,step2,step3]),order


	def kegg_enrichment2(self,species='rno',lncRNA='no'):
		order = 'order kegg_enrichment_%s after adjusting_parameter' % (self.diffName)
		if lncRNA=='yes':
			mod1='.gene_level'
			mod2='.mRNA'
			#step0='perl %s %s %s %s/mRNA_gene.fa' % (os.path.join(self.softwares['bin'],'extractcDNAfromFA.pl'),os.path.join(self.Root,'Assembly','lncRNA_filter','Novel_mRNA.gtf'),self.databases['fasta'],self.diffRoot)
			step1='perl %s %s %s %s/mRNA_gene.fa\n/PUBLIC/software/public/Alignment/ncbi-blast-2.2.28+/bin/blastx -query %s/mRNA_gene.fa -db %s -evalue 1e-5 -outfmt 5 -max_target_seqs 1 -num_threads 1 -out %s/kobas_blast.xml' % (os.path.join(self.softwares['bin'],'extractcDNAfromFA.pl'),os.path.join(self.Root,'Assembly','lncRNA_filter','mRNA.gtf'),self.databases['fasta'],self.diffRootself.diffRoot,self.databases['kobas_blast_database'],self.diffRoot)
		else:
			mod1=''
			mod2=''
			#step1='ln -s %s %s/kobas_blast.xml' % (self.databases['kobas_blast_xml'],self.diffRoot)
			step1='if [ ! -f %s ];then ln -s %s %s;fi' % (os.path.join(self.diffRoot,'kobas_blast.xml'),self.databases['kobas_blast_xml'],os.path.join(self.diffRoot,'kobas_blast.xml'))

		step2='cd %s' % os.path.join(self.diffDir,'KEGG')
		step3='cut -f1 %s|sed \'1d\' > %s' % (os.path.join(self.diffDir,self.diffName+mod2+'.diffgene.xls'),os.path.join(self.diffDir,self.diffName+mod2+'.diffgene.id'))
		step4='perl %s -id %s -out-dir %s -species %s -blast-result %s -sample-names %s > %s/KEGG/run.sh' % (os.path.join(self.softwares['bin'],'KEGG_step2_enrich_v2.pl'),os.path.join(self.diffDir,self.diffName+mod2+'.diffgene.id'),os.path.join(self.diffDir,'KEGG'),species,os.path.join(self.diffRoot,'kobas_blast.xml'),self.diffName,self.diffDir)
		step5='sh %s/KEGG/run.sh' % self.diffDir
		step6='python %s --table %s/KEGG/add.%s.identify.xls --abbr %s' % (os.path.join(self.softwares['bin'],'pathway_annotation_flow_parallel_simple_tolerant.py'),self.diffDir,self.diffName,species)
		return ' && \\\n'.join([step1,step2,step3,step4,step5,step6]),order 
	def reactome_enrichment(self,species='hsa',lncRNA='no'):
		order = 'order reactome_enrichment_%s after adjusting_parameter' % (self.diffName)
		if lncRNA == 'yes':
			mod1='.gene_level'
			mod2='.mRNA'
			#order = 'order reactome_enrichment_%s after adjusting_parameter' % (self.diffName)
		elif lncRNA == 'no':
			mod1=''
			mod2=''
			#order = 'order reactome_enrichment_%s after diff_expression_%s' % (self.diffName,self.diffName)
		step1 = '\\\n\t'.join([os.path.join(self.softwares['bin'],'Reactome.R'),
			'--diffgene %s' % os.path.join(self.diffDir,self.diffName+mod2+'.diffgene.xls'),
			'--diffresult %s' % os.path.join(self.diffDir,self.diffName+mod2+mod1+'.Differential_analysis_results.xls'),
			'--species %s --prefix %s' % (species,os.path.join(self.diffDir,'Reactome',self.diffName))])
		step2 = '\\\n\t'.join([os.path.join(self.softwares['bin'],'EnrichBar.R'),
			'--enrich %s' % os.path.join(self.diffDir,'Reactome',self.diffName+'.Reactomeenrich.xls'),
			'--nbar 20 --type Reactome --title %s --prefix %s' % (self.diffName,os.path.join(self.diffDir,'Reactome',self.diffName))])
		step3 = '\\\n\t'.join([os.path.join(self.softwares['bin'],'EnrichDot.R'),
			'--enrich %s' % os.path.join(self.diffDir,'Reactome',self.diffName+'.Reactomeenrich.xls'),
			'--ndot 20 --title %s --prefix %s' % (self.diffName,os.path.join(self.diffDir,'Reactome',self.diffName))])
		return ' && \\\n'.join([step1,step2,step3]),order
	

	def cluster(self,compare,mod='both',lncRNA='yes'):
		order = 'order cluster after adjusting_parameter'
		step_list=[]
		if lncRNA=='yes':
			if mod=='gene' or mod=='both':
				#lncRNA_compare_result_gene=[os.path.join(self.diffRoot,'*','%s.lncRNA.diffgene.xls') % elem for elem in compare.split(',')]
				lncRNA_compare_result_gene=[os.path.join(self.diffRoot,elem,elem+'.lncRNA.diffgene.xls') for elem in compare.split(',')]
				step1 = 'if [ -d %s ];then rm -r %s && mkdir -p %s;else mkdir -p %s;fi' % (os.path.join(self.diffRoot,'cluster','gene','lncRNA'),os.path.join(self.diffRoot,'cluster','gene','lncRNA'),os.path.join(self.diffRoot,'cluster','gene','lncRNA'),os.path.join(self.diffRoot,'cluster','gene','lncRNA'))
				#step1 = 'mkdir -p %s' % os.path.join(self.diffRoot,'cluster','gene','lncRNA')
				step_list.append(step1)
				step2 = '\\\n\t'.join(['python %s '% os.path.join(self.softwares['bin'],'get_cluster.py'),
					'%s' % ','.join(lncRNA_compare_result_gene),
					'%s' % os.path.join(self.expDir,'Summary','anno_genes.FPKM.xls'),
					'%s' % os.path.join(self.diffRoot,'cluster','gene','lncRNA','lncRNA_Diff_gene.FPKM.xls')])
				step_list.append(step2)
				step2_add='less %s|awk -v FS=\'\\t\' -v OFS=\'\\t\' \'{$3="";print $0}\' |sed \'s/\\t\\t/\\t/\'|sed \'s/\\t/:/\'>%s' % (os.path.join(self.diffRoot,'cluster','gene','lncRNA','lncRNA_Diff_gene.FPKM.xls'),os.path.join(self.diffRoot,'cluster','gene','lncRNA','cluster_lncRNA_Diff_gene.FPKM.xls'))
				step_list.append(step2_add)
				step3 = '\\\n\t'.join(['/home/weiwenting/anaconda3/envs/R3.3.2/bin/Rscript %s '% os.path.join(self.softwares['bin'],'cluster.R'),
					'%s' % os.path.join(self.diffRoot,'cluster','gene','lncRNA','cluster_lncRNA_Diff_gene.FPKM.xls'),
					'%s' % os.path.join(self.diffRoot,'cluster','gene','lncRNA')])
				step_list.append(step3)
				step3_add='rm %s' % os.path.join(self.diffRoot,'cluster','gene','lncRNA','cluster_lncRNA_Diff_gene.FPKM.xls')
				step_list.append(step3_add)
	
				#mRNA_compare_result_gene=[os.path.join(self.diffRoot,'*','%s.mRNA.diffgene.xls') % elem for elem in compare.split(',')]
				mRNA_compare_result_gene=[os.path.join(self.diffRoot,elem,elem+'.mRNA.diffgene.xls')  for elem in compare.split(',')]
				step4 = 'if [ -d %s ];then rm -r %s && mkdir -p %s;else mkdir -p %s;fi' % (os.path.join(self.diffRoot,'cluster','gene','mRNA'),os.path.join(self.diffRoot,'cluster','gene','mRNA'),os.path.join(self.diffRoot,'cluster','gene','mRNA'),os.path.join(self.diffRoot,'cluster','gene','mRNA'))
				#step4 = 'mkdir -p %s' % os.path.join(self.diffRoot,'cluster','gene','mRNA')
				step_list.append(step4)
				step5 = '\\\n\t'.join(['python %s '% os.path.join(self.softwares['bin'],'get_cluster.py'),
					'%s' % ','.join(mRNA_compare_result_gene),
					'%s' % os.path.join(self.expDir,'Summary','anno_genes.FPKM.xls'),
					'%s' % os.path.join(self.diffRoot,'cluster','gene','mRNA','mRNA_Diff_gene.FPKM.xls')])
				step_list.append(step5)
				step5_add='less %s|awk -v FS=\'\\t\' -v OFS=\'\\t\' \'{$3="";print $0}\' |sed \'s/\\t\\t/\\t/\'|sed \'s/\\t/:/\'>%s' %(os.path.join(self.diffRoot,'cluster','gene','mRNA','mRNA_Diff_gene.FPKM.xls'),os.path.join(self.diffRoot,'cluster','gene','mRNA','cluster_mRNA_Diff_gene.FPKM.xls'))
				step_list.append(step5_add)
				step6 = '\\\n\t'.join(['/home/weiwenting/anaconda3/envs/R3.3.2/bin/Rscript %s '% os.path.join(self.softwares['bin'],'cluster.R'),
					'%s' % os.path.join(self.diffRoot,'cluster','gene','mRNA','cluster_mRNA_Diff_gene.FPKM.xls'),
					'%s' % os.path.join(self.diffRoot,'cluster','gene','mRNA')])
				step_list.append(step6)
				step6_add='rm %s' % os.path.join(self.diffRoot,'cluster','gene','mRNA','cluster_mRNA_Diff_gene.FPKM.xls')
				step_list.append(step6_add)
	
			if mod=='isoform' or mod=='both':
				#lncRNA_compare_result_transcript=[os.path.join(self.diffRoot,'*','%s.lncRNA.difftranscript.xls') % elem for elem in compare.split(',')]
				lncRNA_compare_result_transcript=[os.path.join(self.diffRoot,elem,elem+'.lncRNA.difftranscript.xls')  for elem in compare.split(',')]
				step7 = 'if [ -d %s ];then rm -r %s && mkdir -p %s;else mkdir -p %s;fi' % (os.path.join(self.diffRoot,'cluster','transcript','lncRNA'),os.path.join(self.diffRoot,'cluster','transcript','lncRNA'),os.path.join(self.diffRoot,'cluster','transcript','lncRNA'),os.path.join(self.diffRoot,'cluster','transcript','lncRNA'))
				#step7 = 'mkdir -p %s' % os.path.join(self.diffRoot,'cluster','transcript','lncRNA')
				step_list.append(step7)
				step8 = '\\\n\t'.join(['python %s '% os.path.join(self.softwares['bin'],'get_cluster.py'),
					'%s' % ','.join(lncRNA_compare_result_transcript),
					'%s' % os.path.join(self.expDir,'Summary','anno_transcripts.FPKM.xls'),
					'%s' % os.path.join(self.diffRoot,'cluster','transcript','lncRNA','lncRNA_Diff_transcript.FPKM.xls')])
				step_list.append(step8)
				step8_add = 'less %s|awk -v FS=\'\\t\' -v OFS=\'\\t\' \'{$3="";$4="";$5="";print $0}\' |sed \'s/\\t\\t\\t\\t/\\t/\'|sed \'s/\\t/:/\'>%s' %(os.path.join(self.diffRoot,'cluster','transcript','lncRNA','lncRNA_Diff_transcript.FPKM.xls'),os.path.join(self.diffRoot,'cluster','transcript','lncRNA','cluster_lncRNA_Diff_transcript.FPKM.xls'))
				step_list.append(step8_add)
				step9 = '\\\n\t'.join(['/home/weiwenting/anaconda3/envs/R3.3.2/bin/Rscript %s '% os.path.join(self.softwares['bin'],'cluster.R'),
					'%s' % os.path.join(self.diffRoot,'cluster','transcript','lncRNA','cluster_lncRNA_Diff_transcript.FPKM.xls'),
					'%s' % os.path.join(self.diffRoot,'cluster','transcript','lncRNA')])
				step_list.append(step9)
				step9_add='rm %s' % os.path.join(self.diffRoot,'cluster','transcript','lncRNA','cluster_lncRNA_Diff_transcript.FPKM.xls')
				step_list.append(step9_add)
				#mRNA_compare_result_transcript=[os.path.join(self.diffRoot,'*','%s.mRNA.difftranscript.xls') % elem for elem in compare.split(',')]
				mRNA_compare_result_transcript=[os.path.join(self.diffRoot,elem,elem+'.mRNA.difftranscript.xls') for elem in compare.split(',')]
				step10 = 'if [ -d %s ];then rm -r %s && mkdir -p %s;else mkdir -p %s;fi' % (os.path.join(self.diffRoot,'cluster','transcript','mRNA'),os.path.join(self.diffRoot,'cluster','transcript','mRNA'),os.path.join(self.diffRoot,'cluster','transcript','mRNA'),os.path.join(self.diffRoot,'cluster','transcript','mRNA'))
				#step10 = 'mkdir -p %s' % os.path.join(self.diffRoot,'cluster','transcript','mRNA')
				step_list.append(step10)
				step11 = '\\\n\t'.join(['python %s '% os.path.join(self.softwares['bin'],'get_cluster.py'),
					'%s' % ','.join(mRNA_compare_result_transcript),
					'%s' % os.path.join(self.expDir,'Summary','anno_transcripts.FPKM.xls'),
					'%s' % os.path.join(self.diffRoot,'cluster','transcript','mRNA','mRNA_Diff_transcript.FPKM.xls')])
				step_list.append(step11)
				step11_add = 'less %s|awk -v FS=\'\\t\' -v OFS=\'\\t\' \'{$3="";$4="";$5="";print $0}\' |sed \'s/\\t\\t\\t\\t/\\t/\'|sed \'s/\\t/:/\'>%s' % (os.path.join(self.diffRoot,'cluster','transcript','mRNA','mRNA_Diff_transcript.FPKM.xls'),os.path.join(self.diffRoot,'cluster','transcript','mRNA','cluster_mRNA_Diff_transcript.FPKM.xls'))
				step_list.append(step11_add)
				step12 = '\\\n\t'.join(['/home/weiwenting/anaconda3/envs/R3.3.2/bin/Rscript %s '% os.path.join(self.softwares['bin'],'cluster.R'),
					'%s' % os.path.join(self.diffRoot,'cluster','transcript','mRNA','cluster_mRNA_Diff_transcript.FPKM.xls'),
					'%s' % os.path.join(self.diffRoot,'cluster','transcript','mRNA')])
				step_list.append(step12)
				step12_add='rm %s' % os.path.join(self.diffRoot,'cluster','transcript','mRNA','cluster_mRNA_Diff_transcript.FPKM.xls')
				step_list.append(step12_add)
			
		else:
			#compare_result=[os.path.join(self.diffRoot,'*','%s.diffgene.xls') % elem for elem in compare.split(',')]
			compare_result=[os.path.join(self.diffRoot,elem,elem+'.diffgene.xls')  for elem in compare.split(',')]
			step13 = 'if [ -d %s ];then rm -r %s && mkdir -p %s;else mkdir -p %s;fi' % (os.path.join(self.diffRoot,'cluster','gene'),os.path.join(self.diffRoot,'cluster','gene'),os.path.join(self.diffRoot,'cluster','gene'),os.path.join(self.diffRoot,'cluster','gene'))
			#step13 = 'mkdir %s' % os.path.join(self.diffRoot,'cluster','gene')
			step_list.append(step13)
			step14 = '\\\n\t'.join(['python %s '% os.path.join(self.softwares['bin'],'get_cluster.py'),
				'%s' % ','.join(compare_result),
				'%s' % os.path.join(self.expDir,'Summary','anno_genes.FPKM.xls'),
				'%s' % os.path.join(self.diffRoot,'cluster','gene','Diff_gene.FPKM.xls')])
			step_list.append(step14)
			step14_add='less %s|awk -v FS=\'\\t\' -v OFS=\'\\t\' \'{$3="";print $0}\' |sed \'s/\\t\\t/\\t/\'|sed \'s/\\t/:/\'>%s' % (os.path.join(self.diffRoot,'cluster','gene','Diff_gene.FPKM.xls'),os.path.join(self.diffRoot,'cluster','gene','cluster_Diff_gene.FPKM.xls'))
			step_list.append(step14_add)
			step15 = '\\\n\t'.join(['/home/weiwenting/anaconda3/envs/R3.3.2/bin/Rscript %s '% os.path.join(self.softwares['bin'],'cluster.R'),
				'%s' % os.path.join(self.diffRoot,'cluster','gene','cluster_Diff_gene.FPKM.xls'),
				'%s' % os.path.join(self.diffRoot,'cluster','gene')])
			step_list.append(step15)
			step15_add = 'rm %s' % os.path.join(self.diffRoot,'cluster','gene','cluster_Diff_gene.FPKM.xls')
			step_list.append(step15_add)
		return ' && \\\n'.join(step_list),order

	def PCA(self,mod='both',lncRNA='yes'):
		order = 'order PCA after quantification_summary'
		step_list=[]
		if lncRNA=='yes':
			if mod=='gene' or mod=='both':
				step1 = 'if [ -d %s ];then rm -r %s && mkdir -p %s;else mkdir -p %s;fi' % (os.path.join(self.diffRoot,'PCA','gene'),os.path.join(self.diffRoot,'PCA','gene'),os.path.join(self.diffRoot,'PCA','gene'),os.path.join(self.diffRoot,'PCA','gene'))
				step_list.append(step1)
				step2 = '\\\n\t'.join(['%s '% os.path.join(self.softwares['bin'],'pca3D.R'),
					'--matrix %s' % os.path.join(self.expDir,'Summary','genes.FPKM.xls'),
					'--condition %s' % os.path.join(self.diffRoot,'condition.txt'),
					'--outdir %s' % os.path.join(self.diffRoot,'PCA','gene')])
				step_list.append(step2)
				step3 = '\\\n\t'.join(['%s '% os.path.join(self.softwares['bin'],'pca2D.R'),
					'--matrix %s' % os.path.join(self.expDir,'Summary','genes.FPKM.xls'),
					'--condition %s' % os.path.join(self.diffRoot,'condition.txt'),
					'--outdir %s' % os.path.join(self.diffRoot,'PCA','gene')])
				step_list.append(step3)
			if mod=='isoform' or mod=='both':
				step1 = 'if [ -d %s ];then rm -r %s && mkdir -p %s;else mkdir -p %s;fi' % (os.path.join(self.diffRoot,'PCA','transcript'),os.path.join(self.diffRoot,'PCA','transcript'),os.path.join(self.diffRoot,'PCA','transcript'),os.path.join(self.diffRoot,'PCA','transcript'))
				step_list.append(step1)
				step2 = '\\\n\t'.join(['%s '% os.path.join(self.softwares['bin'],'pca3D.R'),
					'--matrix %s' % os.path.join(self.expDir,'Summary','transcripts.FPKM.xls'),
					'--condition %s' % os.path.join(self.diffRoot,'condition.txt'),
					'--outdir %s' % os.path.join(self.diffRoot,'PCA','transcript')])
				step_list.append(step2)
				step3 = '\\\n\t'.join(['%s '% os.path.join(self.softwares['bin'],'pca2D.R'),
					'--matrix %s' % os.path.join(self.expDir,'Summary','transcripts.FPKM.xls'),
					'--condition %s' % os.path.join(self.diffRoot,'condition.txt'),
					'--outdir %s' % os.path.join(self.diffRoot,'PCA','transcript')])
				step_list.append(step3)
		else:
			step1 = 'if [ -d %s ];then rm -r %s && mkdir -p %s;else mkdir -p %s;fi' % (os.path.join(self.diffRoot,'PCA'),os.path.join(self.diffRoot,'PCA'),os.path.join(self.diffRoot,'PCA'),os.path.join(self.diffRoot,'PCA'))
			step_list.append(step1)
			step2 = '\\\n\t'.join(['%s '% os.path.join(self.softwares['bin'],'pca3D.R'),
				'--matrix %s' % os.path.join(self.expDir,'Summary','genes.FPKM.xls'),
				'--condition %s' % os.path.join(self.diffRoot,'condition.txt'),
				'--outdir %s' % os.path.join(self.diffRoot,'PCA')])
			step_list.append(step2)
			step3 = '\\\n\t'.join(['%s '% os.path.join(self.softwares['bin'],'pca2D.R'),
				'--matrix %s' % os.path.join(self.expDir,'Summary','genes.FPKM.xls'),
				'--condition %s' % os.path.join(self.diffRoot,'condition.txt'),
				'--outdir %s' % os.path.join(self.diffRoot,'PCA')])
			step_list.append(step3)
		return ' && \\\n'.join(step_list),order
