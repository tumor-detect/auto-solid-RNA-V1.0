#!usr/bin/python
## -*- coding:utf-8 -*-

import os,sys
import string
import re
import time

times=time.strftime('%Y%m%d',time.localtime(time.time()))

class RNAVariants:

	def __init__(self,sampleID,qcDir,alignDir,mutDir,softwares,databases,tumor_type):
		self.sampleID = sampleID
		self.qcDir = qcDir
		self.alignDir = alignDir
		self.mutDir = mutDir
		self.bam = os.path.join(alignDir, self.sampleID+'.bam')
		self.softwares = softwares
		self.annovar_dbs = 'GeneName,refGene,Gencode,avsnp142,cosmic70,clinvar_20150330,gwasCatalog,1000g2014oct_eas,1000g2014oct_all,esp6500siv2_all,exac03_ALL_EAS,ljb26_sift,ljb26_pp2hvar,ljb26_pp2hdiv,ljb26_mt,gerp++gt2'
		self.databases = databases
		self.tumor_type = tumor_type

	def construct_fq(self):
		fq1 = os.path.join(self.qcDir,'_'.join([self.sampleID,'1.clean.fq.gz']))
		fq2 = os.path.join(self.qcDir,'_'.join([self.sampleID,'2.clean.fq.gz']))
		return fq1,fq2

	def process_bam(self):
		order = 'order process_bam_%s after mapping_report' % self.sampleID
		step1 = '\\\n\t'.join(['%s MarkDuplicates'%self.softwares['picard'],
			'-I %s' % self.bam,
			'-M %s' % os.path.join(self.mutDir,self.sampleID+'.marked_dup_metrics.txt'),
			'-O %s' % os.path.join(self.mutDir,self.sampleID+'.nodup.bam')])
		step1_2 = ' \\\n\t'.join(['%s markdup' % os.path.join('sambamba'),
			'--tmpdir %s --overflow-list-size 800000 --remove-duplicates --nthreads=4' % self.mutDir,self.bam, os.path.join(self.mutDir,self.sampleID+'.nodup.bam')])
		step2 = 'samtools index %s\n' % os.path.join(self.mutDir,self.sampleID+'.nodup.bam')
		return ' && \\\n'.join([step1_2,step2]),order

	def mutation_calling(self,TR,soft='VarDict'):
		order = 'order mutation_calling_%s after process_bam_%s' % (self.sampleID,self.sampleID)
		if soft == 'samtools':
			step1 = '\\\n\t'.join(['%s mpileup' % os.path.join(self.softwares['samtools']),
				'-q 1 -C 50 -t DP,SP,DV -m 2 -F 0.002 -ugf',
				self.databases['fasta'], 
				os.path.join(self.mutDir,self.sampleID+'.nodup.bam'),
				'|%s call -vmO v' % os.path.join(self.softwares['bcftools']),
				'-o %s' % os.path.join(self.mutDir,self.sampleID+'.raw.vcf')])
			step2 = 'sed -i \'s#\\t*,#\\t#g\' %s' % os.path.join(self.mutDir,self.sampleID+'.raw.vcf')
			step3 = 'sed -i \'s#,*\\t#\\t#g\' %s' % os.path.join(self.mutDir,self.sampleID+'.raw.vcf')
		elif soft == 'gatk':
			step1 = '\\\n\t'.join(['%s HaplotypeCaller' % self.softwares['gatk'],
				'--java-options \"-Xmx4G\"',
				'-L %s' % TR,
				'-R %s' % self.databases['fasta'],
				'-I %s' % os.path.join(self.mutDir,self.sampleID+'.nodup.bam'),
				'-O %s' % os.path.join(self.mutDir,self.sampleID+'.raw.vcf')])
			step2 = 'sed -i \'s#\\*,##g\' %s' % os.path.join(self.mutDir,self.sampleID+'.raw.vcf')
			step3 = 'sed -i \'s#,\\*##g\' %s' % os.path.join(self.mutDir,self.sampleID+'.raw.vcf')
			step4 = 'python %s -i %s -s %s\n' % (self.softwares['filter_vcf'],os.path.join(self.mutDir,self.sampleID+'.raw.vcf'),os.path.join(self.mutDir,self.sampleID+'_'+self.tumor_type[self.sampleID]+'.filter.somatic.vcf'))
		elif soft == 'VarDict':
			if re.search('RA',self.sampleID):
				vcf_bed = '%s/Vardict_bed/RA.merged.porbe.extend20bp.sort.bed' % self.softwares['BED']
			elif re.search('RH',self.sampleID):
				vcf_bed = '%s/Vardict_bed/RH.merge.add.porbe.sort.merge.bed' % self.softwares['BED']
			elif re.search('RS',self.sampleID):
				vcf_bed = '%s/Vardict_bed/RS.merge.add.porbe.sort.merge.bed ' % self.softwares['BED']
			elif re.search('RT',self.sampleID):
				vcf_bed = '%s/Vardict_bed/RT.merge.add.porbe.sort.merge.bed ' % self.softwares['BED']
			else:
				vcf_bed = '/public/pipeline/Blood_tumor/script_new/BED/AR.merge.add.porbe.sort.merge.bed'
			step1 = 'cd %s\n%s -G %s -f 0.005 -N %s -b %s -c 1 -S 2 -E 3 -g 4 -th 10 --nosv %s |%s %s |%s -N %s -E -f 0.005 > %s.vars.vcf'%(self.mutDir,self.softwares['vardict-java'],self.databases['fasta'],self.sampleID,os.path.join(self.mutDir,self.sampleID+'.nodup.bam'),vcf_bed,self.softwares['R4.0.3'],self.softwares['teststrandbias.R'],self.softwares['var2vcf_valid.pl'],self.sampleID,os.path.join(self.mutDir,self.sampleID))
			step4 = 'python %s -i %s.vars.vcf -s %s_%s.filter.somatic.vcf -a %s\nln -sf %s_%s.filter.somatic.vcf %s_%s.filter.germline.vcf\n'%(self.softwares['filter_VarDict_single'],os.path.join(self.mutDir,self.sampleID),os.path.join(self.mutDir,self.sampleID),self.tumor_type[self.sampleID],'0.005',os.path.join(self.mutDir,self.sampleID),self.tumor_type[self.sampleID],os.path.join(self.mutDir,self.sampleID),self.tumor_type[self.sampleID])
		return ' && \\\n'.join([step1,step4]),order

	####### 20250905 add pindel
	def somatic_Pindel(self,mypindeldir,sampleID,mutDir):
		order = 'order somatic_Pindel_%s after mutation_calling_%s' % (sampleID,sampleID)
		step0 = 'printf \"%s.nodup.bam\t250\t%s\n\">%s/%s.config' % (os.path.join(mutDir,sampleID,sampleID),sampleID,mypindeldir,sampleID)
		step1 = '\ncd %s' % mypindeldir
		step2 = '%s -f %s -i %s/%s.config -o %s -j %s ' % (self.softwares['pindel'],self.databases['fasta'],mypindeldir,sampleID,os.path.join(mypindeldir,sampleID),self.softwares['FLT3.KMT2A.bed'])
		step3 = ' && \\\n'.join([
				'%s -p %s/%s_D -r %s -R hg19 -d %s -v %s/%s_D.vcf' % (self.softwares['pindel2vcf'],mypindeldir,sampleID,self.databases['fasta'],times,mypindeldir,sampleID),
				'%s -p %s/%s_SI -r %s -R hg19 -d %s -v %s/%s_SI.vcf' % (self.softwares['pindel2vcf'],mypindeldir,sampleID,self.databases['fasta'],times,mypindeldir,sampleID),
				'%s -p %s/%s_TD -r %s -R hg19 -d %s -v %s/%s_TD.vcf' % (self.softwares['pindel2vcf'],mypindeldir,sampleID,self.databases['fasta'],times,mypindeldir,sampleID),
				'%s -p %s/%s_INV -r %s -R hg19 -d %s -v %s/%s_INV.vcf' % (self.softwares['pindel2vcf'],mypindeldir,sampleID,self.databases['fasta'],times,mypindeldir,sampleID),
		])
		step4 = 'rm -rf  %s/%s.merge.vcf' %(mypindeldir,sampleID)
		step5 = 'cat %s/%s_SI.vcf |grep \'^#\' > %s/%s.merge.vcf && cat %s/%s_*.vcf |grep -v \'^#\' >> %s/%s.merge.vcf' %(mypindeldir,sampleID,mypindeldir,sampleID,mypindeldir,sampleID,mypindeldir,sampleID)
		step6 = 'python %s %s/%s.merge.vcf\n' % (self.softwares['get_FLT3.v2.AD.py'],mypindeldir,sampleID)
		return ' && \\\n'.join([step0,step1,step2,step3,step4,step5,step6]),order


	def somatic_Pindel_annovar(self,mypindeldir,sampleID):
		order = 'order somatic_Pindel_annovar_%s after somatic_Pindel_%s ' % (sampleID,sampleID)
		step1 = 'export PATH=/public/software/anaconda3_new/envs/perl_5.26.2/bin:/public/software/anaconda3_new/envs/perl_5.26.2/bin:/public/software/anaconda3/envs/python27/bin/:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin\ncd %s' % mypindeldir 
		step2 = 'python %s %s/%s.merge.filter.vcf %s/%s.vcf.cvt' % (self.softwares['get_cvt.py'],mypindeldir,sampleID,mypindeldir,sampleID)
		step3 = '###\ncvtfile_size=%s/%s.vcf.cvt' % (mypindeldir,sampleID)
		step4 =  ''.join([
			'if [ -s $cvtfile_size ];then\n',
			'perl %s --step step1 --input %s/%s.merge.filter.vcf --db iHMAPdb\n'  % (self.softwares['Analyses.2024.pl'],mypindeldir,sampleID),
			'python %s -i %s/%s.merge.filter.vcf.anno.txt -o %s/%s.ITD.KMT2A.xls\n' % (self.softwares['pindel_format.py'],mypindeldir,sampleID,mypindeldir,sampleID),
			'else\nls -s %s/%s.vcf.cvt |awk \'{print $10 "阴性"}\' > %s/%s.ITD.KMT2A.xls\nfi' % (mypindeldir,sampleID,mypindeldir,sampleID),
			])
		step5 = 'perl %s %s/%s.ITD.KMT2A.xls %s/%s.ITD.KMT2A.xlsx\n' % (self.softwares['combine.Pindel'],mypindeldir,sampleID,mypindeldir,sampleID)
		return ' && \\\n'.join([step1,step2,step3,step4,step5]),order

    ### 20250908 add 
	def run_conpair(self,myConpair,NP,eachsample,mutdir):
		order = 'order run_conpair_%s after mutation_calling_%s' % (eachsample,eachsample)
		step1 = 'cd %s' % myConpair
		step2 = '##prepare\npython %s -n %s -rb %s -s %s -o %s/samplelist_standard_%s' % (self.softwares['prepare_conpair'],NP,os.path.join(mutdir,eachsample,eachsample+".nodup.bam"),eachsample,myConpair,eachsample)
		step3 = '##run conpair\npython %s -i %s/samplelist_standard_%s -o %s' % (self.softwares['run_conpair'],myConpair,eachsample,myConpair)
		step4 = '##result\nif [ -s %s/concordance.pairwise.xls ];then\npython %s %s/concordance.pairwise.xls > %s/concordance.result\nfi\n' % (myConpair,self.softwares['conpair_result'],myConpair,myConpair)
		return ' && \\\n'.join([step1,step2,step3,step4]),order

	#def summary_conpair(self,conpairdir,myConpair,eachsample):
	#	order = 'order summary_conpair after run_conpair_%s' % eachsample
	#	step1 = 'cd %s' % conpairdir
	#	step2 = 'cat %s/*/concordance.pairwise.xls > conpair.summaryResult.out' % (conpairdir)
	#	return ' && \\\n'.join([step1,step2]),order

	### 20251111 TTMV_screen
	def TTMVscreen(self,eachsample,mymapdir,TTMVscreendir):
		order = 'order run_TTMVscreen after mapping_report' 
		cmd = '### %s\n' % eachsample
		cmd += 'perl %s -i %s/%s.bam -o %s -g hg19 --taxtype ttmv &&\\\n' % (self.softwares['ttmv_rara_SR'],mymapdir,eachsample,TTMVscreendir)
		cmd += 'python3 %s -i %s_ttmv.out -o %s-ttmv.out.xls &&\\\n' % (self.softwares['ttmv_rara_results'],eachsample,eachsample) 
		cmd += 'perl %s %s-ttmv.out.xls %s-ttmv.out.xlsx \n' % (self.softwares['ttmv_rara_excel'],eachsample,eachsample)
		return cmd,order 


	def filter_annotation(self,soft='samtools',genome='human_B37'):
		order = 'order filter_annotation_%s after mutation_calling_%s' % (self.sampleID,self.sampleID)
		if soft == 'samtools':
			step1 = '\\\n\t'.join(['%s filter' % os.path.join(self.softwares['bcftools']),
				'-s FLTER -i "%QUAL>20 && DP>4 && MQ>30"',
				os.path.join(self.mutDir,self.sampleID+'.raw.vcf'),
				'> %s' % os.path.join(self.mutDir,self.sampleID+'.filter.vcf')])
			step2 = 'awk -F "\\t" \'{if(/^#/){print}else{if($7=="PASS" && !/INDEL;/){print}}}\' \\\n\t%s \\\n\t> %s ' % \
				(os.path.join(self.mutDir,self.sampleID+'.filter.vcf'),os.path.join(self.mutDir,self.sampleID+'.snp.vcf'))
			step3 = 'awk -F "\\t" \'{if(/^#/){print}else{if($7=="PASS" && /INDEL;/){print}}}\' \\\n\t%s \\\n\t> %s ' % \
				(os.path.join(self.mutDir,self.sampleID+'.filter.vcf'),os.path.join(self.mutDir,self.sampleID+'.indel.vcf'))
		else:
			step1 = '\\\n\t'.join([self.softwares['java'],
				'-Xmx1g -jar %s' % self.softwares['GATK'],
				'-T VariantFiltration -window 35 -cluster 3 ',
				'-filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0"',
				'-R %s' % self.databases['fasta'],
				'-V %s' % os.path.join(self.mutDir,self.sampleID+'.raw.vcf'),
				'-o %s' % os.path.join(self.mutDir,self.sampleID+'.filter.vcf')])
			step2 = '\\\n\t'.join([self.softwares['java'],
				'-Xmx1g -jar %s' % self.softwares['GATK'],
				'-T SelectVariants -selectType SNP',
				'-R %s' % self.databases['fasta'],
				'-V %s' % os.path.join(self.mutDir,self.sampleID+'.filter.vcf'),
				'-o %s' % os.path.join(self.mutDir,self.sampleID+'.snp.vcf')])
			step3 = '\\\n\t'.join([self.softwares['java'],
				'-Xmx1g -jar %s' % self.softwares['GATK'],
				'-T SelectVariants -selectType INDEL',
				'-R %s' % self.databases['fasta'],
				'-V %s' % os.path.join(self.mutDir,self.sampleID+'.filter.vcf'),
				'-o %s' % os.path.join(self.mutDir,self.sampleID+'.indel.vcf')])
		z=''
		annovar_db_in=''
		if genome in ['human_B37','human_B38']:
			z=self.annovar_dbs
		elif genome in ['mm10']:
			z='GeneName,refGene,genomicSuperDups,gff3,snp142chr'
		step4 = '\\\n\t'.join([self.softwares['annovar'],
			'-b SNP -u %s' % (self.databases['annovar_version']),
			'-r %s' % self.databases['fasta'],
			'-d %s' % self.databases['annovar_db'],
			'-z %s' % z,
			'%s %s' % (os.path.join(self.mutDir,self.sampleID+'.snp.vcf'),self.sampleID)])
		step5 = '\\\n\t'.join([self.softwares['annovar'],
			'-b INDEL -u %s' % (self.databases['annovar_version']),
			'-r %s' % self.databases['fasta'],
			'-d %s' % self.databases['annovar_db'],
			'-z %s' % z,
			'%s %s' % (os.path.join(self.mutDir,self.sampleID+'.indel.vcf'),self.sampleID)])
		step6 = 'rm %s' % os.path.join(self.mutDir,self.sampleID+'.filter.vcf*')
		step7 = 'gzip -f %s' % os.path.join(self.mutDir,self.sampleID+'.raw.vcf')
		step8 = 'bgzip -f %s && \\\nbgzip -f %s' % (os.path.join(self.mutDir,self.sampleID+'.snp.vcf'),os.path.join(self.mutDir,self.sampleID+'.indel.vcf'))
		step9 = 'tabix -p vcf %s.gz && \\\ntabix -p vcf %s.gz' % (os.path.join(self.mutDir,self.sampleID+'.snp.vcf'), \
			os.path.join(self.mutDir,self.sampleID+'.indel.vcf'))
		step10 = 'rm %s' % os.path.join(self.mutDir,self.sampleID+'.nodup.bam')
		step10_2 = 'rm %s' % os.path.join(self.mutDir,self.sampleID+'.nodup.bai')
		if z=='':
			return ' && \\\n'.join([step1,step2,step3,step6,step7,step8,step9,step10,step10_2]),order
		else:
			return ' && \\\n'.join([step1,step2,step3,step4,step5,step6,step7,step8,step9,step10,step10_2]),order

	def zhushi(self,zhushidir,analydir,t_type,i,TR):
		order=[]
		for each in t_type[i]:
			order.append('order zhushi_%s after mutation_calling_%s' %(i,each))
		sample_tmp = []
		dirtmp = analydir+'/Mutation'
		vcf_file=open("%s/%s.vcf.conf"%(zhushidir,i),'w')
		sample_file=open("%s/%s.sample.conf"%(zhushidir,i),'w')
		for each in t_type[i]:
			if not re.search("NEG",each):
				vcf_file.write("%s/%s/%s_%s.filter.somatic.vcf\t%s\n"%(dirtmp,each,each,i,each))
				sample_file.write("case\t%s\t%s/%s/%s_%s.filter.somatic.vcf\n"%(each,dirtmp,each,each,i))
			else:
				vcf_file.write("%s/%s/%s.vars.vcf\t%s\n"%(dirtmp,each,each,each))
				sample_file.write("case\t%s\t%s/%s/%s.vars.vcf\n"%(each,dirtmp,each,each))
		vcf_file.close()
		sample_file.close()
		################################
		step1 = ' && \\\n'.join(['%s\ncd %s\nperl %s %s.vcf.conf %s.vcf %s'%(self.softwares['perl_library'],zhushidir,self.softwares['vcfmerge.pl'],i,i,self.databases['fasta']),
			'sed -i \'s#1/0#0/1#g\' %s.vcf'%(i),
			'perl %s --step step1 --input %s.vcf --db iHMALLdb'%(self.softwares['Analyses.2024.pl'],i)])
		step2_1 = 'perl %s --step step2 --Func splicing,exonic --input %s.vcf --printall true --amplicon %s --yourlib %s > %s.vcf.focus.xls 2> %s.vcf.ignore.xls'%(self.softwares['Analyses.2024.pl'],i,TR,self.softwares['Somatic_Lib_iHMAP'],i,i)
		step2_2 = 'perl %s --step step2 --Func splicing,exonic --input %s.vcf --printall true --amplicon %s --yourlib %s > %s.vcf.focus.xls 2> %s.vcf.ignore.xls'%(self.softwares['Analyses.2024.pl'],i,TR,self.softwares['Somatic_Lib_phlike'],i,i)
		step2_3 = 'perl %s --step step2 --Func splicing,exonic --input %s.vcf --printall true --amplicon %s --yourlib %s > %s.vcf.focus.xls 2> %s.vcf.ignore.xls'%(self.softwares['Analyses.2024.pl'],i,TR,self.softwares['Somatic_Lib_iHMLYMP'],i,i)
		step3 = 'perl %s --step step3 --input %s.vcf.focus.xls --outdir %s --yourgroup %s.sample.conf --somatic true --depth 200 --db iHMALLdb --selfdb %s > %s.vcf.final.xls\n'%(self.softwares['Analyses.2024.pl'],i,zhushidir,i,self.softwares['RNA.selfdb'],i)

		if i == 'phlike':
			return ' && \\\n'.join([step1,step2_2,step3]),order
		elif i == "iHMLYMP":
			return ' && \\\n'.join([step1,step2_3,step3]),order
		else:
			return ' && \\\n'.join([step1,step2_1,step3]),order

	def zhushi_solid(self,zhushidir,analydir,t_type,i,TR):
		order=[]
		for each in t_type[i]:
			order.append('order zhushi_%s after mutation_calling_%s' %(i,each))
		sample_tmp = []
		dirtmp = analydir+'/Mutation'
		vcf_file=open("%s/%s.vcf.conf"%(zhushidir,i),'w')
		sample_file=open("%s/%s.sample.conf"%(zhushidir,i),'w')
		for each in t_type[i]:
			if not re.search("NEG",each):
				vcf_file.write("%s/%s/%s_%s.filter.somatic.vcf\t%s\n"%(dirtmp,each,each,i,each))
				sample_file.write("case\t%s\t%s/%s/%s_%s.filter.somatic.vcf\n"%(each,dirtmp,each,each,i))
			else:
				vcf_file.write("%s/%s/%s.vars.vcf\t%s\n"%(dirtmp,each,each,each))
				sample_file.write("case\t%s\t%s/%s/%s.vars.vcf\n"%(each,dirtmp,each,each))
		vcf_file.close()
		sample_file.close()
		################################
		cmd1 = ' && \\\n'.join(['export PATH=/public/software/anaconda3_new/envs/perl_5.26.2/bin:$PATH\ncd %s\nperl %s %s.vcf.conf %s.vcf %s'%(zhushidir,self.softwares['vcfmerge.pl'],i,i,self.softwares['refData']),
			'sed -i \'s#1/0#0/1#g\' %s.vcf'%(i),
			'perl %s --step step1 --AF 0.01 --input %s.vcf --db iHMALLdb'%(self.softwares['Analyses.2024.pl'],i),
			'perl %s --step step2 --AF 0.01 --Func splicing,exonic --input %s.vcf --printall true --amplicon %s --yourlib %s/XHHM_%s.Lib.txt > %s.vcf.focus.xls 2> %s.vcf.ignore.xls'%(self.softwares['Analyses.2024.pl'],i,TR,self.softwares['XP_784_Somatic_Lib'],i,i,i),
			'perl %s --step step3 --AF 0.01 --input %s.vcf.focus.xls --outdir %s --yourgroup %s.sample.conf --somatic true --depth 200 --db iHMALLdb --selfdb %s > %s.vcf.final.xls\n'%(self.softwares['Analyses.2024.pl'],i,zhushidir,i,self.softwares['solid_selfdb.txt'],i)])
		return cmd1,order


	def filter_cmd(self,zhushidir,i,each,disease,analydir):
		order=[]
		#sample_type = 'organization'
		order.append('order filter_cmd after zhushi_%s' %(i))
		if disease[each]=="phlike":
			disease_each="phlikeRNA"
		else:
			disease_each=disease[each]
		cmd1 = ' && \\\n'.join(['perl %s Sample_%s.step3.xls %s ' % (self.softwares['iHMAP_XHHM'],each,self.softwares['AA_position_iHMAP']),
			'python %s Sample_%s.filter.tmp.result'%(self.softwares['anno.transvar'],each),
			'perl %s Sample_%s.filter.tmp.result.transvar'%(self.softwares['filter.result.sort.iHM303'],each),
			'python3 %s -i Sample_%s.filter.tmp.result.transvar -s Sample_%s.%s.filter.tmp.result -p %s -n %s'%(self.softwares['filter.genelist.false.postive'],each,each,disease[each],disease_each,each),
			'python %s -i Sample_%s.%s.filter.tmp.result -o Sample_%s.%s.filter.tmp.result.filter'%(self.softwares['result.filter'],each,disease[each],each,disease[each]), 
			'perl %s Sample_%s.%s.filter.tmp.result.filter'%(self.softwares['filter.result.sort.iHM303'],each,disease[each]),
			#'file=Sample_%s.%s.filter.result.sort.xls\nif [ -f "$file" ] && [ -s "$file" ] && [ $(tail -n +2 "$file" | grep -c \'[^[:space:]]\') -eq 0 ]; then\necho "阴性" >> "$file"\nfi'% (each,disease[each]),
            'if [[ $(wc -l Sample_%s.%s.filter.tmp.result.filter |awk \'{print $1}\') -eq 1 ]];then ls -s Sample_%s.%s.filter.tmp.result.filter |awk \'{print $10 "阴性"}\' >> Sample_%s.%s.filter.result.sort.xls;fi' % (each,disease[each],each,disease[each],each,disease[each]),
			'perl %s Sample_%s.step3.xlsx Sample_%s.filter.result.sort.xls Sample_%s.%s.filter.result.sort.xls %s Sample_%s.out.xlsx\n\n'%(self.softwares['excel.combine.panel.pl'],each,each,each,disease[each],disease[each],each)])
		cmd2 = ' && \\\n'.join(['perl %s Sample_%s.step3.xls %s ' % (self.softwares['iHMLYMP_XHHM'],each,self.softwares['AA_position_phlike']),
			'python %s Sample_%s.filter.tmp.result'%(self.softwares['anno.transvar'],each),
			'perl %s Sample_%s.filter.tmp.result.transvar'%(self.softwares['filter.result.sort.iHM303'],each),
			'python3 %s -i Sample_%s.filter.tmp.result.transvar -s Sample_%s.%s.filter.tmp.result -p %s -n %s'%(self.softwares['filter.genelist.false.postive'],each,each,disease[each],disease_each,each),
			'python %s -i Sample_%s.%s.filter.tmp.result -o Sample_%s.%s.filter.tmp.result.filter'%(self.softwares['result.filter'],each,disease[each],each,disease[each]),
			'perl %s Sample_%s.%s.filter.tmp.result.filter'%(self.softwares['filter.result.sort.iHM303'],each,disease[each]),
			#'file=Sample_%s.%s.filter.result.sort.xls\nif [ -f "$file" ] && [ -s "$file" ] && [ $(tail -n +2 "$file" | grep -c \'[^[:space:]]\') -eq 0 ]; then\necho "阴性" >> "$file"\nfi'% (each,disease[each]),
            'if [[ $(wc -l Sample_%s.%s.filter.tmp.result.filter |awk \'{print $1}\') -eq 1 ]];then ls -s Sample_%s.%s.filter.tmp.result.filter |awk \'{print $10 "阴性"}\' >> Sample_%s.%s.filter.result.sort.xls;fi' % (each,disease[each],each,disease[each],each,disease[each]),
			'perl %s Sample_%s.step3.xlsx Sample_%s.filter.result.sort.xls Sample_%s.%s.filter.result.sort.xls %s Sample_%s.out.xlsx\n\n'%(self.softwares['excel.combine.panel.pl'],each,each,each,disease[each],disease[each],each)])
		cmd3 = ' && \\\n'.join(['perl %s Sample_%s.step3.xls %s ' % (self.softwares['iHMLYMP_XHHM'],each,self.softwares['AA_position_iHMLYMP']),
			'python %s Sample_%s.filter.tmp.result'%(self.softwares['anno.transvar'],each),
			'perl %s Sample_%s.filter.tmp.result.transvar'%(self.softwares['filter.result.sort.iHM303'],each),
			'python3 %s -i Sample_%s.filter.tmp.result.transvar -s Sample_%s.%s.filter.tmp.result -p %s -n %s'%(self.softwares['filter.genelist.false.postive'],each,each,disease[each],disease_each,each),
			'python %s -i Sample_%s.%s.filter.tmp.result -o Sample_%s.%s.filter.tmp.result.filter'%(self.softwares['result.filter'],each,disease[each],each,disease[each]),
			'perl %s Sample_%s.%s.filter.tmp.result.filter'%(self.softwares['filter.result.sort.iHM303'],each,disease[each]),
			#'file=Sample_%s.%s.filter.result.sort.xls\nif [ -f "$file" ] && [ -s "$file" ] && [ $(tail -n +2 "$file" | grep -c \'[^[:space:]]\') -eq 0 ]; then\necho "阴性" >> "$file"\nfi'% (each,disease[each]),
            'if [[ $(wc -l Sample_%s.%s.filter.tmp.result.filter |awk \'{print $1}\') -eq 1 ]];then ls -s Sample_%s.%s.filter.tmp.result.filter |awk \'{print $10 "阴性"}\' >> Sample_%s.%s.filter.result.sort.xls;fi' % (each,disease[each],each,disease[each],each,disease[each]),
			'perl %s Sample_%s.step3.xlsx Sample_%s.filter.result.sort.xls Sample_%s.%s.filter.result.sort.xls %s Sample_%s.out.xlsx\n\n'%(self.softwares['excel.combine.panel.pl'],each,each,each,disease[each],disease[each],each)])
		##########
		if i == 'phlike':
			return cmd2,order
		elif i == "iHMLYMP":
			return cmd3,order
		else:
			return cmd1,order

	def filter_cmd_solid(self,zhushidir,i,each,disease,analydir,type_tmp):
		order=[]
		if type_tmp == 'PM':
			sample_type = 'cfDNA'
		elif type_tmp == 'UMI':
			sample_type = 'UMI'
		else:
			sample_type = 'organization'
		order.append('order filter_cmd after zhushi_%s' %(i))
		cmd1 = ' && \\\n'.join([
			'python %s --infile Sample_%s.step3.xls --outfile Sample_%s.filter.tmp.result --tumor_type %s --sample_type %s'%(self.softwares['tumor_XHHM.XP.filter'],each,each,i,sample_type),
			'python %s Sample_%s.filter.tmp.result'%(self.softwares['anno_transvar_solid'],each),
			'perl %s Sample_%s.filter.tmp.result.transvar'%(self.softwares['filter.result.sort_solid'],each),
			'python %s -i Sample_%s.filter.tmp.result.transvar -s Sample_%s.%s.filter.tmp.result -p %s -n %s'%(self.softwares['filter.genelist_solid'],each,each,disease[each],disease[each],each),
			'python %s -i Sample_%s.%s.filter.tmp.result -o Sample_%s.%s.filter.tmp.result.filter' % (self.softwares['result.filter'],each,disease[each],each,disease[each]),
			'perl %s Sample_%s.%s.filter.tmp.result.filter' % (self.softwares['filter.result.sort_solid'],each,disease[each]),
			'if [[ $(wc -l Sample_%s.%s.filter.tmp.result.filter |awk \'{print $1}\') -eq 1 ]];then ls -s Sample_%s.%s.filter.tmp.result.filter |awk \'{print $10 "阴性"}\' >> Sample_%s.%s.filter.result.sort.xls;fi' % (each,disease[each],each,disease[each],each,disease[each]),
			'perl %s Sample_%s.step3.xlsx Sample_%s.filter.result.sort.xls Sample_%s.%s.filter.result.sort.xls %s Sample_%s.out.xlsx\n\n'%(self.softwares['excel.combine.panel.pl'],each,each,each,disease[each],disease[each],each)])

		return cmd1,order

	def filter_cmd_NEG(self,zhushidir,i,each,disease,analydir):
		order=[]
		order.append('order filter_cmd after zhushi_%s' %(i))
		cmd1 = '\n'.join(['\nif [[ -f Sample_%s.step3.xls && -s Sample_%s.step3.xls ]]; then' % (each,each),
				'perl %s Sample_%s.step3.xls %s ' % (self.softwares['iHMAP_XHHM'],each,self.softwares['AA_position_iHMAP']),
				'python %s Sample_%s.filter.tmp.result'%(self.softwares['anno.transvar'],each),
				'perl %s Sample_%s.filter.tmp.result.transvar'%(self.softwares['filter.result.sort.iHM303'],each),
				'python3 %s -i Sample_%s.filter.tmp.result.transvar -s Sample_%s.%s.filter.tmp.result -p %s -n %s'%(self.softwares['filter.genelist.false.postive'],each,each,'iHMAP','iHMAP',each),
				'python %s -i Sample_%s.%s.filter.tmp.result -o Sample_%s.%s.filter.tmp.result.filter'%(self.softwares['result.filter'],each,'iHMAP',each,'iHMAP'),
				'perl %s Sample_%s.%s.filter.tmp.result.filter'%(self.softwares['filter.result.sort.iHM303'],each,'iHMAP'),
				'if [[ $(wc -l Sample_%s.%s.filter.tmp.result.filter |awk \'{print $1}\') -eq 1 ]];then ls -s Sample_%s.%s.filter.tmp.result.filter |awk \'{print $10 "阴性"}\' >> Sample_%s.%s.filter.result.sort.xls;fi' % (each,'iHMAP',each,'iHMAP',each,'iHMAP'),
				'perl %s Sample_%s.step3.xlsx Sample_%s.filter.result.sort.xls Sample_%s.%s.filter.result.sort.xls %s Sample_%s.out.xlsx'%(self.softwares['excel.combine.panel.pl'],each,each,each,'iHMAP','iHMAP',each),
				'fi\n\n'])
		return cmd1,order



class Fusiongene:

	def __init__(self,sampleID,qcDir,alignDir,fusionDir,softwares,databases,tumor_type,analydir,NP,disease):
		self.sampleID = sampleID
		self.qcDir = qcDir
		self.alignDir = alignDir
		self.fusionDir = fusionDir
		self.Root = os.path.dirname(fusionDir)
		self.bam = os.path.join(alignDir, self.sampleID+'.bam')
		self.softwares = softwares
		self.databases = databases
		self.tumor_type = tumor_type
		self.analydir = analydir
		self.NP = NP
		self.disease = disease


	def construct_fq(self):
		fq1 = os.path.join(self.qcDir,'_'.join([self.sampleID,'1.clean.fq.gz']))
 		fq2 = os.path.join(self.qcDir,'_'.join([self.sampleID,'2.clean.fq.gz']))
		return fq1,fq2

	def fusion_step1(self):
		order = 'order fusion_step1_%s after mapping_report' % self.sampleID
		clean_fq1,clean_fq2 = self.construct_fq()
                step0 = '%s\n' % self.softwares['perl_library']
		step1 = '\\\n\t'.join(['%s' % self.softwares['starfusion'],
			'--genome_lib_dir %s' % self.databases['starfusion'],
			'--left_fq %s' % clean_fq1,
			'--right_fq %s' % clean_fq2, 
			'--output_dir %s' % self.fusionDir])
		step2 = 'less %s/star-fusion.fusion_predictions.abridged.tsv|grep -v LeftBreakpoint|awk -F \'\\t\' \'{print $6"\\n"$8}\'|awk -F \':\' \'{print $1"\\t"$2"\\t"$2"\\t0\\t0\\t%s"}\'|sort -u |awk -F "\\t" \'($1~/^chr[XY0-9]+$/){sub("^chr","",$0);print $0}\' |sort -k1,1n -k2,2n|awk \'{print "chr"$0}\' |sed \'/^#/d\' > %s/star-fusion.fusion_predictions.abridged.tsv.cvt'%(self.fusionDir,self.sampleID,self.fusionDir)
		step3 = '%s %s/star-fusion.fusion_predictions.abridged.tsv.cvt -out %s/star-fusion.fusion_predictions.abridged.tsv.anno %s -buildver hg19 -otherinfo -remove -nastring . -protocol refGene,gff3,dgvMerged,exon_intron,cytoBand,interproscan --gff3dbfile hg19_rmsk.gff -operation g,r,r,r,r,f '%(self.softwares['table_annovar'],self.fusionDir,self.fusionDir,self.softwares['iHMALLdb'])
		return step0+' && \\\n'.join([step1,step2,step3]),order


	def fusion_step2(self,version,sample_db):
		#### The casevar database used for target and mRNA is different
		if re.search("^RH|^RS|^RA",self.sampleID):
			self.casevar = "targetRNA.star-fusion.casevar.list"
			self.false_positive = "targetRNA.Starfusion_False_positive_Site.xls"
		elif re.search("^RT",self.sampleID):
			self.casevar = "solid.targetRNA.star-fusion.casevar.list"
			self.false_positive = "solid.targetRNA.Starfusion_False_positive_Site.xls"
		elif re.search("^R",self.sampleID) and not re.search("^RH|^RS|^RA|^RT",self.sampleID):
			self.casevar = "mRNA.star-fusion.casevar.list"
			self.false_positive = "mRNA.Starfusion_False_positive_Site.xls"
		else:
			print "please ",self.sampleID

		#if self.sampleID.startswith('RH') and self.disease[self.sampleID] == 'phlike':  
		#	fusion_hotspot = self.softwares['fusion_hotspot_phlike']
		if self.sampleID.startswith('RT'):
			cancer_type = sample_db[self.sampleID]
			fusion_hotspot = self.softwares['fusion_hotspot_solid'] + "Fusion_" + cancer_type +".txt"
		else:
			fusion_hotspot = self.softwares['fusion_hotspot']

		if self.sampleID.startswith('RT'):
			fusion_whitelist = self.softwares['script'] + "/RNAfusion/fusiongene.whitelist_solid"
		else:
			fusion_whitelist = self.softwares['script'] + "/RNAfusion/fusiongene.whitelist_blood"
	
		order = 'order fusion_step2_%s after fusion_step1_%s' % (self.sampleID,self.sampleID)
		step1 = '%s\ncd %s' % (self.softwares['perl_library'],self.fusionDir)
		step2 = 'python %s -i star-fusion.fusion_predictions.abridged.tsv -a star-fusion.fusion_predictions.abridged.tsv.anno.hg19_multianno.txt -hot %s -o %s.star-fusion.fusion_predictions.abridged.tsv.out -b %s -case %s' % (self.softwares['fusion_anno_blood'],fusion_hotspot,self.sampleID,os.path.join(self.alignDir,self.sampleID+'.bam'),os.path.join(self.softwares['Fusion_Casevar_path'],self.casevar))
		####
		step3 = '### filter\npython %s -i %s.star-fusion.fusion_predictions.abridged.tsv.out -fp %s -o %s.star-fusion.fusion_predictions.abridged.tsv.filter.result -hot %s -wl %s' % (self.softwares['fusion_filter'],self.sampleID,os.path.join(self.softwares['False_positive_path'],self.false_positive),self.sampleID,fusion_hotspot,fusion_whitelist)
		step4 = 'perl %s %s.star-fusion.fusion_predictions.abridged.tsv.filter.result' % (self.softwares['filter.result.sort'],self.sampleID)
		step5 = '### \nperl %s %s/Mapping/%s/%s_qc_map.stat %s.star-fusion.fusion_predictions.abridged.tsv.out %s.star-fusion.fusion_predictions.abridged.tsv.filter.result.sort.xls %s.star-fusion.xlsx' % (self.softwares['excel.combine.fusion'],self.analydir,self.sampleID,self.sampleID,self.sampleID,self.sampleID,self.sampleID)
		####
		step6 = '###\nsamtools sort -@ 2 Aligned.sortedByCoord.out.bam %s.Starfusion.sort && \\\nsamtools index %s.Starfusion.sort.bam ' % (self.sampleID,self.sampleID)
		step7 = 'if [ $(tail -c 28 %s.Starfusion.sort.bam | /public/software/tools/xxd -p) = "1f8b08040000000000ff0600424302001b0003000000000000000000" ]; then\n\trm Aligned.sortedByCoord.out.bam\nfi\n'  % self.sampleID

		return ' && \\\n'.join([step1,step2,step3,step4,step5,step6,step7]),order

	def soapfuse_prepare(self,version,paired='paired',length='150',mode='normal'):
		order = 'order soapfuse_prepare after mapping_report'
		if mode=='normal':
			step1 = 'awk -F "\\t" \'{print $3"\\t%s\\t%s"}\' %s/qc_list > %s/sample.txt' % (paired,length,self.analydir,self.fusionDir)
			step2 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['SOAPfuse']),
				'--version %s' % version,
				'--database %s' % self.databases['soapfuse'],
				'--pipline %s' % os.path.join(self.softwares['bin'],'soapfuse_pipeline'),
				'--rootdir %s' % self.analydir,
				'--NP %s '% self.NP])
		else:
			step1 = 'awk -F "\\t" \'{print $7"\\t"$2"\\t"$3"\\t%s\\t%s"}\' %s/qc_list > %s/sample.txt' % (paired,length,self.analydir,self.fusionDir)
			step2 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['SOAPfuse']),
				'--version %s' % version,
				'--database %s' % self.databases['soapfuse'],
				'--pipline %s' % os.path.join(self.softwares['bin'],'soapfuse_pipeline'),
				'--mode somatic',
				'--rootdir %s' % self.analydir,
				'--NP %s '% self.NP])
		return ' && \\\n'.join([step1,step2]),order


	def arriba_step1(self,sampleID):
		order = 'order arriba_step1_%s after mapping_report'%sampleID
		step1 = '%s\ncd %s\nln -sf %s/*.gz %s/' % (self.softwares['perl_library'],self.fusionDir,self.qcDir,self.fusionDir)
		step2 = '\\\n\t'.join(['\n## STAR + arriba\n%s ' % self.softwares['STAR'],
			'--runThreadN 8 ',
			'--genomeDir %s --genomeLoad NoSharedMemory ' % self.softwares['STAR-index'],
			'--readFilesIn %s_1.clean.fq.gz %s_2.clean.fq.gz --readFilesCommand zcat '% (sampleID,sampleID),
			'--outStd BAM_Unsorted --outSAMtype BAM Unsorted --outSAMunmapped Within --outBAMcompression 0 ',
			'--outFilterMultimapNmax 50 --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 --alignSJstitchMismatchNmax 5 -1 5 5 ',
			'--chimSegmentMin 10 --chimOutType WithinBAM HardClip --chimJunctionOverhangMin 10 --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 ',
			'--chimSegmentReadGapMax 3 --chimMultimapNmax 50 | ',
			'%s -x /dev/stdin -f blacklist -o fusions.tsv -O fusions.discarded.tsv -a %s -g %s' % (self.softwares['Arriba'],self.databases['fasta'],self.databases['gtf'])])
		step3 = '\\\n'.join([
            '\n### integrate\npython %s -f fusions.discarded.tsv -o retrieve.fusions.tsv &&' % self.softwares['retrieve.fusions'],
            'cat fusions.tsv retrieve.fusions.tsv > total.fusions.tsv &&\\\n'
			'## anno\ncat total.fusions.tsv |awk -F\'\\t\' \'$5 !~ /breakpoint/{print $5"\\n"$6}\' | sed \'s/:/\\t/\' | sort -u | sort -k1V,1 -k2,2n > mutations.all.sort.lst &&',
			'cat mutations.all.sort.lst  |awk \'BEGIN {OFS="\\t"} {print $1,$2,$2+1,"0","0"}\'   > mutations.all.sort.tsv &&',
			'%s mutations.all.sort.tsv -out mutations.all.sort.tsv.anno %s -buildver hg19 -otherinfo -remove -nastring . -protocol refGene,gff3,dgvMerged,exon_intron,cytoBand,interproscan --gff3dbfile hg19_rmsk.gff -operation g,r,r,r,r,f\n ' % (self.softwares['table_annovar'],self.softwares['iHMALLdb'])])

		return ' && \\\n'.join([step1,step2,step3]),order

	def arriba_step2(self,sample_db,sampleID):
		if re.search("^RH|^RS|^RA",sampleID):
			self.casevar = "targetRNA.Arriba.casevar.list"
			self.false_positive = "targetRNA.Arriba_False_positive_Site.xls"
		elif re.search("^RT",sampleID):
			self.casevar = "solid.targetRNA.Arriba.casevar.list"
			self.false_positive = "solid.targetRNA.Arriba_False_positive_Site.xls"
		elif re.search("^R",sampleID) and not re.search("^RH|^RS|^RA|^RT",sampleID):
			self.casevar = "mRNA.Arriba.casevar.list"
			self.false_positive = "mRNA.Arriba_False_positive_Site.xls"
		else:
			print "please ",sampleID
		
		#if self.sampleID.startswith('RH') and self.disease[self.sampleID] == 'phlike':
		#	fusion_hotspot = self.softwares['fusion_hotspot_phlike']
		if sampleID.startswith('RT'):
			fusion_hotspot = self.softwares['fusion_hotspot_solid'] + "Fusion_" + sample_db[sampleID] +".txt"
		else:
			fusion_hotspot = self.softwares['fusion_hotspot']

		if sampleID.startswith('RT'):
			fusion_whitelist = self.softwares['script'] + "/RNAfusion/fusiongene.whitelist_solid"
			rearrangement_genelist = self.softwares['solid_rearrangement_genelist']
		else:            
			fusion_whitelist = self.softwares['script'] + "/RNAfusion/fusiongene.whitelist_blood"
			rearrangement_genelist = self.softwares['rearrangement_genelist']

		
		order = 'order arriba_step2_%s after arriba_step1_%s'%(sampleID,sampleID)
		step1 = '%s\ncd %s' % (self.softwares['perl_library'],self.fusionDir)
		step2 = '\\\n'.join([
			'### hot anno\npython %s -f total.fusions.tsv -a mutations.all.sort.tsv.anno.hg19_multianno.txt -case %s -hot %s &&' % (self.softwares['arriba_hotanno'],os.path.join(self.softwares['Fusion_Casevar_path'],self.casevar),fusion_hotspot),
			'python %s mutations.all.sort.tsv.gene_hotspot.anno.raw.txt mutations.all.sort.tsv.gene_hotspot.anno.raw.sort.txt &&' % self.softwares['arriba_hotanno_sort'],
			'## filter\npython %s -i mutations.all.sort.tsv.gene_hotspot.anno.raw.sort.txt -o mutations.all.sort.tsv.gene_hotspot.anno.filter.sort.txt -hot %s -fp %s -rg %s -wl %s &&' % (self.softwares['arriba_filter'],fusion_hotspot,os.path.join(self.softwares['False_positive_path'],self.false_positive),rearrangement_genelist,fusion_whitelist),
			'## class\npython %s -i mutations.all.sort.tsv.gene_hotspot.anno.filter.sort.txt -o breakpoint.fusiongene.result  &&' % self.softwares['arriba_fusion'],
			'perl %s breakpoint.fusiongene.result &&' % self.softwares['result_sort'],
			'python %s -i mutations.all.sort.tsv.gene_hotspot.anno.filter.sort.txt -o rearrangement.fusiongene.result -g %s &&' % (self.softwares['arriba_rearrangement'],rearrangement_genelist),
			'perl %s rearrangement.fusiongene.result &&' % self.softwares['result_sort'],
			'## combine\nperl %s %s/Mapping/%s/%s_qc_map.stat mutations.all.sort.tsv.gene_hotspot.anno.raw.sort.txt breakpoint.fusiongene.result.sort.xls rearrangement.fusiongene.result.sort.xls %s.arriba.final.xlsx\n' % (self.softwares['arriba_comb'],self.analydir,sampleID,sampleID,sampleID)])
		
		return ' && \\\n'.join([step1,step2]),order



class ASdifferential:

	def __init__(self,groupTs,groupNs,groupTname,groupNname,grp_name,alignDir,asdiffDir,softwares,databases):
		self.groupTs = groupTs
		self.groupNs = groupNs
		self.groupTname = groupTname
		self.groupNname = groupNname
		self.grp_name = grp_name
		self.alignDir = alignDir
		self.asdiffDir = asdiffDir
		self.diffRoot = os.path.dirname(asdiffDir)
		self.diffName = '%s_vs_%s' % (groupTname,groupNname)
		self.softwares = softwares
 		self.databases = databases

	def get_parameter(self):

		if len(self.groupTs) >= 3 and len(self.groupTs)==len(self.groupNs):
			analysis='P'
		else:
			analysis='U'

		insertmeanT=[]
		deviationsT=[]

		insertmeanC=[]
		deviationsC=[]
		for eachS in self.groupTs:
			tmp=open(os.path.join(self.alignDir,eachS,eachS+'.insertsize.stat')).read().split('\t')
			eachI=tmp[1].split('.')[0]
			insertmeanT.append(eachI)

			eachD=tmp[2].split('.')[0]
			deviationsT.append(eachD)

		insertmeanT=','.join(insertmeanT)
		deviationsT=','.join(deviationsT)
		for eachS in self.groupNs:
			tmp=open(os.path.join(self.alignDir,eachS,eachS+'.insertsize.stat')).read().split('\t')
			eachI=tmp[1].split('.')[0]
			insertmeanC.append(eachI)

			eachD=tmp[2].split('.')[0]
			deviationsC.append(eachD)

		insertmeanC=','.join(insertmeanC)
		deviationsC=','.join(deviationsC)
		return insertmeanT,insertmeanC,deviationsT,deviationsC,analysis
	
	def rnaseq_mats(self,analysis,t='paired',len='150',libtype='fr-firststrand'):
		#order = ['order rnaseq_mats_%s after rnaseq_mapping_%s' % (self.diffName,each) for each in self.groupNs+self.groupTs]
		order = 'order rnaseq_mats_%s after mapping_report' % self.grp_name
		step1_old = '''
T=(%s)
for i in ${T[*]}
do
	eachIT=$(cut -f2 %s/$i/$i.insertsize.stat|awk -F '.' '{print $1}')
	#insertmeanT=(${insertmeanT[@]} $eachI)
	insertmeanT=${insertmeanT}","${eachIT}
	eachDT=$(cut -f3 %s/$i/$i.insertsize.stat|awk -F '.' '{print $1}')
	#deviationsT=(${deviationsT[@]} $eachD)
	deviationsT=${deviationsT}","${eachDT}
done
N=(%s)			
for i in ${N[*]}
do
	eachIN=$(cut -f2 %s/$i/$i.insertsize.stat|awk -F '.' '{print $1}')
	insertmeanN=${insertmeanN}","${eachIN}
	eachDN=$(cut -f3 %s/$i/$i.insertsize.stat|awk -F '.' '{print $1}')
	deviationsN=${deviationsN}","${eachDN}
done

insertmeanT=`echo ${insertmeanT:1}`
deviationsT=`echo ${deviationsT:1}`
insertmeanN=`echo ${insertmeanN:1}`
deviationsN=`echo ${deviationsN:1}`

			''' % (' '.join(self.groupTs),self.alignDir,self.alignDir,' '.join(self.groupNs),self.alignDir,self.alignDir)
		### 20250408 : reads len 
		egExample = os.path.join(self.alignDir,self.groupTs[0],self.groupTs[0]+'.bam')
		step1 = "readsseq=`samtools view %s | head -1 |cut -f 10`\nlength=${#readsseq}\n" % egExample

		step2_old = '\\\n\t'.join(['python %s ' % self.softwares['rmats'],
			'-b1 %s' % ','.join([os.path.join(self.alignDir,each,each+'.bam') for each in self.groupTs]),
			'-b2 %s' % ','.join([os.path.join(self.alignDir,each,each+'.bam') for each in self.groupNs]),
			'-gtf %s' % self.databases['gtf'],
			'-t %s -len %s -a 8 -r1 %s -r2 %s -sd1 %s -sd2 %s -c 0.01 -analysis %s -libType %s' % (t,len,'$insertmeanT','$insertmeanN','$deviationsT','$deviationsN',analysis,libtype),
			'-o %s' % self.asdiffDir])
        #####
		step2 = '\\\n\t'.join(['python %s ' % self.softwares['rmats'],
			'-b1 %s' % ','.join([os.path.join(self.alignDir,each,each+'.bam') for each in self.groupTs]),
			'-b2 %s' % ','.join([os.path.join(self.alignDir,each,each+'.bam') for each in self.groupNs]),
			'-gtf %s' % self.databases['gtf'],
			'-t %s -len $length -a 8 -c 0.01 -analysis %s -libType %s' % (t,analysis,libtype),
			'-o %s' % self.asdiffDir])
		step3 = 'mv %s/MATS_output/SE.MATS.ReadsOnTargetAndJunctionCounts.txt %s/%s.SE.xls' % (self.asdiffDir,self.asdiffDir,self.grp_name)
		step4 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['rMATSresult']),
			'--infile %s/%s.SE.xls' % (self.asdiffDir,self.grp_name),
			'--cutoff 0.01',
			'--outfile %s/%s.SE.significant.xls' % (self.asdiffDir,self.grp_name)])
		step4_plot_0='if [ -d %s ];then rm -r %s && mkdir -p %s;else mkdir -p %s;fi' % (os.path.join(self.asdiffDir,'SE_plot'),os.path.join(self.asdiffDir,'SE_plot'),os.path.join(self.asdiffDir,'SE_plot'),os.path.join(self.asdiffDir,'SE_plot'))
		step4_plot_1='number=$(wc -l %s/%s.SE.significant.xls|cut -d \' \' -f 1)' % (self.asdiffDir,self.grp_name)
		step4_plot_2='if [[ $number -gt 6 ]]; then\n    head -n 6 %s/%s.SE.significant.xls > %s/SE_plot/%s.SE.plot.txt\nelif [[ $number -eq 1 ]];then\n    echo \'%s.SE.significant.xls is empty \'\nelse ln -s %s/%s.SE.significant.xls %s/SE_plot/%s.SE.plot.txt\nfi' % (self.asdiffDir,self.grp_name,self.asdiffDir,self.grp_name,self.grp_name,self.asdiffDir,self.grp_name,self.asdiffDir,self.grp_name)
		step4_plot_3='\\\n\t'.join([self.softwares['rmats2sashimiplot'],
			'--b1 %s' % ','.join([os.path.join(self.alignDir,each,each+'.bam') for each in self.groupTs]),
			'--b2 %s' % ','.join([os.path.join(self.alignDir,each,each+'.bam') for each in self.groupNs]),
			'-t SE',
			'-e %s/SE_plot/%s.SE.plot.txt' % (self.asdiffDir,self.grp_name),
			'--l1 %s' % self.groupTname,
			'--l2 %s' % self.groupNname,
			'-o %s' % os.path.join(self.asdiffDir,'SE_plot')])
		step4_plot_4='rm -rf %s/Sashimi_index*' % os.path.join(self.asdiffDir,'SE_plot')
		step5 = 'mv %s/MATS_output/RI.MATS.ReadsOnTargetAndJunctionCounts.txt %s/%s.RI.xls' % (self.asdiffDir,self.asdiffDir,self.grp_name)
		step6 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['rMATSresult']),
			'--infile %s/%s.RI.xls' % (self.asdiffDir,self.grp_name),
			'--cutoff 0.01',
			'--outfile %s/%s.RI.significant.xls' % (self.asdiffDir,self.grp_name)])
		step6_plot_0='if [ -d %s ];then rm -r %s && mkdir -p %s;else mkdir -p %s;fi' % (os.path.join(self.asdiffDir,'RI_plot'),os.path.join(self.asdiffDir,'RI_plot'),os.path.join(self.asdiffDir,'RI_plot'),os.path.join(self.asdiffDir,'RI_plot'))
		step6_plot_1='number=$(wc -l %s/%s.RI.significant.xls|cut -d \' \' -f 1)' % (self.asdiffDir,self.grp_name)
		step6_plot_2='if [[ $number -gt 6 ]]; then\n    head -n 6 %s/%s.RI.significant.xls > %s/RI_plot/%s.RI.plot.txt\nelif [[ $number -eq 1 ]];then\n    echo \'%s.RI.significant.xls is empty \'\nelse ln -s %s/%s.RI.significant.xls %s/RI_plot/%s.RI.plot.txt\nfi' % (self.asdiffDir,self.grp_name,self.asdiffDir,self.grp_name,self.grp_name,self.asdiffDir,self.grp_name,self.asdiffDir,self.grp_name)
		step6_plot_3='\\\n\t'.join([self.softwares['rmats2sashimiplot'],
			'--b1 %s' % ','.join([os.path.join(self.alignDir,each,each+'.bam') for each in self.groupTs]),
			'--b2 %s' % ','.join([os.path.join(self.alignDir,each,each+'.bam') for each in self.groupNs]),
			'-t RI',
			'-e %s/RI_plot/%s.RI.plot.txt' % (self.asdiffDir,self.grp_name),
			'--l1 %s' % self.groupTname,
			'--l2 %s' % self.groupNname,
			'-o %s' % os.path.join(self.asdiffDir,'RI_plot')])
		step6_plot_4='rm -rf %s/Sashimi_index*' % os.path.join(self.asdiffDir,'RI_plot')
		step7 = 'mv %s/MATS_output/MXE.MATS.ReadsOnTargetAndJunctionCounts.txt %s/%s.MXE.xls' % (self.asdiffDir,self.asdiffDir,self.grp_name)
		step8 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['rMATSresult']),
			'--infile %s/%s.MXE.xls' % (self.asdiffDir,self.grp_name),
			'--cutoff 0.01',
			'--outfile %s/%s.MXE.significant.xls' % (self.asdiffDir,self.grp_name)])
		step8_plot_0='if [ -d %s ];then rm -r %s && mkdir -p %s;else mkdir -p %s;fi' % (os.path.join(self.asdiffDir,'MXE_plot'),os.path.join(self.asdiffDir,'MXE_plot'),os.path.join(self.asdiffDir,'MXE_plot'),os.path.join(self.asdiffDir,'MXE_plot'))
		step8_plot_1='number=$(wc -l %s/%s.MXE.significant.xls|cut -d \' \' -f 1)' % (self.asdiffDir,self.grp_name)
		step8_plot_2='if [[ $number -gt 6 ]]; then\n    head -n 6 %s/%s.MXE.significant.xls > %s/MXE_plot/%s.MXE.plot.txt\nelif [[ $number -eq 1 ]];then\n    echo \'%s.MXE.significant.xls is empty \'\nelse ln -s %s/%s.MXE.significant.xls %s/MXE_plot/%s.MXE.plot.txt\nfi' % (self.asdiffDir,self.grp_name,self.asdiffDir,self.grp_name,self.grp_name,self.asdiffDir,self.grp_name,self.asdiffDir,self.grp_name)
		step8_plot_3='\\\n\t'.join([self.softwares['rmats2sashimiplot'],
			'--b1 %s' % ','.join([os.path.join(self.alignDir,each,each+'.bam') for each in self.groupTs]),
			'--b2 %s' % ','.join([os.path.join(self.alignDir,each,each+'.bam') for each in self.groupNs]),
			'-t MXE',
			'-e %s/MXE_plot/%s.MXE.plot.txt' % (self.asdiffDir,self.grp_name),
			'--l1 %s' % self.groupTname,
			'--l2 %s' % self.groupNname,
			'-o %s' % os.path.join(self.asdiffDir,'MXE_plot')])
		step8_plot_4='rm -rf %s/Sashimi_index*' % os.path.join(self.asdiffDir,'MXE_plot')

		step9 = 'mv %s/MATS_output/A5SS.MATS.ReadsOnTargetAndJunctionCounts.txt %s/%s.A5SS.xls' % (self.asdiffDir,self.asdiffDir,self.grp_name)
		step10 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['rMATSresult']),
			'--infile %s/%s.A5SS.xls' % (self.asdiffDir,self.grp_name),
			'--cutoff 0.01',
			'--outfile %s/%s.A5SS.significant.xls' % (self.asdiffDir,self.grp_name)])
		step10_plot_0='if [ -d %s ];then rm -r %s && mkdir -p %s;else mkdir -p %s;fi' % (os.path.join(self.asdiffDir,'A5SS_plot'),os.path.join(self.asdiffDir,'A5SS_plot'),os.path.join(self.asdiffDir,'A5SS_plot'),os.path.join(self.asdiffDir,'A5SS_plot'))
		step10_plot_1='number=$(wc -l %s/%s.A5SS.significant.xls|cut -d \' \' -f 1)' % (self.asdiffDir,self.grp_name)
		step10_plot_2='if [[ $number -gt 6 ]]; then\n    head -n 6 %s/%s.A5SS.significant.xls > %s/A5SS_plot/%s.A5SS.plot.txt\nelif [[ $number -eq 1 ]];then\n    echo \'%s.A5SS.significant.xls is empty \'\nelse ln -s %s/%s.A5SS.significant.xls %s/A5SS_plot/%s.A5SS.plot.txt\nfi' % (self.asdiffDir,self.grp_name,self.asdiffDir,self.grp_name,self.grp_name,self.asdiffDir,self.grp_name,self.asdiffDir,self.grp_name)
		step10_plot_3='\\\n\t'.join([self.softwares['rmats2sashimiplot'],
			'--b1 %s' % ','.join([os.path.join(self.alignDir,each,each+'.bam') for each in self.groupTs]),
			'--b2 %s' % ','.join([os.path.join(self.alignDir,each,each+'.bam') for each in self.groupNs]),
			'-t A5SS',
			'-e %s/A5SS_plot/%s.A5SS.plot.txt' % (self.asdiffDir,self.grp_name),
			'--l1 %s' % self.groupTname,
			'--l2 %s' % self.groupNname,
			'-o %s' % os.path.join(self.asdiffDir,'A5SS_plot')])
                step10_plot_4='rm -rf %s/Sashimi_index*' % os.path.join(self.asdiffDir,'A5SS_plot')

		step11 = 'mv %s/MATS_output/A3SS.MATS.ReadsOnTargetAndJunctionCounts.txt %s/%s.A3SS.xls' % (self.asdiffDir,self.asdiffDir,self.grp_name)
		step12 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['rMATSresult']),
			'--infile %s/%s.A3SS.xls' % (self.asdiffDir,self.grp_name),
			'--cutoff 0.01',
			'--outfile %s/%s.A3SS.significant.xls' % (self.asdiffDir,self.grp_name)])
		step12_plot_0='if [ -d %s ];then rm -r %s && mkdir -p %s;else mkdir -p %s;fi' % (os.path.join(self.asdiffDir,'A3SS_plot'),os.path.join(self.asdiffDir,'A3SS_plot'),os.path.join(self.asdiffDir,'A3SS_plot'),os.path.join(self.asdiffDir,'A3SS_plot'))
		step12_plot_1='number=$(wc -l %s/%s.A3SS.significant.xls|cut -d \' \' -f 1)' % (self.asdiffDir,self.grp_name)
		step12_plot_2='if [[ $number -gt 6 ]]; then\n    head -n 6 %s/%s.A3SS.significant.xls > %s/A3SS_plot/%s.A3SS.plot.txt\nelif [[ $number -eq 1 ]];then\n    echo \'%s.A3SS.significant.xls is empty \'\nelse ln -s %s/%s.A3SS.significant.xls %s/A3SS_plot/%s.A3SS.plot.txt\nfi' % (self.asdiffDir,self.grp_name,self.asdiffDir,self.grp_name,self.grp_name,self.asdiffDir,self.grp_name,self.asdiffDir,self.grp_name)
		step12_plot_3='\\\n\t'.join([self.softwares['rmats2sashimiplot'],
			'--b1 %s' % ','.join([os.path.join(self.alignDir,each,each+'.bam') for each in self.groupTs]),
			'--b2 %s' % ','.join([os.path.join(self.alignDir,each,each+'.bam') for each in self.groupNs]),
			'-t A3SS',
			'-e %s/A3SS_plot/%s.A3SS.plot.txt' % (self.asdiffDir,self.grp_name),
			'--l1 %s' % self.groupTname,
			'--l2 %s' % self.groupNname,
			'-o %s' % os.path.join(self.asdiffDir,'A3SS_plot')])
                step12_plot_4='rm -rf %s/Sashimi_index*' % os.path.join(self.asdiffDir,'A3SS_plot')


  
		#return step1+'\n'+' && \\\n\n'.join([step2,step3,step4,step4_plot_0,step4_plot_1,step4_plot_2,step4_plot_3,step4_plot_4,step5,step6,step6_plot_0,step6_plot_1,step6_plot_2,step6_plot_3,step6_plot_4,step7,step8,step8_plot_0,step8_plot_1,step8_plot_2,step8_plot_3,step8_plot_4,step9,step10,step10_plot_0,step10_plot_1,step10_plot_2,step10_plot_3,step10_plot_4,step11,step12,step12_plot_0,step12_plot_1,step12_plot_2,step12_plot_3,step12_plot_4]),order
                # cut plot 
                return step1+'\n'+' && \\\n\n'.join([step2,step3,step4,step4_plot_0,step4_plot_1,step4_plot_2,step5,step6,step6_plot_0,step6_plot_1,step6_plot_2,step7,step8,step8_plot_0,step8_plot_1,step8_plot_2,step9,step10,step10_plot_0,step10_plot_1,step10_plot_2,step11,step12,step12_plot_0,step12_plot_1,step12_plot_2+"\n"]),order


	def rnaseq_mats_anno(self):
		order = 'order rnaseq_mats_anno_%s after rnaseq_mats_%s' % (self.grp_name,self.grp_name)
                step0 = '%s\ncd %s' % (self.softwares['perl_library'],self.asdiffDir)
		step1 = 'cat %s/%s.SE.xls|grep -v \'IncLevelDifference\'|awk -F \'\\t\' \'{print $4"\\t"$6"\\t"$7"\\t0\\t0\\t%s\\n"$4"\\t"$8"\\t"$9"\\t0\\t0\\t%s\\n"$4"\\t"$10"\\t"$11"\\t0\\t0\\t%s"}\'|sort -u |awk -F "\\t" \'($1~/^chr[XY0-9]+$/){sub("^chr","",$0);print $0}\' |sort -k1,1n -k2,2n|awk \'{print "chr"$0}\' |sed \'/^#/d\' >%s.SE.cvt' %(self.asdiffDir,self.grp_name,self.grp_name,self.grp_name,self.grp_name,self.grp_name)
		step2 = '%s %s/%s.SE.cvt -out %s/%s.SE.cvt.anno %s -buildver hg19 -otherinfo -remove -nastring . -protocol refGene,gff3,dgvMerged,exon_intron,cytoBand,interproscan --gff3dbfile hg19_rmsk.gff -operation g,r,r,r,r,f' % (self.softwares['table_annovar'],self.asdiffDir,self.grp_name,self.asdiffDir,self.grp_name,self.softwares['iHMALLdb'])
		step3 = 'python %s -i %s/%s.SE.xls -o %s/%s.SE.altersplice.xls -a %s/%s.SE.cvt.anno.hg19_multianno.txt' % (self.softwares['altersplice_anno'],self.asdiffDir,self.grp_name,self.asdiffDir,self.grp_name,self.asdiffDir,self.grp_name)
		##blood
		step4_1 = 'cat %s/%s.SE.altersplice.xls|grep -w -E \'geneSymbol|IKZF1\' > %s/%s.SE.altersplice.filter.xls' % (self.asdiffDir,self.grp_name,self.asdiffDir,self.grp_name)
		step5_1 = 'python %s %s/%s.SE.altersplice.filter.xls %s/%s.SE.altersplice.filter.IKZF1.xls' % (self.softwares['altersplice.filter.blood'],self.asdiffDir,self.grp_name,self.asdiffDir,self.grp_name)
		step6_1 = 'perl  %s %s/%s.SE.altersplice.xls %s/%s.SE.altersplice.filter.xls %s/%s.SE.altersplice.filter.IKZF1.xls IKZF1 %s/%s.SE.altersplice.xlsx\n' % (self.softwares['excel.combine.altersplice'],self.asdiffDir,self.grp_name,self.asdiffDir,self.grp_name,self.asdiffDir,self.grp_name,self.asdiffDir,self.grp_name)
		##solid
		step4_2 = 'cat %s/%s.SE.altersplice.xls|grep -w -E \'geneSymbol|MET\' > %s/%s.SE.altersplice.filter.xls' % (self.asdiffDir,self.grp_name,self.asdiffDir,self.grp_name)
		step5_2 = 'python %s %s/%s.SE.altersplice.filter.xls %s/%s.SE.altersplice.filter.MET.xls' % (self.softwares['altersplice.filter.solid'],self.asdiffDir,self.grp_name,self.asdiffDir,self.grp_name)
		step6_2 = 'perl  %s %s/%s.SE.altersplice.xls %s/%s.SE.altersplice.filter.xls %s/%s.SE.altersplice.filter.MET.xls MET %s/%s.SE.altersplice.xlsx\n' % (self.softwares['excel.combine.altersplice'],self.asdiffDir,self.grp_name,self.asdiffDir,self.grp_name,self.asdiffDir,self.grp_name,self.asdiffDir,self.grp_name)
		
		if self.grp_name.startswith('RT'):
			return ' && \\\n'.join([step0,step1,step2,step3,step4_2,step5_2,step6_2]),order
		else:
			return ' && \\\n'.join([step0,step1,step2,step3,step4_1,step5_1,step6_1]),order


