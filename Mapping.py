## -*- coding:utf-8 -*-
import os,sys
import string

class Mapping:

	def __init__(self,sampleID,qcDir,alignDir,softwares,databases,TR):
		self.sampleID = sampleID
		self.qcDir = qcDir
		self.alignDir = alignDir
		self.TR = TR
		self.bam = os.path.join(alignDir, self.sampleID+'.bam')
		self.sam = os.path.join(alignDir, self.sampleID+'.sam')
		self.softwares = softwares
		self.databases = databases

	def construct_fq(self):
		fq1_gz = os.path.join(self.qcDir,'_'.join([self.sampleID,'1.clean.fq.gz']))
		fq2_gz = os.path.join(self.qcDir,'_'.join([self.sampleID,'2.clean.fq.gz']))
		return fq1_gz,fq2_gz


	def mapping_rnaseq(self,soft='histat2',libtype='fr-firststrand',assemble_soft='cufflinks'):
		order = 'order mapping_rnaseq_%s after qc_%s' % (self.sampleID,self.sampleID)
		clean_fq1,clean_fq2 = self.construct_fq()
		if soft == 'tophat2':
			cmd = '''cd %s
# export PATH=/PUBLIC/software/public/VarCall/samtools/samtools-0.1.18:$PATH
# export PATH=/PUBLIC/software/public/Alignment/bowtie2-2.0.6:$PATH
# export PATH=/PUBLIC/software/RNA/tophat-2.0.9:$PATH
/PUBLIC/software/RNA/tophat-2.0.9/tophat2 \\
  -p 8 -N 2 --read-edit-dist 2 --rg-id %s --rg-sample %s --rg-library %s --rg-platform illumina --rg-platform-unit PU --library-type %s \\
  -G %s \\
  -o %s \\
  %s \\
  %s \\
  %s && \\
#mv %s/accepted_hits.bam %s && \\
/PUBLIC/software/RNA/samtools-1.2/samtools sort %s/accepted_hits.bam %s/%s && \\
/PUBLIC/software/RNA/samtools-1.2/samtools index %s 
#rm unmapped.bam insertions.bed deletions.bed junctions.bed
''' % (self.alignDir,self.sampleID,self.sampleID,self.sampleID,libtype,self.databases['gtf'],self.alignDir,self.databases['fasta'],clean_fq1,clean_fq2,self.alignDir,self.bam,self.alignDir,self.alignDir,self.sampleID,self.bam)
		elif soft == 'hisat2':
			if libtype=='fr-firststrand':
				libchoice='--rna-strandness RF'
			elif libtype=='fr-secondstrand':
				libchoice='--rna-strandness FR'
			else:
				libchoice=''
			if assemble_soft=='cufflinks':
				assembly_type='--dta-cufflinks'
			elif assemble_soft=='stringtie':
				assembly_type='--dta'
			cmd = '''cd %s 
%s \\
  -p 4 --no-unal -t --phred33 --rg-id %s --rg "SM:%s" --rg "PL:illumina" --rg "LB:%s" %s %s \\
  -x %s \\
  -1 %s \\
  -2 %s \\
  -S %s \\
  --un-conc-gz %s && \\\n
samtools view -bS  %s > %s && \\\n
samtools sort %s %s && \\\n
samtools index %s  &&\\\n
rm %s && \\\n
rm %s
''' % (self.alignDir,self.softwares['hisat2'],self.sampleID,self.sampleID,self.sampleID,libchoice,assembly_type,self.databases['fasta'],clean_fq1,clean_fq2,self.sam,os.path.join(self.alignDir,self.sampleID+'.unmapped.fq.gz'),self.sam,self.bam+'.temp.bam',self.bam+'.temp.bam',self.sampleID,self.bam,self.sam,self.bam+'.temp.bam')
		elif soft == 'star':
			rg = 'ID:%s SM:%s LB:%s PL:illumina'%(self.sampleID,self.sampleID,self.sampleID)
			cmd = '''cd %s
rm -rf %s/star
/PUBLIC/software/RNA/STAR/STAR2.5/bin/Linux_x86_64/STAR --readFilesCommand zcat \\
   --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --limitGenomeGenerateRAM 10000000000 \\
  --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFilterMismatchNmax 2 \\
  --runThreadN 8  --outSAMtype BAM SortedByCoordinate --outFileNamePrefix %s. \\
  --outSAMattrRGline %s \\
  --genomeDir %s \\
  --outTmpDir %s/star \\
  --readFilesIn %s %s && \\
#mv %s/%s.Aligned.sortedByCoord.out.bam %s && \\
/PUBLIC/software/RNA/samtools-1.2/samtools sort %s/%s.Aligned.sortedByCoord.out.bam %s/%s && \\
/PUBLIC/software/RNA/samtools-1.2/samtools index %s
''' % (self.alignDir,self.alignDir,self.sampleID,rg,self.databases['star_index'],self.alignDir,clean_fq1,clean_fq2,self.alignDir,self.sampleID,self.bam,self.alignDir,self.sampleID,self.alignDir,self.sampleID,self.bam)
		return cmd,order

	def pollution(self):
		rm_unmapped = 'rm %s\nrm %s' % (os.path.join(self.alignDir,self.sampleID+'.unmapped.fq.1.gz'),os.path.join(self.alignDir,self.sampleID+'.unmapped.fq.2.gz'))

		order = 'order pollution_%s after mapping_summary_%s' % (self.sampleID,self.sampleID)
		step1 = 'export PERL5LIB=/home/weiwenting/anaconda3/lib/site_perl/5.26.2/x86_64-linux-thread-multi/:$PERL5LIB\nperl %s %s 650 random.record > %s' % (os.path.join(self.softwares['script'],'Mapping/pollution/random_extract_fq.pl'),os.path.join(self.alignDir,self.sampleID+'.unmapped.fq.1.gz'),os.path.join(self.alignDir,'pollution','extracted.'+self.sampleID+'.unmapped.fq.1'))
		step2 = '%s %s 650 1e-1 random.record > %s' % (os.path.join(self.softwares['script'],'Mapping/pollution/generate_contamination_QC_sh.pl'),os.path.join(self.alignDir,'pollution',self.sampleID+'.unmapped.fq.1'),os.path.join(self.alignDir,'pollution','pollution.sh'))
		step3 = 'rate=`sed -n \'3p\' %s.stat | awk -F \'\\t\' \'{print $2}\' | awk -F \'%%\' \'{print $1}\' | awk -F \'(\' \'{print $2}\'`' % os.path.join(self.alignDir,self.sampleID)
		step_rm1 = 'rm extracted.%s.unmapped.fq.1' % self.sampleID
		step_rm2 = 'rm random.record'
		step_rm3 = 'rm extracted.%s.unmapped.fq.1.readxml' % self.sampleID
		step_rm4 = 'rm parsed.extracted.%s.unmapped.fq.1.readxml' % self.sampleID

		step4 = 'if [ $(echo "$rate < 85.0"|bc) -eq 1 ];then\n%s\n%s\ncd %s\nsh %s\n%s\n%s\n%s\n%s\nelse\n%s\nfi' % (step1,step2,os.path.join(self.alignDir,'pollution'),os.path.join(self.alignDir,'pollution','pollution.sh'),step_rm1,step_rm2,step_rm3,step_rm4,rm_unmapped)
		return step3+'\n'+step4+'\n',order


	def cal_readcount(self,libtype):
		order = 'order cal_readcount_%s after mapping_rnaseq_%s' % (self.sampleID,self.sampleID)
		if libtype=='fr-unstranded':
			s='0'
		elif libtype=='fr-firststrand':
			s='2'
		elif libtype=='fr-secondstrand':
			s='1'
		step1 = '\\\n\t'.join([self.softwares['featurecount'],
			'-T 8 -F GTF -t exon -g gene_id -s %s -Q 10 -C -p' % s,
			'-a %s' % self.databases['intron_gtf'],
			'-o %s' % os.path.join(self.alignDir,self.sampleID+'.readcount'),
			self.bam])
		step2 = '\\\n\t'.join([os.path.join(self.softwares['script'],'Mapping/regionstat'),
			'--count %s' % os.path.join(self.alignDir,self.sampleID+'.readcount'),
			'--summary %s.summary' % os.path.join(self.alignDir,self.sampleID+'.readcount'),
			'--outfile %s' % os.path.join(self.alignDir,self.sampleID+'.mapregion.xls')])
		step3 = '\\\n\t'.join([os.path.join(self.softwares['script'],'Mapping/mapregion.R'),
			'--region %s' % os.path.join(self.alignDir,self.sampleID+'.mapregion.xls'),
			'--prefix %s' % os.path.join(self.alignDir,self.sampleID),
			'--sampleid %s' % self.sampleID])
		step4 = 'rm %s' % os.path.join(self.alignDir,self.sampleID+'.readcount')
		
		return ' && \\\n'.join([step1,step2,step3,step4]),order

	def mapping_summary(self):
		order = 'order mapping_summary_%s after mapping_rnaseq_%s' % (self.sampleID,self.sampleID)
		step0 = "%s\ncd %s" % (self.softwares['mapping_summary_lib'],self.alignDir)
		step1 = '\\\n\t'.join(['perl %s' % os.path.join(self.softwares['script'],'Mapping/bam_stat.pl'),
			os.path.join(self.qcDir,self.sampleID+'.stat'),
			self.bam, os.path.join(self.alignDir,self.sampleID)])
		## insert size
		step2 = '\\\n\t'.join(['\njava -Xmx6G -jar %s' %self.softwares['picard'],
			'I=%s' % self.bam,
			'O=%s' % os.path.join(self.alignDir,self.sampleID+'.insertsize.tmp'),
			'H=%s' % os.path.join(self.alignDir,self.sampleID+'.insertsize.pdf')])
		step3 = '\\\n\t'.join(["\nsed -n '8p' %s" % os.path.join(self.alignDir,self.sampleID+'.insertsize.tmp'),
			"|cut -f 5-6 | awk '{print \"sample1\\t\"$0}' > %s" % os.path.join(self.alignDir,self.sampleID+'.insertsize.stat')])
		step4 = '\\\n\t'.join(["sed '1,12d' %s" % os.path.join(self.alignDir,self.sampleID+'.insertsize.tmp'),
			">%s" % os.path.join(self.alignDir,self.sampleID+'.insertsize.txt')])
		step5 = '\\\n\t'.join(['\n%s %s' % (self.softwares['R4.0.3'],os.path.join(self.softwares['script'],'Mapping/insertsize.R')),
			os.path.join(self.alignDir,self.sampleID+'.insertsize.txt'),
			os.path.join(self.alignDir,self.sampleID),
			self.sampleID])
		if self.sampleID.startswith('RS'):
			step6 = '\\\n\t'.join(['%s'%self.softwares['bamdst'],
				'-p %s/RS.merge.add.porbe.sort.merge.bed ' % self.softwares['BED'],
				'-o %s %s'%(self.alignDir,self.bam)])
		elif self.sampleID.startswith('RA'):
			step6 = '\\\n\t'.join(['%s'%self.softwares['bamdst'],
				'-p %s/RA.loci.bed ' % self.softwares['BED'],
				'-o %s %s'%(self.alignDir,self.bam)])
		elif self.sampleID.startswith('RT'):
			step6 = '\\\n\t'.join(['%s'%self.softwares['bamdst'],
				'-p %s/RT.bed ' % self.softwares['BED'],
				'-o %s %s'%(self.alignDir,self.bam)])    
		else:
			step6 = '\\\n\t'.join(['%s'%self.softwares['bamdst'],
				'-p %s '%(self.TR),
				'-o %s %s'%(self.alignDir,self.bam)])
		#step7 = 'cat %s %s >%s_qc_map.stat\n'%(os.path.join(self.qcDir,self.sampleID+'.stat'),os.path.join(self.alignDir,self.sampleID+'.stat'),os.path.join(self.alignDir,self.sampleID))
		step7 = 'python %s --analydir %s --sample %s --outfile %s_qc_map.stat\n' % (os.path.join(self.softwares['script'],'Mapping/merge.qc_mapping.py'),self.alignDir.split("Mapping")[0],self.sampleID,self.sampleID)
		return ' && \\\n'.join([step0,step1,step2,step3,step4,step5,step6,step7]),order

class Assembly:

	def __init__(self,sampleID,alignDir,assemDir,assemRootDir,softwares,databases):
		self.sampleID = sampleID
		self.alignDir = alignDir
		self.assemDir = assemDir
		self.assemRootDir = assemRootDir
		self.bam = os.path.join(alignDir, self.sampleID+'.bam')
		self.gtf = os.path.join(assemDir, self.sampleID+'.gtf')
		self.softwares = softwares
		self.databases = databases

	def assembly(self,soft='cufflinks',libtype='fr-firststrand'):
		order = 'order assembly_%s after mapping_report' % self.sampleID
		if soft == 'cufflinks':
			step1 = '\\\n\t'.join([self.softwares['cufflinks'],
				'-p 4 -u --library-type %s' % libtype,
				'-g %s' % self.databases['gtf'],
				'-o %s' % self.assemDir,
				self.bam])
			step2 = 'mv %s %s'%(os.path.join(self.assemDir,'transcripts.gtf'), self.gtf)
		elif soft == "stringtie":
			step1 = '\\\n\t'.join([self.softwares['stringtie'],
				'-p 4 -l %s' % self.sampleID,
				'-G %s' % self.databases['gtf'],
				'-o %s' % self.assemDir,
				self.bam])
			step2 = 'sed -i \'1,2d\' %s' % os.path.join(self.assemDir,self.sampleID+'.gtf')
		return ' && \\\n'.join([step1,step2]),order

	def merge_assembly(self,samples,soft='cufflinks',mode='transcript'):
		order = ['order merge_assembly after assembly_%s' % each  for each in samples]
		open(os.path.join(self.assemRootDir,'merge_gtf','gtf.list'),'w').write( \
			'\n'.join([os.path.join(self.assemRootDir,each,each+'.gtf') for each in samples])+'\n')
		step0 = 'export PATH=/PUBLIC/software/RNA/cufflinks-2.1.1:$PATH\n'
		step1 = '\\\n\t'.join([self.softwares['cuffmerge'],
			'-o %s' % os.path.join(self.assemRootDir,'merge_gtf'),
			'-p 4', 
			'-g %s' % self.databases['gtf'],
			'-s %s' % self.databases['fasta'],
			'%s' % os.path.join(self.assemRootDir,'merge_gtf','gtf.list')])
		step2 = '\\\n\t'.join([self.softwares['cuffcompare'],	 
			'-T -C -r %s' % self.databases['gtf'],
			'-R %s' % os.path.join(self.assemRootDir,'merge_gtf','merged.gtf'),
			'-o %s' % os.path.join(self.assemRootDir,'merge_gtf','compare')])
		step3 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['bin'],'compare_gtf.py'),
			self.databases['gtf'],
			os.path.join(self.assemRootDir,'merge_gtf','compare.combined.gtf'),
			os.path.join(self.assemRootDir,'merge_gtf')])
		if mode == 'gene':
			step4 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['bin'],'gff_merge_sort.py'),
				os.path.join(self.assemRootDir,'merge_gtf','novel_genes.gtf'),
				self.databases['gtf'],
				'> %s' % os.path.join(self.assemRootDir,'merge_gtf','known_novel.merge.gtf')])
		else:
			step4 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['bin'],'gff_merge_sort.py'),
				os.path.join(self.assemRootDir,'merge_gtf','novel_transcripts.gtf'),
				self.databases['gtf'],
				'> %s' % os.path.join(self.assemRootDir,'merge_gtf','known_novel.merge.gtf')])
		step5 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['bin'],'extractInfo.py'),
			os.path.join(self.assemRootDir,'merge_gtf','known_novel.merge.gtf'),
			'%s' % os.path.join(self.assemRootDir,'merge_gtf','gene_info.txt')])
		step6 = 'gzip -f %s' % ' '.join([os.path.join(self.assemRootDir,each,each+'.gtf') for each in samples])
		return step0+' && \\\n'.join([step1,step2,step3,step4,step5,step6]),order

	def novel_gene_annotation(self):
		order = 'order novel_gene_annotation after '
		step1 = ''
		return step1,order

	def cuffquant_quantification(self,libtype='fr-firststrand'):
		order = 'order cuffquant_quantification_%s after merge_assembly' %(self.sampleID)
		gtf_quant = os.path.join(self.assemRootDir,'merge_gtf','known_novel.merge.gtf')
		step1 = '\\\n\t'.join([self.softwares['cuffquant'],
			'--library-type %s -p 2 --max-bundle-frags 1165754' % libtype,
			'-o %s' % os.path.join(self.assemRootDir,'lncRNA_filter',self.sampleID+'.cuffquant'),
			gtf_quant, self.bam])
		return step1,order

	def expression_filter(self,samples,mode='transcript',libtype='fr-firststrand'):
		order = ['order expression_filter after cuffquant_quantification_%s' % each for each in samples]
		cxbs = [os.path.join(self.assemRootDir,'lncRNA_filter',self.sampleID+'.cuffquant','abundances.cxb') for each in samples]
		gtf_norm = os.path.join(self.assemRootDir,'merge_gtf','known_novel.merge.gtf')
		gtf4filter = os.path.join(self.assemRootDir,'merge_gtf','novel_transcripts.gtf')
		if mode == 'gene':
			gtf4filter = os.path.join(self.assemRootDir,'merge_gtf','novel_genes.gtf')
		step1 = '\\\n\t'.join([self.softwares['cuffnorm'],
			'--library-type %s --num-threads 2' % libtype,
			'-o %s' % os.path.join(self.assemRootDir,'lncRNA_filter','cuffnorm'),
			'-l %s' % ','.join(samples),
			gtf_norm] + cxbs)
		step2 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['bin'],'expression_filter.py'),
			'--exp_table %s' % os.path.join(self.assemRootDir,'lncRNA_filter','cuffnorm','isoforms.fpkm_table'),
			'--gtf %s' % gtf4filter,
			'--tr_info %s' % os.path.join(self.assemRootDir,'lncRNA_filter','assembled_tr_info.json'),
			'--output %s' % os.path.join(self.assemRootDir,'lncRNA_filter','novel.exp_filter.gtf')])
		step3 = 'rm %s' % gtf_norm
		step4 = 'rm -rf %s' % ' '.join([os.path.join(self.assemRootDir,'lncRNA_filter',self.sampleID+'.cuffquant') for each in samples])
		return ' && \\\n'.join([step1,step2,step3]),order
			
	def extract_fasta(self):
		order = 'order extract_fasta after merge_assembly'
		step1 = '\\\n\t'.join(['perl %s' % os.path.join(self.softwares['bin'],'extractTranscriptsfromFA.pl'),
			os.path.join(self.assemRootDir,'merge_gtf','novel_transcripts.gtf'),
			self.databases['fasta'],
			os.path.join(self.assemRootDir,'lncRNA_filter','novel.exp_filter.fasta')])
		return step1,order
	def split_fasta(self,fa_split_num=20):
		order = 'order split_fasta after extract_fasta'
		step1 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['bin'],'split_fa.py'),
			os.path.join(self.assemRootDir,'lncRNA_filter','novel.exp_filter.fasta'),
			'%s' % fa_split_num,
			os.path.join(self.assemRootDir,'lncRNA_filter','split')])
		return step1,order
	def lncRNA_cpc(self,eachfa,eachoutdir,n):
		order = 'order lncRNA_cpc_%s after split_fasta' % n
		step1 = '\\\n\t'.join(['perl %s' % os.path.join(self.softwares['bin'],'runCPC.pl'),
			'-seq %s' %  eachfa,
			'-db sp -strand plus',
			'-outdir %s' % eachoutdir])
		return step1,order

	def merge_cpc(self):
		order = ['order merge_cpc after %s' % each for each in ['lncRNA_cpc_'+str(elem) for elem in range(1,21)]]
		step1 = 'cat %s/lncRNA_filter/CPC/*/*cpc.txt|sed \'1itranscript_id\\ttranscript_length\\ttype\\tscore\' > %s/lncRNA_filter/CPC/CPC.result.xls' % (self.assemRootDir,self.assemRootDir)
		step2 = 'awk \'{if($3 == "noncoding"){print $1}}\' %s/lncRNA_filter/CPC/CPC.result.xls > %s/lncRNA_filter/CPC/CPC.noncoding.id.txt'%(self.assemRootDir,self.assemRootDir)
		step3 = 'awk \'{if($3 == "coding"){print $1}}\' %s/lncRNA_filter/CPC/CPC.result.xls > %s/lncRNA_filter/CPC/CPC.coding.id.txt\n'%(self.assemRootDir,self.assemRootDir)
		return ' && \\\n'.join([step1,step2,step3]),order

	def lncRNA_cnci(self,eachfa,eachoutdir,n):
		order = 'order lncRNA_cnci_%s after split_fasta' % n
		step1 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['bin'],'CNCI.py'),
			'-f %s' %  eachfa,
			'-m ve -p 1',
			'-o %s' % eachoutdir])
		return step1,order
	
	def merge_cnci(self):
		order = ['order merge_cnci after %s' % each for each in ['lncRNA_cnci_'+str(elem) for elem in range(1,21)]]	
		step1 = 'cat %s/lncRNA_filter/CNCI/*/CNCI.index > %s/lncRNA_filter/CNCI/CNCI.result.xls' % (self.assemRootDir,self.assemRootDir)
		step2 = 'awk \'{if($2 == "noncoding"){print $1}}\' %s/lncRNA_filter/CNCI/CNCI.result.xls > %s/lncRNA_filter/CNCI/CNCI.noncoding.id.txt'%(self.assemRootDir,self.assemRootDir)
		step3 = 'awk \'{if($2 == "coding"){print $1}}\' %s/lncRNA_filter/CNCI/CNCI.result.xls > %s/lncRNA_filter/CNCI/CNCI.coding.id.txt\n'%(self.assemRootDir,self.assemRootDir)
		return ' && \\\n'.join([step1,step2,step3]),order

	def lncRNA_pfam(self,eachfa,eachoutdir,n):
		order = 'order lncRNA_pfam_%s after split_fasta' % n
		step1 = 'PATH=/PUBLIC/software/public/Annotation/hmmer-3.1b1/bin:$PATH\nexport PERL5LIB=/PUBLIC/software/RNA/perl5/lib/perl5:/PUBLIC/software/RNA/perl5/lib/perl5/x86_64-linux:/PUBLIC/software/RNA/perl5/lib/perl5/x86_64-linux-thread-multi:/PUBLIC/software/public/System/Perl-5.18.2/lib/perl5/5.18.2/:$SOFT/System/Perl-5.18.2/lib/perl5/site_perl/5.18.2/x86_64-linux/:/PUBLIC/software/RNA/PfamScan\n'
		step2 = 'perl %s %s > %s' % (os.path.join(self.softwares['bin'],'seq2protein.pl'),eachfa,os.path.join(self.assemRootDir,'lncRNA_filter','PFAM',eachoutdir,'protein.fa'))
		step3 = '\\\n\t'.join(['perl %s' % os.path.join(self.softwares['bin'],'pfam_scan.pl'),
			'-fasta %s' %  os.path.join(self.assemRootDir,'lncRNA_filter','PFAM',eachoutdir,'protein.fa'),
			'-dir %s' % self.databases['pfam_hmm'],
			'-outfile %s' % os.path.join(eachoutdir,'pfam_scan.out'),
			'-pfamB -cpu 2'])
		return step1+' && \\\n'.join([step2,step3]),order

	def merge_pfam(self):
		order = ['order merge_pfam after %s' % each for each in ['lncRNA_pfam_'+str(elem) for elem in range(1,21)]]
		step1 = 'cat %s/lncRNA_filter/PFAM/*/pfam_scan.out > %s/lncRNA_filter/PFAM/pfam_scan.out' % (self.assemRootDir,self.assemRootDir)
		step2 = 'cat %s/lncRNA_filter/PFAM/pfam_scan.out |sed -e \'/^#/d\' -e \'/^$/d\' |awk \'OFS="\\t"{print ($1,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15)}\'|sed \'1iseq id\\thmm acc\\thmm name\\ttype\\thmm start\\thmm end\\thmm length\\tbit score\\tE-value\\tsignificance\tclan\' > %s/lncRNA_filter/PFAM/PFAM.result.xls'%(self.assemRootDir,self.assemRootDir)
		step3 = 'python %s %s %s %s\n' % (os.path.join(self.softwares['bin'],'select_pfam_id.py'),os.path.join(self.assemRootDir,'merge_gtf','novel_transcripts.gtf'),os.path.join(self.assemRootDir,'lncRNA_filter','PFAM','PFAM.result.xls'),os.path.join(self.assemRootDir,'lncRNA_filter','PFAM'))
		return ' && \\\n'.join([step1,step2,step3]),order

	def lncRNA_identification(self):
		order = ['order lncRNA_identification after %s' % each for each in ['merge_cpc','merge_cnci','merge_pfam']]
		step1 = '\\\n\t'.join(['perl %s' % os.path.join(self.softwares['bin'],'CPC_CNCI_PFAM_CSF.venn_v2.pl'),
			'-gtf %s' % os.path.join(self.assemRootDir,'merge_gtf','novel_transcripts.gtf'),
			'-outdir %s' % os.path.join(self.assemRootDir,'lncRNA_filter'),
			'-CNCI %s' % os.path.join(self.assemRootDir,'lncRNA_filter','CNCI','CNCI.noncoding.id.txt'),
			'-CPC %s' % os.path.join(self.assemRootDir,'lncRNA_filter','CPC','CPC.noncoding.id.txt'),
			'-PFAM %s' % os.path.join(self.assemRootDir,'lncRNA_filter','PFAM','PFAM.noncoding.id.txt')])
		step2 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['bin'],'get_new_gtf.py'),
			os.path.join(self.assemRootDir,'lncRNA_filter','noncoding.result.id'),
			os.path.join(self.assemRootDir,'lncRNA_filter','coding.result.id'),
			os.path.join(self.assemRootDir,'merge_gtf','novel_transcripts.gtf'),
			os.path.join(self.assemRootDir,'lncRNA_filter')])
		step3 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['bin'],'lncRNA_classify.py'),
			'--lnc_gtf %s' % os.path.join(self.assemRootDir,'lncRNA_filter','Novel_lncRNA.gtf'),
			'--out_dir %s' % os.path.join(self.assemRootDir,'lncRNA_filter')])
		step4 = '\\\n\t'.join(['Rscript %s' % os.path.join(self.softwares['bin'],'lncRNA_classify.R'),
			'%s %s' % (os.path.join(self.assemRootDir,'lncRNA_filter','lncRNA_classification.txt'),os.path.join(self.assemRootDir,'lncRNA_filter'))])
		step5 = 'cat %s %s > %s' % (self.databases['lncRNA_gtf'],os.path.join(self.assemRootDir,'lncRNA_filter','Novel_lncRNA.gtf'),os.path.join(self.assemRootDir,'lncRNA_filter','lncRNA.gtf'))
		step6 = 'cat %s %s > %s\n' % (self.databases['mRNA_gtf'],os.path.join(self.assemRootDir,'lncRNA_filter','Novel_mRNA.gtf'),os.path.join(self.assemRootDir,'lncRNA_filter','mRNA.gtf'))
		return ' && \\\n'.join([step1,step2,step3,step4,step5,step6]),order

	def lncRNA_signatures(self):
		order = 'order lncRNA_signatures after lncRNA_identification'
		step1 = '\\\n\t'.join(['perl %s' % os.path.join(self.softwares['bin'],'extractTranscriptsfromFA.pl'),
			'%s' % os.path.join(self.assemRootDir,'lncRNA_filter','Novel_lncRNA.gtf'),
			'%s' % self.databases['fasta'],
			'%s' % os.path.join(self.assemRootDir,'lncRNA_filter','Novel_lncRNA.fa')])
		step2 = '\\\n\t'.join([self.softwares['getorf'],
			'-sequence %s' % os.path.join(self.assemRootDir,'lncRNA_filter','Novel_lncRNA.fa'),
			'-outseq %s' % os.path.join(self.assemRootDir,'lncRNA_filter','Novel_lncRNA.orf.fa')])
		step3 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['bin'],'get_transcript_info.py'),
			'--tr_type Novel_lncRNA',
			'--gtf %s' % os.path.join(self.assemRootDir,'lncRNA_filter','Novel_lncRNA.gtf'),
			'--orf_info %s' % os.path.join(self.assemRootDir,'lncRNA_filter','Novel_lncRNA.orf.fa'),
			'--out_dir %s' % os.path.join(self.assemRootDir,'lncRNA_filter')])
		step4 = '\\\n\t'.join(['perl %s' % os.path.join(self.softwares['bin'],'extractTranscriptsfromFA.pl'),
			'%s' % os.path.join(self.assemRootDir,'lncRNA_filter','Novel_mRNA.gtf'),
			'%s' % self.databases['fasta'],
			'%s' % os.path.join(self.assemRootDir,'lncRNA_filter','Novel_mRNA.fa')])
		step5 = '\\\n\t'.join([self.softwares['getorf'],
			'-sequence %s' % os.path.join(self.assemRootDir,'lncRNA_filter','Novel_mRNA.fa'),
			'-outseq %s' % os.path.join(self.assemRootDir,'lncRNA_filter','Novel_mRNA.orf.fa')])
		step6 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['bin'],'get_transcript_info.py'),
			'--tr_type Novel_mRNA',
			'--gtf %s' % os.path.join(self.assemRootDir,'lncRNA_filter','Novel_mRNA.gtf'),
			'--orf_info %s' % os.path.join(self.assemRootDir,'lncRNA_filter','Novel_mRNA.orf.fa'),
			'--out_dir %s' % os.path.join(self.assemRootDir,'lncRNA_filter')])
		step7 = 'ln -s %s %s' % (self.databases['Annotated_lncRNA_information.xls'],os.path.join(self.assemRootDir,'lncRNA_filter','Annotated_lncRNA_information.xls'))
		step8 = 'ln -s %s %s' % (self.databases['Annotated_lncRNA_transcript_info.json'],os.path.join(self.assemRootDir,'lncRNA_filter','Annotated_lncRNA_transcript_info.json'))
		step9 = 'ln -s %s %s' % (self.databases['Annotated_mRNA_information.xls'],os.path.join(self.assemRootDir,'lncRNA_filter','Annotated_mRNA_information.xls'))
		step10 = 'ln -s %s %s' % (self.databases['Annotated_mRNA_transcript_info.json'],os.path.join(self.assemRootDir,'lncRNA_filter','Annotated_mRNA_transcript_info.json'))
		step11 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['bin'],'merge_transcript_info.py'),
		'--tr_info_dir %s' % os.path.join(self.assemRootDir,'lncRNA_filter')])
		step12 = '\\\n\t'.join(['Rscript %s' % os.path.join(self.softwares['bin'],'lncRNA_feature_plot.R'),
			os.path.join(self.assemRootDir,'lncRNA_filter','lncRNA_features.txt'),
			'lncRNA',
			os.path.join(self.assemRootDir,'lncRNA_filter')])
		
		return ' && \\\n'.join([step1,step2,step3,step4,step5,step6,step7,step8,step9,step10,step11,step12]),order

class Quantification:

	def __init__(self,sampleID,qcDir,alignDir,expDir,softwares,databases,dic_disease):
		self.sampleID = sampleID
		self.qcDir = qcDir
		self.alignDir = alignDir
		self.alignRoot = os.path.dirname(alignDir)
		self.expDir = expDir
		self.expRoot = os.path.dirname(expDir)
		self.Root=os.path.dirname(os.path.dirname(expDir))
		self.bam = os.path.join(alignDir, self.sampleID+'.bam')
		self.softwares = softwares
		self.databases = databases
		self.dic_disease = dic_disease
	def quant_prepare(self,assembly=True,assemdir='',analysis_type='RNAseq',sample_lis=[]):
		if assembly:
			if analysis_type == 'lncRNA':
				order = 'order quant_prepare after lncRNA_identification'
				step1 = 'cat %s %s %s > %s' %(os.path.join(assemdir,'lncRNA_filter','Novel_lncRNA.gtf'),os.path.join(assemdir,'lncRNA_filter','Novel_mRNA.gtf'),os.path.join(assemdir,'lncRNA_filter','Unclassified.gtf'),os.path.join(self.expRoot,'tmp_all_novel.gtf'))
				step2 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['bin'],'gff_merge_sort.py'),
					os.path.join(self.expRoot,'tmp_all_novel.gtf'),
					self.databases['gtf'],
					'> %s' % os.path.join(self.expRoot,'all_transcripts.gtf')])
				step3 = 'rm %s' % os.path.join(self.expRoot,'tmp_all_novel.gtf')
				step4 = 'ln -s %s %s' % (os.path.join(assemdir,'lncRNA_filter','lncRNA.gtf'),os.path.join(self.expRoot,'lncRNA.gtf'))
				step5 = 'ln -s %s %s' % (os.path.join(assemdir,'lncRNA_filter','mRNA.gtf'),os.path.join(self.expRoot,'mRNA.gtf'))
				step6 = 'ln -s %s %s' % (os.path.join(assemdir,'lncRNA_filter','Unclassified.gtf'),os.path.join(self.expRoot,'Unclassified.gtf'))
				step7 = '\\\n\t'.join(['perl %s '% os.path.join(self.softwares['script'],'get_gene_transcript_id.pl'),
					'%s' % os.path.join(self.expRoot,'lncRNA.gtf'),
					'%s' % os.path.join(self.expRoot,'lncRNA_gene.id'),
					'%s' % os.path.join(self.expRoot,'lncRNA_transcript.id')])
				step8 = '\\\n\t'.join(['perl %s '% os.path.join(self.softwares['script'],'get_gene_transcript_id.pl'),
					'%s' % os.path.join(self.expRoot,'mRNA.gtf'),
					'%s' % os.path.join(self.expRoot,'mRNA_gene.id'),
					'%s' % os.path.join(self.expRoot,'mRNA_transcript.id')])
				step9 = '\\\n\t'.join(['perl %s '% os.path.join(self.softwares['script'],'get_gene_transcript_id.pl'),
					'%s' % os.path.join(self.expRoot,'Unclassified.gtf'),
					'%s' % os.path.join(self.expRoot,'Unclassified_gene.id'),
					'%s' % os.path.join(self.expRoot,'Unclassified_transcript.id')])

                		step10 = 'python %s %s %s %s' % (os.path.join(self.softwares['script'],'Quantification/gene_description.py'),os.path.join(self.expRoot,'all_transcripts.gtf'),os.path.join(self.expRoot,'trans.description'),os.path.join(self.expRoot,'gene.description'))

				return ' && \\\n'.join([step1,step2,step3,step4,step5,step6,step7,step8,step9,step10]),order
			elif analysis_type == 'RNAseq':
				order = 'order quant_prepare after merge_assembly'
				step1 = 'ln -s %s %s' % (os.path.join(assemdir,'merge_gtf','known_novel.merge.gtf'),os.path.join(self.expRoot,'all_transcripts.gtf'))
				step2 = 'python %s %s %s %s' % (os.path.join(self.softwares['script'],'Quantification/gene_description.py'),os.path.join(self.expRoot,'all_transcripts.gtf'),os.path.join(self.expRoot,'trans.description'),os.path.join(self.expRoot,'gene.description'))
				return ' && \\\n'.join([step1,step2]),order
		else:
			order = 'order quant_prepare after mapping_report'
			step1 = 'ln -s %s %s' % (self.databases['gtf'],os.path.join(self.expRoot,'all_transcripts.gtf'))
			step2 = 'python %s %s %s %s' % (os.path.join(self.softwares['script'],'Quantification/gene_description.py'),os.path.join(self.expRoot,'all_transcripts.gtf'),os.path.join(self.expRoot,'trans.description'),os.path.join(self.expRoot,'gene.description'))
			return ' && \\\n'.join([step1,step2+'\n']),order
	def construct_fq(self):
		fq1_gz = os.path.join(self.qcDir,'_'.join([self.sampleID,'1.clean.fq.gz']))
		fq2_gz = os.path.join(self.qcDir,'_'.join([self.sampleID,'2.clean.fq.gz']))
		return fq1_gz,fq2_gz

	def prepare_rsem_index(self,map='bowtie2'):
		order = 'order prepare_rsem_index after quant_prepare'
		if map == 'star':
			step1 = 'cd %s' % os.path.join(self.expRoot,'RSEM_INDEX')
			step2 = 'ln -s %s %s' % (self.databases['fasta'],os.path.join(self.expRoot,'RSEM_INDEX','genome.fa'))
			step3 = 'ln -s %s %s' % (os.path.join(self.expRoot,'all_transcripts.gtf'),os.path.join(self.expRoot,'RSEM_INDEX','all_transcripts.gtf'))
			step4 = ' '.join(['perl %s' % os.path.join(self.softwares['bin'],'getid_transcript.pl'),
				'all_transcripts.gtf','all_transcripts.txt'])
			step5 = ' '.join(['/PUBLIC/software/RNA/RSEM/RSEM-1.3.0//rsem-prepare-reference',
				'--gtf all_transcripts.gtf',
				'--transcript-to-gene-map all_transcripts.txt',
				'--star --star-path /PUBLIC/software/RNA/STAR/STAR2.5/bin/Linux_x86_64 -p 8',
				'genome.fa all_transcripts','\n'])
			return ' && \\\n'.join([step1,step2,step3,step4,step5]),order
		elif map == 'bowtie2':
			step1 = 'cd %s' % os.path.join(self.expRoot,'RSEM_INDEX')
			step2 = 'ln -s %s %s' % (self.databases['fasta'],os.path.join(self.expRoot,'RSEM_INDEX','genome.fa'))
			step3 = 'ln -s %s %s' % (os.path.join(self.expRoot,'all_transcripts.gtf'),os.path.join(self.expRoot,'RSEM_INDEX','all_transcripts.gtf'))
			step4 = ' '.join(['perl %s' % os.path.join(self.softwares['bin'],'getid_transcript.pl'),
				'all_transcripts.gtf','all_transcripts.txt'])
			step5 = ' '.join(['/PUBLIC/software/RNA/RSEM/RSEM-1.3.0//rsem-prepare-reference',
				'--gtf all_transcripts.gtf',
				'--transcript-to-gene-map all_transcripts.txt',
				'--bowtie2 --bowtie2-path /PUBLIC/software/public/Alignment/bowtie2-2.0.6 -p 8 ',
				'genome.fa all_transcripts','\n'])
			return ' && \\\n'.join([step1,step2,step3,step4,step5]),order
		
	def rsem_quantification(self,assembly=False,libtype='fr-firststrand',map='bowtie2'):
		strand = {'fr-firststrand':'0','fr-unstranded':'0.5','fr-secondstrand':'1'}

		order = ['order quantification_%s after quant_prepare' % self.sampleID]
		if map=='bowtie2':
			rsem_index = self.databases['rsem_index_bowtie2']
		else:
			rsem_index = self.databases['rsem_index_star']
		if assembly:
			order += ['order quantification_%s after prepare_rsem_index' % self.sampleID]
			rsem_index = os.path.join(self.expRoot,'RSEM_INDEX','all_transcripts')
		clean_fq1,clean_fq2 = self.construct_fq()
		opt = '--bowtie2 --bowtie2-path %s' % self.softwares['bowtie2']
		if map == 'star':
			opt = '--star --star-path %s --star-gzipped-read-file' % self.softwares['star_bin']
		step0 = 'cd %s\n' % self.expDir
		step1 = '\\\n\t'.join([self.softwares['rsem_exp'],
			'-p 8 --paired-end --estimate-rspd --forward-prob %s --output-genome-bam --sort-bam-by-coordinate --time' % strand[libtype],
			opt,
			'--temporary-folder %s' % os.path.join(self.expDir,'temp'),
			clean_fq1,
			clean_fq2,
			rsem_index,
			os.path.join(self.expDir,self.sampleID)])
		step2 = 'cut -f1,5  %s.genes.results|sed \'1d\' > %s.genes.readcount' % (os.path.join(self.expDir,self.sampleID),os.path.join(self.expDir,self.sampleID))
		step3 = 'cut -f1,7  %s.genes.results|sed \'1d\' > %s.genes.fpkm' % (os.path.join(self.expDir,self.sampleID),os.path.join(self.expDir,self.sampleID))
		step4 = 'cut -f1,5  %s.isoforms.results|sed \'1d\' > %s.transcripts.readcount' % (os.path.join(self.expDir,self.sampleID),os.path.join(self.expDir,self.sampleID))
		step5 = 'cut -f1,7  %s.isoforms.results|sed \'1d\' > %s.transcripts.fpkm' % (os.path.join(self.expDir,self.sampleID),os.path.join(self.expDir,self.sampleID))
		step6 = 'rm %s' % os.path.join(self.expDir,self.sampleID+'.transcript.bam')
		return step0 + ' && \\\n'.join([step1,step2,step3,step4,step5,step6]),order

	def htseq_quantification(self,libtype='fr-firststrand'):
		order = 'order quantification_%s after quant_prepare' % self.sampleID
		strand = {'fr-firststrand':'reverse','fr-unstranded':'no','fr-secondstrand':'yes'}
		sort_bam = os.path.join([])	
		step1 = '%s\nln -s %s %s' % (self.softwares['perl_library'],self.bam,os.path.join(self.expDir,self.sampleID+'.bam'))
		step2 = '\\\n\t'.join(['python -m HTSeq.scripts.count -m union -s %s -t exon -f bam' % strand[libtype],
			'%s' % os.path.join(self.expDir,self.sampleID+'.bam'),
			'%s' % os.path.join(self.expRoot,'all_transcripts.gtf'),
			' > %s' % os.path.join(self.expDir,self.sampleID+'.genes.readcount')])
		step3 = '\\\n\t'.join(['perl %s' % os.path.join(self.softwares['script'],'Quantification/calRPKM.pl'),
			'%s' % os.path.join(self.expDir,self.sampleID+'.genes.readcount'),
			'%s' % os.path.join(self.expRoot,'all_transcripts.gtf'),
			'%s' % os.path.join(self.expDir,self.sampleID+'.genes.fpkm') ])
		step4 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['script'],'Quantification/anno_fpkm_sample.py'),
			'--fpkm %s' % os.path.join(self.expRoot,self.sampleID,self.sampleID+'.genes.fpkm'),
			'--genedes %s' % os.path.join(self.expRoot,'gene.description'),
			'--outfile %s' % os.path.join(self.expRoot,self.sampleID,self.sampleID+'.anno_genes.FPKM.xls'),
			'--sample %s' % self.sampleID])
		step5 = '\\\n\t'.join(['kepping=`cat %s|grep ABL1|awk \'{print $NF}\'`\npython %s' % (os.path.join(self.expRoot,self.sampleID,self.sampleID+'.anno_genes.FPKM.xls'),os.path.join(self.softwares['script'],'Quantification/get_FPKM_Ratio.py')),
			'--sample %s' % self.sampleID,
			'--fpkm %s' % os.path.join(self.expRoot,self.sampleID,self.sampleID+'.anno_genes.FPKM.xls'),
			'--mean %s' % os.path.join(self.softwares['script'],'Quantification/FPKM_Ratio/'),
			'--panel %s' % (self.dic_disease[self.sampleID]),
			'--outfile %s' % (os.path.join(self.expRoot,self.sampleID,self.sampleID+'.anno_genes.FPKM_ratio'))])
		step6 = 'cd %s\nperl %s %s %s %s %s %s %s' % (os.path.join(self.expRoot,self.sampleID),os.path.join(self.softwares['script'],'Quantification/excel.combine.pl'),self.sampleID+'.anno_genes.FPKM_ratio.ABL1',self.sampleID+'.anno_genes.FPKM_ratio.GAPDH',self.sampleID+'.anno_genes.FPKM_ratio.ABL1.filter',self.sampleID+'.anno_genes.FPKM_ratio.GAPDH.filter',self.sampleID+'.anno_genes.FPKM_ratio.overlap.filter',self.sampleID+'.FPKM_ratio.out.xlsx')
		step7 = 'cat %s|grep CRLF2|cut -f 1-7|awk \'{if($NF>0.15) print $0"\\t表达增高";else if($NF<0.15) print $0"\\t表达正常"}\' > %s.CRLF2.FPKM.txt' % (os.path.join(self.expRoot,self.sampleID,self.sampleID+'.anno_genes.FPKM_ratio.ABL1'),self.sampleID)

		return ' && \\\n'.join([step1,step2,step3,"echo no analysis\n"]),order


	def htseq_quantification_FPKM_ratio(self,libtype='fr-firststrand'):
		order = 'order quantification_FPKM_ratio_%s after quantification_%s' % (self.sampleID,self.sampleID)
		### 20250923 基线分开用各自对应的基线文件
		if self.sampleID.startswith("RS") or self.sampleID.startswith("RH") or self.sampleID.startswith("RA") or self.sampleID.startswith("RT"):
			get_FPKM_ratio = os.path.join(self.softwares['script'],'Quantification/FPKM_Ratio/targetRNA.get_FPKM_ratio.py')
		else:
			get_FPKM_ratio = os.path.join(self.softwares['script'],'Quantification/FPKM_Ratio/mRNA.get_FPKM_ratio.py')
		
		step4 = '\\\n\t'.join(['%s\ncd %s\npython %s' % (self.softwares['perl_library'],os.path.join(self.expRoot,self.sampleID),os.path.join(self.softwares['script'],'Quantification/anno_fpkm_sample.py')),
			'--fpkm %s' % os.path.join(self.expRoot,self.sampleID,self.sampleID+'.genes.fpkm'),
			'--genedes %s' % os.path.join(self.expRoot,'gene.description'),
			'--outfile %s' % os.path.join(self.expRoot,self.sampleID,self.sampleID+'.anno_genes.FPKM.xls'),
			'--sample %s' % self.sampleID])
		step5 = '\\\n\t'.join(['python %s' % (get_FPKM_ratio),
			'--sample %s' % self.sampleID,
			'--fpkm %s' % os.path.join(self.expRoot,self.sampleID,self.sampleID+'.anno_genes.FPKM.xls'),
			'--mean %s' % os.path.join(self.softwares['script'],'Quantification/FPKM_Ratio/'),
			'--panel %s' % (self.dic_disease[self.sampleID]),
			'--outfile %s' % (os.path.join(self.expRoot,self.sampleID,self.sampleID+'.anno_genes.FPKM_ratio'))])
		### add step6
		step6 = '\n###extract specific gene\ncd %s\npython %s -g %s -i %s -o %s' % (os.path.join(self.expRoot,self.sampleID),os.path.join(self.softwares['script'],'Quantification/extract.specific.gene.py'),os.path.join(self.softwares['script'],'Quantification/specific.genelist'),self.sampleID+'.anno_genes.FPKM_ratio.ABL1',self.sampleID+'.anno_genes.FPKM_ratio.specificGene')
		step7 = 'perl %s %s %s %s' % (os.path.join(self.softwares['script'],'Quantification/excel.combine.pl'),self.sampleID+'.anno_genes.FPKM_ratio.ABL1',self.sampleID+'.anno_genes.FPKM_ratio.specificGene',self.sampleID+'.FPKM_ratio.out.xlsx')
		#step7 = '\ncd %s\nperl %s %s %s %s %s %s %s' % (os.path.join(self.expRoot,self.sampleID),os.path.join(self.softwares['script'],'Quantification/excel.combine.pl'),self.sampleID+'.anno_genes.FPKM_ratio.ABL1',self.sampleID+'.anno_genes.FPKM_ratio.GAPDH',self.sampleID+'.anno_genes.FPKM_ratio.ABL1.filter',self.sampleID+'.anno_genes.FPKM_ratio.GAPDH.filter',self.sampleID+'.anno_genes.FPKM_ratio.overlap.filter',self.sampleID+'.FPKM_ratio.out.xlsx')
		step8 = '\n###CRLF2(0.15) \ncat %s|grep gene_id > %s.CRLF2.FPKM.txt \ncat %s|grep CRLF2|cut -f 1-7|awk \'{if($NF>0.15) print $0"\\t表达增高";else if($NF<0.15) print $0"\\t表达正常"}\' >> %s.CRLF2.FPKM.txt' % (os.path.join(self.expRoot,self.sampleID,self.sampleID+'.anno_genes.FPKM_ratio.ABL1'),self.sampleID,os.path.join(self.expRoot,self.sampleID,self.sampleID+'.anno_genes.FPKM_ratio.ABL1'),self.sampleID)
		step9 = 'perl %s %s %s' % (os.path.join(self.softwares['script'],'Quantification/convert.txt_to_xlsx.pl'),self.sampleID+'.CRLF2.FPKM.txt',self.sampleID+'.CRLF2.FPKM.xlsx')
		#### hemoglobin
		step10 = '\n##calculate hemoglobin percentage\npython %s --fpkm %s.anno_genes.FPKM.xls --readcount %s.genes.readcount --out %s.hemoglobin.ratio.xls'% (os.path.join(self.softwares['script'],'Quantification/count.hemoglobin.ratio.py'),self.sampleID,self.sampleID,self.sampleID)    
	
		### 20250923 RH:只生成CRLF2.FPKM.txt； 非RH生成FPKM_ratio.out.xlsx(5个sheet改成2个sheet)
		if self.sampleID.startswith("RH"):
			return ' && \\\n'.join([step4,step5,step8,step9,step10,"echo done\n"]),order
		else:
			return ' && \\\n'.join([step4,step5,step6,step7,step10,"echo done\n"]),order



	def cuffquant(self,libtype='fr-firststrand'):
		order = 'order cuffquant_%s after quant_prepare' % self.sampleID
		step1 = ' '.join([self.softwares['cuffquant'],
			'--library-type %s -p 2' % libtype,
			'-o %s' % self.expDir,
			'%s' % os.path.join(self.expRoot,'all_transcripts.gtf'),
			self.bam])
		return step1,order
	def cuffdiff_quantification(self,sample_lis,groupname,groups,libtype='fr-firststrand',assembly=True,analysis_type="RNAseq",repeat_flag=True):
		order = ['order cuffdiff_quantification after quant_prepare']
		gtf = self.databases['gtf']
		samples=','.join(sample_lis)

		if assembly:
			gtf = os.path.join(self.expRoot,'all_transcripts.gtf')
				
		if repeat_flag:
			order += ['order cuffdiff_quantification after cuffquant_%s' % each for each in sample_lis]
			cxb_tmp=[]
			cxbs=[]
			for each in [elem.split(',') for elem in [each.replace(':',',') for each in groups]]:
				for every in each:
					cxb_tmp.append(os.path.join(self.expRoot,every,'abundances.cxb'))
				cxbs.append(','.join(cxb_tmp))
				cxb_tmp=[]
			step1 = '\\\n\t'.join([self.softwares['cuffdiffv2'],
				'-o %s' % os.path.join(self.expRoot,'cuffdiff'),
				'--labels %s ' % groupname,
				'%s' % gtf,
				'%s' % '\\\n\t'.join(cxbs),
				'--FDR 0.05 --library-type %s --library-norm-method geometric --dispersion-method pooled --num-threads 2 ' % libtype])
			return step1,order
		else :
			order += ['order cuffdiff_quantification after mapping_report']
			bam_tmp=[]
			bams=[]
			for each in [elem.split(',') for elem in [each.replace(':',',') for each in groups]]:
				for every in each:
					bam_tmp.append(os.path.join(self.alignRoot,every,every+'.bam'))
				bams.append(','.join(bam_tmp))
				bam_tmp=[]
			step1 = '\\\n\t'.join([self.softwares['cuffdiffv1'],
				'-o %s' % os.path.join(self.expRoot,'cuffdiff'),
				'--labels %s' % groupname,
				'%s' % gtf,
				'%s' % '\\\n\t'.join(bams),
				' --FDR 0.05 --library-type %s --num-threads 8' % libtype ])
			return step1,order
	def quantification(self,gtf=None,assembly=False,soft='htseq',analysis_type="RNAseq",libtype='fr-firststrand',repeat_flag=''):
		gtf_quant = self.databases['gtf']
		if gtf:
			gtf_quant = gtf
		if soft == 'htseq':
			script,order = self.htseq_quantification(libtype=libtype)
		elif soft == 'rsem':
			script,order = self.rsem_quantification(assembly=assembly,libtype=libtype)
		elif soft == 'cuffdiff' and repeat_flag == True:
			script,order = self.cuffquant()
		return script,order

	def quantification_summary(self,sample_lis,groups,groupname,mod='gene',seqtype='RNAseq',exp_soft='cuffdiff',repeat_flag=True):
		order = ['order quantification_summary after quantification_%s'%(each) \
			for each in sample_lis]
		step_lis=[]
		if exp_soft == 'cuffdiff':
			if repeat_flag:
				replication='replication'
			else:
				replication='non_replication'
			order = 'order quantification_summary after cuffdiff_quantification'
			step1 = '\\\n\t'.join(['python %s ' % os.path.join(self.softwares['script'],'quantification_'+replication+'.py'),
				'--out-dir %s' % self.expRoot,
				'--Quantification_dir %s' % os.path.join(self.expRoot,'cuffdiff'),
				'--labels %s' % ','.join(sample_lis),
				'--groups %s' % groups,
				'--groupnames %s' % groupname])
			step_lis.append(step1)
			step1_add1='\\\n\t'.join(['python %s' % os.path.join(self.softwares['script'],'Quantification/anno_fpkm.py'),
				'--fpkm %s' % os.path.join(self.expRoot,'Summary','genes.FPKM.xls'),
				'--genedes %s' % os.path.join(self.expRoot,'gene.description'),
				'--outfile %s' % os.path.join(self.expRoot,'Summary','anno_genes.FPKM.xls')])
			step_lis.append(step1_add1)
			step1_add2='\\\n\t'.join(['python %s' % os.path.join(self.softwares['script'],'Quantification/anno_fpkm.py'),
				'--fpkm %s' % os.path.join(self.expRoot,'Summary','transcripts.FPKM.xls'),
				'--genedes %s' % os.path.join(self.expRoot,'trans.description'),
				'--outfile %s' % os.path.join(self.expRoot,'Summary','anno_transcripts.FPKM.xls'),
				'--tr'])
			step_lis.append(step1_add2)

		else:
			if mod == 'gene' or mod == 'both':
				cmd_genereadcount = ','.join([os.path.join(self.expRoot,eachsample,eachsample+'.genes.readcount') for eachsample in sample_lis])
				cmd_genefpkm = ','.join([os.path.join(self.expRoot,eachsample,eachsample+'.genes.fpkm') for eachsample in sample_lis])
				step1 = '\\\n\t'.join(['\npython %s %s gene %s' % (os.path.join(self.softwares['script'],'Quantification/get_merge.py'),cmd_genereadcount,os.path.join(self.expRoot,'Summary','genes.readcount.xls'))])
				step_lis.append(step1)
				step2 = '\\\n\t'.join(['\npython %s %s gene %s' % (os.path.join(self.softwares['script'],'Quantification/get_merge.py'),cmd_genefpkm,os.path.join(self.expRoot,'Summary','genes.FPKM.xls'))])
				step_lis.append(step2)
				step2_add = '\\\n\t'.join(['\npython %s' % os.path.join(self.softwares['script'],'Quantification/anno_fpkm.py'),
					'--fpkm %s' % os.path.join(self.expRoot,'Summary','genes.FPKM.xls'),
					'--genedes %s' % os.path.join(self.expRoot,'gene.description'),
					'--outfile %s' % os.path.join(self.expRoot,'Summary','anno_genes.FPKM.xls')])
				step_lis.append(step2_add)
			if mod == 'isoform' or mod == 'both':
				cmd_isoformreadcount = ','.join([os.path.join(self.expRoot,eachsample,eachsample+'.transcripts.readcount') for eachsample in sample_lis])
				cmd_isoformfpkm = ','.join([os.path.join(self.expRoot,eachsample,eachsample+'.transcripts.fpkm') for eachsample in sample_lis])
				step1 = '\\\n\t'.join(['\npython %s %s transcript %s' % (os.path.join(self.softwares['script'],'Quantification/get_merge.py'),cmd_isoformreadcount,os.path.join(self.expRoot,'Summary','transcripts.readcount.xls'))])
				step_lis.append(step1)
				step2 = '\\\n\t'.join(['\npython %s %s transcript %s' % (os.path.join(self.softwares['script'],'Quantification/get_merge.py'),cmd_isoformfpkm,os.path.join(self.expRoot,'Summary','transcripts.FPKM.xls'))])
				step_lis.append(step2)
				step2_add = '\\\n\t'.join(['\npython %s' % os.path.join(self.softwares['script'],'Quantification/anno_fpkm.py'),
					'--fpkm %s' % os.path.join(self.expRoot,'Summary','transcripts.FPKM.xls'),
					'--genedes %s' % os.path.join(self.expRoot,'trans.description'),
					'--outfile %s' % os.path.join(self.expRoot,'Summary','anno_transcripts.FPKM.xls'),
					'--tr'])
				step_lis.append(step2_add)
		if seqtype=='lncRNA':
			step1 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['script'],'filter_FPKM.py'),
				'%s' % os.path.join(self.expRoot,'lncRNA_transcript.id'),
				'%s' % os.path.join(self.expRoot,'Summary','transcripts.FPKM.xls'),
				'%s' % os.path.join(self.expRoot,'Summary','lncRNA.transcripts.FPKM.xls')])
			step_lis.append(step1)
			step2 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['script'],'filter_FPKM.py'),
				'%s' % os.path.join(self.expRoot,'mRNA_transcript.id'),
				'%s' % os.path.join(self.expRoot,'Summary','transcripts.FPKM.xls'),
				'%s' % os.path.join(self.expRoot,'Summary','mRNA.transcripts.FPKM.xls')])
			step_lis.append(step2)
			step3 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['script'],'filter_FPKM.py'),
				'%s' % os.path.join(self.expRoot,'Unclassified_transcript.id'),
				'%s' % os.path.join(self.expRoot,'Summary','transcripts.FPKM.xls'),
				'%s' % os.path.join(self.expRoot,'Summary','Unclassified.transcripts.FPKM.xls')])
			step_lis.append(step3)
			step4 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['script'],'filter_FPKM.py'),
				'%s' % os.path.join(self.expRoot,'lncRNA_gene.id'),
				'%s' % os.path.join(self.expRoot,'Summary','genes.FPKM.xls'),
				'%s' % os.path.join(self.expRoot,'Summary','lncRNA.genes.FPKM.xls')])
			step_lis.append(step4)
			step5 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['script'],'filter_FPKM.py'),
				'%s' % os.path.join(self.expRoot,'mRNA_gene.id'),
				'%s' % os.path.join(self.expRoot,'Summary','genes.FPKM.xls'),
				'%s' % os.path.join(self.expRoot,'Summary','mRNA.genes.FPKM.xls')])
			step_lis.append(step5)
			step6 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['script'],'filter_FPKM.py'),
				'%s' % os.path.join(self.expRoot,'Unclassified_gene.id'),
				'%s' % os.path.join(self.expRoot,'Summary','genes.FPKM.xls'),
				'%s' % os.path.join(self.expRoot,'Summary','Unclassified.genes.FPKM.xls')])
			step_lis.append(step6)
		strpP_prepare = 'perl %s -i %s -o %s' % (os.path.join(self.softwares['script'],'Quantification/filter_FPKM_one_column_gt1.pl'),os.path.join(self.expRoot,'Summary','genes.FPKM.xls'),os.path.join(self.expRoot,'Summary','filter.genes.FPKM.xls'))
		step_lis.append(strpP_prepare)
		stepP = '%s %s %s %s \n' % (self.softwares['R4.0.3'],os.path.join(self.softwares['script'],'Quantification/plot_den_box.R'),os.path.join(self.expRoot,'Summary','filter.genes.FPKM.xls'),os.path.join(self.expRoot,'Summary'))
		step_lis.append(stepP)
		return ' && \\\n'.join(step_lis),order
	
	def wgcna(self):
		order = 'order wgcna after quantification_summary'
		step = '\\\n\t'.join(['/PUBLIC/software/CANCER/Software/R/R3.4.0/bin/Rscript %s' % os.path.join(self.softwares['script'],'WGCNA_coExpNet.R'),
			'--rpkm %s' % os.path.join(self.expRoot,'Summary','genes.FPKM.xls'),
			'--outdir %s' % os.path.join(self.expRoot,'WGCNA')])
		return step,order
