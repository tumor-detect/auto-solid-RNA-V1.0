import os,sys

class QC:

	def __init__(self,sampleID,rawdir,qcDir,fastqs,softwares,databases):
		self.sampleID = sampleID
		self.rawdir = rawdir
		self.qcDir = qcDir
		self.fastqs = fastqs
		self.softwares = softwares
		self.databases = databases
		self.rootdir = os.path.dirname(os.path.dirname(qcDir))

	def rm_rRNA(self,rRNA_rate=10,species='hsa'):
		rawfqs1 = [self.fastqs[alib][alane][0] for alib in self.fastqs for alane in self.fastqs[alib]]
		rawfqs2 = [self.fastqs[alib][alane][1] for alib in self.fastqs for alane in self.fastqs[alib]]
		rawfq1 = rawfqs1[0].replace("_L1_1.fq.gz", "_L*_1.fq.gz")
		rawfq2 = rawfqs2[0].replace("_L1_2.fq.gz", "_L*_2.fq.gz")
		adapts1 = []
		'''for each in rawfqs1:
			if os.path.exists(each[:-5]+'adapter.list.gz'):
				adapts1.append(each[:-5]+'adapter.list.gz')
			else:
				print '\033[1;37;41m Warning: ',rawfqs1,'or it\'s adapter sequence may not exists\033[0m'
				adapts1.append(' ')

		#adapts2 = [each[:-5]+'adapter.list.gz' for each in rawfqs2 if os.path.exists(each[:-5]+'adapter.list.gz')]
		adapts2 = []
		for each in rawfqs2:
			if os.path.exists(each[:-5]+'adapter.list.gz'):
				adapts2.append(each[:-5]+'adapter.list.gz')
			else:
				print '\033[1;37;41m Warning: ',rawfqs2,'or it\'s adapter sequence may not exists\033[0m'
				adapts2.append(' ')

		assert len(adapts1) == len(adapts2)'''
		ln_cmd = []
		if len(rawfqs1) > 1:
			ln_cmd.append('zcat %s |gzip - > %s_1.fq.gz' % (' '.join(rawfqs1),self.sampleID))
			ln_cmd.append('zcat %s |gzip - > %s_2.fq.gz' % (' '.join(rawfqs2),self.sampleID))
		else:
			ln_cmd.append('ln -sf %s %s_1.fq.gz'%(rawfq1,self.sampleID))
			ln_cmd.append('ln -sf %s %s_2.fq.gz'%(rawfq2,self.sampleID))
		'''if len(adapts1) > 1:
			ln_cmd.append('zcat %s |gzip - > %s_1.adapter.list.gz' % (' '.join(adapts1),self.sampleID))
			ln_cmd.append('zcat %s |gzip - > %s_2.adapter.list.gz' % (' '.join(adapts2),self.sampleID))
		else:
			ln_cmd.append('ln -sf %s %s_1.adapter.list.gz'%(adapts1[0],self.sampleID))
			ln_cmd.append('ln -sf %s %s_2.adapter.list.gz'%(adapts2[0],self.sampleID))'''
		if species=='hsa':
			rRNA_lib=self.databases['rRNA_hsa']
		else:
			rRNA_lib=self.databases['rRNA_E']
			
		cmd = '''cd %s

if [[ -a %s_1.fq.gz ]];then
rm %s_1.fq.gz && \\
rm %s_2.fq.gz
fi

%s

gzip -dc %s_1.fq.gz |head -n 40000 > %s_1.test.fq && \\
gzip -dc %s_2.fq.gz |head -n 40000 > %s_2.test.fq && \\
%s \\
  -x %s \\
  -1 %s_1.test.fq -2 %s_2.test.fq --end-to-end --sensitive -p 8 --phred33 --no-mixed -X 600 \\
  -S %s.test.rRNA.sam 2>%s.test.rRNA.stat.txt && \\
rm %s_1.test.fq %s_2.test.fq %s.test.rRNA.sam && \\
maprate=$(tail -n 1 %s.test.rRNA.stat.txt | awk '{print $1}' | awk -F '%%' '{print $1}') && \\
maprate=$(echo $maprate/1|bc) && \\


if [[ $maprate -ge %d ]]; then
  %s \\
    -x %s \\
    -1 %s_1.fq.gz -2 %s_2.fq.gz --end-to-end --sensitive -p 8 --phred33 --no-mixed -X 600 \\
    --un-conc-gz %s.unmap.gz \\
    -S %s.rRNA.sam 2>%s.rRNA.stat.txt && \\
	rm %s.rRNA.sam && \\
    mv %s.unmap.1.gz %s_1.fq.gz && \\
    mv %s.unmap.2.gz %s_2.fq.gz 
 
fi
'''%(os.path.join(self.rawdir,self.sampleID),self.sampleID,self.sampleID,self.sampleID,'  && \\\n'.join(ln_cmd)+'  && \\\n',self.sampleID,self.sampleID,self.sampleID,self.sampleID,self.softwares['bowtie2'],rRNA_lib,self.sampleID,self.sampleID,self.sampleID,self.sampleID,self.sampleID,self.sampleID,self.sampleID,self.sampleID,rRNA_rate,self.softwares['bowtie2'],rRNA_lib,self.sampleID,self.sampleID,self.sampleID,self.sampleID,self.sampleID,self.sampleID,self.sampleID,self.sampleID,self.sampleID,self.sampleID)
		return cmd

	def md5_raw(self):
		order = 'order md5_raw_%s after rm_rRNA_%s' % (self.sampleID,self.sampleID)
		step1 = 'cd %s' % (os.path.join(self.rawdir,self.sampleID))
		step2 = 'md5sum %s_1.fq.gz > %s_1.fq.gz.MD5.txt' % (self.sampleID,self.sampleID)
		step3 = 'md5sum %s_2.fq.gz > %s_2.fq.gz.MD5.txt' % (self.sampleID,self.sampleID)
		return ' && \\\n'.join([step1,step2,step3])+'\n',order

#	def qc(self,opts='',fqcheck=False,rmadapter=False,singlecell=True,cutadapter=False):
        def qc(self,opts='',fqcheck=False,singlecell=False,index=False):

		order = 'order qc_%s after rm_rRNA_%s' % (self.sampleID,self.sampleID)
		rawfq1 = os.path.join(self.rawdir,self.sampleID,'%s_1.fq.gz'%self.sampleID)
		rawfq2 = os.path.join(self.rawdir,self.sampleID,'%s_2.fq.gz'%self.sampleID)
		adapt1 = os.path.join(self.rawdir,self.sampleID,'%s_1.adapter.list.gz'%self.sampleID)
		adapt2 = os.path.join(self.rawdir,self.sampleID,'%s_2.adapter.list.gz'%self.sampleID)
		adapt1 = rawfq1[:-5]+'adapter.list.gz'
		adapt2 = rawfq2[:-5]+'adapter.list.gz'
		fq1 = rawfq1
		fq2 = rawfq2

		cmds = []
		step1 = 'cd %s' % self.qcDir
		cmds.append(step1)
		
		## single cell
		if singlecell == "vazyme":
			step2 = ' &&\\\n'.join(['perl %s -fq %s:%s -index %s -n %s' % (self.softwares['run.generate_adapter'],fq1,fq2,index,self.sampleID),
				'sh generate_adapter.sh'])
			cmds.append(step2)
			fq1 = os.path.join(self.qcDir,'%s_1.fq.gz'%(self.sampleID))
			fq2 = os.path.join(self.qcDir,'%s_2.fq.gz'%(self.sampleID))
			rawfqs1 = [self.fastqs[alib][alane][0] for alib in self.fastqs for alane in self.fastqs[alib]]
			raw_adapts1 = [each[:-5]+'adapter.list.gz' for each in rawfqs1 if os.path.exists(each[:-5]+'adapter.list.gz')]
			step5 = '/data2/Pipeline/RNA/Pipeline/V3/bin/ng_qc -t 4 -i %s,%s -a %s,%s -q 33 -o ./ -L 25 -p 0.4 -N 0.002' % (rawfq1,rawfq2,fq1[:-5]+'adapter.list.gz',fq2[:-5]+'adapter.list.gz')
			cmds.append(step5)
			step5_2 = 'perl %s %s.stat %s' % (os.path.join(self.softwares['bin'],'RawReadsPid3d.pl'),self.sampleID,self.sampleID)
			cmds.append(step5_2)
			step6 = 'echo qc_%s done!' % (self.sampleID)
			cmds.append(step6)

		elif singlecell == "NEB":
			step2 = '\\\n\t'.join(['cutadapt -n 2 -m 60 -O 10 ',
				'-b AAGCAGTGGTATCAACGCAGAGTAC -B AAGCAGTGGTATCAACGCAGAGTAC',
				'-o %s' % os.path.join(self.qcDir,'%s_1.fq.gz'%(self.sampleID)),
				'-p %s' % os.path.join(self.qcDir,'%s_2.fq.gz'%(self.sampleID)),
				'%s %s' % (fq1,fq2)])
			cmds.append(step2)
			fq1 = os.path.join(self.qcDir,'%s_1.fq.gz'%(self.sampleID))
			fq2 = os.path.join(self.qcDir,'%s_2.fq.gz'%(self.sampleID))
			##fqcheck##
			step3 = '\\\n\t'.join(['%s' % os.path.join(self.softwares['fqcheck_adapter']),
				'-a %s' % os.path.join(self.softwares['p7_adapter.fa']),
				'-r %s_1.fq.gz' % (self.sampleID),
				'-l %s_1.adapter.list.gz' % (self.sampleID),
				'-s %s_1.adapter.stat' % (self.sampleID),
				'-c %s_1.adapter.fqcheck' % (self.sampleID)])
			step4 = '\\\n\t'.join(['%s' % os.path.join(self.softwares['fqcheck_adapter']),
				'-a %s' % os.path.join(self.softwares['p5_adapter.fa']),
				'-r %s_2.fq.gz' % (self.sampleID),
				'-l %s_2.adapter.list.gz' % (self.sampleID),
				'-s %s_2.adapter.stat' % (self.sampleID),
				'-c %s_2.adapter.fqcheck' % (self.sampleID)])
			cmds.append(step3)
			cmds.append(step4)
			rawfqs1 = [self.fastqs[alib][alane][0] for alib in self.fastqs for alane in self.fastqs[alib]]
			raw_adapts1 = [each[:-5]+'adapter.list.gz' for each in rawfqs1 if os.path.exists(each[:-5]+'adapter.list.gz')]
			step5 = '/PUBLIC/software/HUMAN/bin/ng_qc -t 4 -i %s,%s -a %s,%s -q 33 -o ./ -L 25 -p 0.4 -N 0.002' % (fq1,fq2,fq1[:-5]+'adapter.list.gz',fq2[:-5]+'adapter.list.gz')
			cmds.append(step5)
			step5_2 = 'perl %s %s.stat %s' % (os.path.join(self.softwares['bin'],'RawReadsPid3d.pl'),self.sampleID,self.sampleID)
			cmds.append(step5_2)
			step6 = 'echo qc_%s done!' % (self.sampleID)
			cmds.append(step6)
		else:
			## qc and adapter
			rawfqs1 = [self.fastqs[alib][alane][0] for alib in self.fastqs for alane in self.fastqs[alib]]
			raw_adapts1 = [each[:-5]+'adapter.list.gz' for each in rawfqs1 if os.path.exists(each[:-5]+'adapter.list.gz')]
			step5 = '\\\n\t'.join(['%s'%self.softwares['fastp'],
				'-i %s'%fq1,
				'-o %s/%s_1.clean.fq.gz'%(self.qcDir,self.sampleID),
				'-I %s'%fq2,
				'-O %s/%s_2.clean.fq.gz'%(self.qcDir,self.sampleID),
				'-j %s/%s.json'%(self.qcDir,self.sampleID),
				'-h %s/%s.html'%(self.qcDir,self.sampleID)])
			cmds.append(step5)
			step6 = 'echo qc_%s done!' % (self.sampleID)
			cmds.append(step6)
		return ' && \\\n'.join(cmds)+'\n',order

	def qc_summary(self,libaryName,analydir):
		order = 'order qc_summary_%s after qc_%s' % (self.sampleID,self.sampleID)
		step1 = 'python %s --infile %s/%s.json --outfile %s/%s.stat --sample %s'%(self.softwares['qc_stat'],self.qcDir,self.sampleID,self.qcDir,self.sampleID,self.sampleID)
		return step1,order

	def md5_clean(self):
		order = 'order md5_clean_%s after qc_%s' % (self.sampleID,self.sampleID)
		step1 = 'cd %s' % self.qcDir
		step3 = 'md5sum %s_1.clean.fq.gz > %s_1.clean.fq.gz.MD5.txt' % (self.sampleID,self.sampleID)
		step5 = 'md5sum %s_2.clean.fq.gz > %s_2.clean.fq.gz.MD5.txt' % (self.sampleID,self.sampleID)
		return ' && \\\n'.join([step1,step3,step5])+'\n',order
		
