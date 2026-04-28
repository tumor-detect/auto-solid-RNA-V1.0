
/public/software/anaconda3/envs/python27/bin/python /public/pipeline/RNA/solid/pipeline/RNAseq_cancer_pipeline.py  \
	--infile samplelist \
	--analydir `pwd` \
	--analy-array 1,2.1,3.2,6.3,7.2,8 \
	--newjob fusion.job  \
	--NP all.NP.txt  \
	--compare-group compare.txt


 
