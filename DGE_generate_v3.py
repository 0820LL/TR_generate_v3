import os
import sys
import os.path
import re
import argparse
import linecache
root_dir=os.getcwd()

##########################################################################
##########################################################################

#parse the arguments
parser = argparse.ArgumentParser(description="DGE pipline v3.0")
parser.add_argument('--project',help='project name, maybe same with the name of root dir, which will be displayed in the final report title, [REQUIRED]',required=True)
parser.add_argument('--sample',help="sample names(sample1,sample2,...)  warning: order sensitive!",required=True)
parser.add_argument('--fq',help="the original directory of the raw fastq reads, [REQUIRED]",default=None)
parser.add_argument('--mapfile',help="mapfile files (mapfile1,mapfile2.The first column of mapfile is the NH number ,the second column of mapfile is the sample name)",default=None)
parser.add_argument('--raw_dir',help="the original directory of the raw fastq reads keep order in line with mapfile files (raw_dir1,raw_dir2)",default=None)
parser.add_argument('--ad',help="the original directory of the adapter list, ",default=None)
parser.add_argument('--generate_adapter',help='whether to generate the adpepter list y or n',choices=['y','n'],default='n')
parser.add_argument('--index',help='the index number of the fastq files,if you want to generate the adapter list,the parameter will be necessary,diffrent samples split by ,warning: order sensitive!',default=None)
parser.add_argument('--ss',help="strand specific, [REQUIRED]",choices=['no','yes','reverse'],default='no')
parser.add_argument('--number',help="chromosome number, [REQUIRED for density plot]",required=True)
parser.add_argument('--length',help="the length of sequenced reads",required=True)
parser.add_argument('--fa',help="the reference FASTA file, [REQUIRED for mapping]",required=True)
parser.add_argument('--gtf',help="the annotation GTF file, [REQUIRED for mapping]",required=True)
parser.add_argument('--group',help="sample classification, e.g. sample1:sample2,sample3, [REQUIRED]",default=None)
parser.add_argument('--groupname',help="group names, [group1,group2,... ]",default=None)
parser.add_argument('--compare',help="group comparison strategy 1:2,1:3,..., [REQUIRED]",default=None)
parser.add_argument('--venn',help="venn and cluster mode, defult is all groups, 1:2_1:3,1:3_2:3,... ",default=None)
parser.add_argument('--goann',help="gene to GO annotations file, [REQUIRED]",default=None)
parser.add_argument('--species',help="abbreviation for species, for kegg pathway. (note: all species: /PUBLIC/database/RNA/kobas2.0-data-20120208/seq_pepi_v2), [defalut=ko]",default='kaas')
parser.add_argument('--ppi_number',help="species code, (ref file: /PUBLIC/database/RNA/string_ppi/species.v9.0.txt)",default=None)
parser.add_argument('--ppi_blast',help="whether to run blast to get the protein-protein interactions",choices=['y','n'],default=None)
parser.add_argument('--genenamefile',help="genenamefile, 1st col is geneID, 2nd col is genename",default=None)
parser.add_argument('--ex',help="the steps you do not wanna perform",default=None)

def checkSample(sample):
  PatChar = r'^(CON|PRN|AUX|CLOCK\$|NUL|COM[1-9]|LPT1|[\W]+)$'
  PatNum = r'^\d+'
  MAXLEN = 8
  message = '''can only use word/digit/underscore in sample, start with word. 
Max length is 8, no windows reserved words like CON PRN...'''
  if re.search(PatChar,sample) or re.search(PatNum,sample) or len(sample) > MAXLEN:
    print '%s<==invalid name' % s
    print message
    exit()

# extract, check and display the parameters
argv = vars(parser.parse_args())
project=argv['project'].strip()
display = open(root_dir + '/' + 'DGE_command.txt','w')
display.write('project: %s\n' % (project))
sample=argv['sample'].strip()
samples=sample.split(',')
samples_tmp=list(set(samples))
for s in samples:
  checkSample(s)
assert len(samples)==len(samples_tmp)
display.write('sample: %s\n' % (sample))
mapfiles=[each.strip() for each in argv['mapfile'].strip().split(',') if each.strip() != '']
mapfile=' '.join(mapfiles)
assert not os.system('cat %s |sort -u |sed \'/^$/d\' >%s' % (mapfile,root_dir+'/libraryID'))
if argv['fq']:
	fq = argv['fq'].strip()
	fq = os.path.abspath(fq)
else:
	assert not os.system('mkdir raw_data')
	assert not os.system('perl /PUBLIC/source/RNA/RefRNA/ln_raw_data.pl %s %s se raw_data' %(argv['mapfile'],argv['raw_dir']))
	fq = root_dir + '/raw_data'
for each in samples:
	fq_tmp1=fq+'/'+each+'.fq.gz'
#	fq_tmp2=fq+'/'+each+'_2.fq.gz'
	assert os.path.isfile(fq_tmp1)
#	assert os.path.isfile(fq_tmp2)
display.write('fq: %s\n' % (fq))
display.write('adapter: ')
if argv['ad']:
	ad = argv['ad'].strip()
	ad = os.path.abspath(ad)
	for each in samples:
		ad_tmp1=ad+'/'+each+'.adapter.list.gz'
#		ad_tmp2=ad+'/'+each+'_2.adapter.list.gz'
		assert os.path.isfile(ad_tmp1)
#		assert os.path.isfile(ad_tmp2)
	display.write('%s' % (ad))
display.write('\n')
if argv['generate_adapter']:
	generate_adapter = argv['generate_adapter'].strip()
else:
	generate_adapter = 'n'
display.write('generate_adapter: %s\n' % (generate_adapter))
display.write('index: \n')
if argv['index']:
	index=argv['index'].strip()
	indexes=index.split(',')
	assert len(samples)==len(indexes)
	for i,index_tmp in enumerate(samples):
		display.write('%s:\t%s\n' % (index_tmp,indexes[i]))
display.write('\n')
if generate_adapter == 'y':
	if argv['ad'] !=None:
		print 'Error:  the parameters --ad and --generate_adapter are not consensus!\n'
		exit()
	if argv['index'] == None:
		print 'Error:  the parameters --index and --generate_adapter are not consensus!\n'
		exit()
else:
	if argv['index'] != None:
		print 'Error:  the parameters --index and --generate_adapter are not consensus!\n'
		exit()
#-----------------------------------------------------------------------
all_content=set([1,2,3,4,5,6])
if argv['ex'] != None:
	excludes=argv['ex'].strip().strip(',').strip().split(',')
	excludes=[int(each.strip()) for each in excludes]
	for each1 in excludes:
		assert each1 in all_content
else:
	excludes=[] #list
includes=all_content-set(excludes) #set
#-----------------------------------------------------------------------
ss = argv['ss'].strip()
display.write('ss: %s\n' % (ss))
number = argv['number'].strip()
display.write('chr number: %s\n' % (number))
if argv['length']:
	length = argv['length'].strip()
else:
	length = '100'
display.write('length: %s\n' % (length))
fa = argv['fa'].strip()
fa = os.path.abspath(fa)
suffix_fa=['.1.bt2','.2.bt2','.3.bt2','.4.bt2','.rev.1.bt2','.rev.2.bt2']
for each in suffix_fa:
	fa_tmp=fa+each
	assert os.path.isfile(fa_tmp)
display.write('fa: %s\n' % (fa))
gtf = argv['gtf'].strip()
gtf = os.path.abspath(gtf)
assert os.path.isfile(gtf)
display.write('gtf: %s\n' % (gtf))
#-----------------------------------------------------------------------
if set([2]).issubset(includes):
	if argv['group'] == None:
		group=samples
		group_iter=samples
		group=sample
		flag_repeat=False
	else:
		groups_iter=[]
		if ':' in argv['group'].strip(':'):
			flag_repeat=True
			groups=[each.strip().split(':') for each in argv['group'].strip().strip(':').split(',') if each.strip() != '']
			for each in groups:
				groups_iter+=each
			group_iter_n=[]
			for each in groups_iter:
				if each not in group_iter_n:
					group_iter_n.append(each)
			groups_iter=group_iter_n
			group=','.join([':'.join(each) for each in groups])
		else:
			flag_repeat=False
			groups=[each.strip() for each in argv['group'].strip().split(',') if each.strip() != '']
			for each in groups:
				groups_iter.append(each)
			group=','.join(groups)
		assert set(groups_iter).issubset(samples)
	group_iter=','.join(groups_iter)
	
	if argv['groupname'] == None:
		if flag_repeat == False:
			groupname_tmp=groups
		else:
			groupname_tmp=['group'+str(k+1) for k in range(len(groups))]
	else:
		groupname_tmp=[each.strip() for each in argv['groupname'].split(',') if each.strip() != '']
		assert len(groupname_tmp) == len(groups)
	groupname=','.join(groupname_tmp)

	compare = argv['compare'].strip()
	compare_samples=[each.strip().split(':') for each in argv['compare'].strip().split(',') if each.strip() != '']
	M=[]
	for each1 in compare_samples:
		temp1=[]
		for each2 in each1:
			assert each2.isdigit()
			M.append(int(each2))
	assert max(M) <= len(groupname_tmp)
	assert min(M) > 0
	compare_sample=','.join([','.join(each) for each in compare_samples])
	temp2=[]
	for each1 in compare_samples:
		temp1=[]
		for each2 in each1:
			temp1.append(groupname_tmp[int(each2)-1])
		temp2.append(':'.join(temp1))
	compares_name=','.join(temp2)
	if argv['venn']:
		venn = argv['venn'].strip()
		com_pairs=compare.split(',')
		venn_clusters=[each.split('_') for each in venn.split(',')]
		temp1=[]
		for each1 in venn_clusters:
			temp2=[]
			for each2 in each1:
				assert each2 in com_pairs
				temp3=each2.split(':')
				assert len(temp3) == 2
				temp2.append(groupname_tmp[int(temp3[0])-1]+':'+groupname_tmp[int(temp3[1])-1])
			temp1.append('_'.join(temp2))
		venn_cluster_name=','.join(temp1)
	else:
		venn = compare.replace(',','_')
		venn_cluster_name=compares_name.replace(',','_')
	groups=group.split(',')
	groupnames=groupname.split(',')
	compares=compare.split(',')
	venns=venn.split(',')
	for each in groups:
		temp=each.split(':')
		assert len(temp)==len(list(set(temp)))
		for each2 in temp:
			assert each2 in samples
	assert len(groups)==len(groupnames)
	display.write('\ngroups:\n')
	assert len(groupnames)==len(list(set(groupnames)))
	for i,each in enumerate(groupnames):
		group_tmp=groups[i]
		group_tmp=group_tmp.replace(':',',')
		display.write('%s: %s\n' % (each,group_tmp))
	display.write('\ncompare:\n')
	compare_name=[]
	for i,each in enumerate(compares):
		display.write('%s: \t' % (each))
		compare_tmp=each
		compare_tmp=compare_tmp.split(':')
		assert len(compare_tmp)==2
		compare_name=groupnames[int(compare_tmp[0])-1]+'vs'+groupnames[int(compare_tmp[1])-1]
		display.write('%s\n' % (compare_name))
	display.write('\nvenn:\n')
	for i,each in enumerate(venns):
		display.write('%s: \t' % (each))
		venn_tmp=each
		venn_tmp=venn_tmp.split('_')
		for each2 in venn_tmp:
			assert each2 in compares
			venn_tmp2=each2.split(':')
			venn_name=groupnames[int(venn_tmp2[0])-1]+'vs'+groupnames[int(venn_tmp2[1])-1]
			display.write('%s,' % (venn_name))
		display.write('\n')
	display.write('\n')
#-----------------------------------------------------------------------
if set([1,2,4]).issubset(includes):
	goann = argv['goann'].strip()
	goann = os.path.abspath(goann)
	assert os.path.isfile(goann)
display.write('goann: %s\n' % (goann))
if set([1,2,5]).issubset(includes):
	if argv['species']:
		species = argv['species'].strip()
	else:
		species = 'kaas'
	display.write('KEGG species: %s\n' % (species))
if set([1,2,6]).issubset(includes):
	display.write('PPI number: ')
	if argv['ppi_number']:
		ppi_number = argv['ppi_number'].strip()
		display.write('%s' % (ppi_number))
	display.write('\n')
	display.write('PPI blast: ')
	if argv['ppi_blast']:
		ppi_blast = argv['ppi_blast'].strip()
		display.write('%s' % (ppi_blast))
	display.write('\n')
	if argv['ppi_blast']:
		if argv['ppi_number'] == None:
			print 'Error:  the parameters --ppi_blast and --ppi_number are not consensus!\n'
			exit()
	else:
		if argv['ppi_number']:
			print 'Error:  the parameters --ppi_blast and --ppi_number are not consensus!\n'
			exit()
display.write('genenamefile: ')
if argv['genenamefile']:
	genenamefile = argv['genenamefile'].strip()
	genenamefile = os.path.abspath(genenamefile)
	assert os.path.isfile(genenamefile)
	display.write('%s' % (genenamefile))
display.write('\n')
display.close()

#################################################################################################
#for QC
code='''
import os

fq='%s'
ss='%s'
fa='%s'
gtf='%s'
sample='%s'
root_dir='%s'
''' % (fq,ss,fa,gtf,sample,root_dir)

if argv['ad']:
	code+='''
ad='%s'
''' % (ad)

if generate_adapter != 'n':
	code+='''
index='%s'
''' % (index)

if argv['mapfile']:
	code+='''
mapfile='%s'
''' % (root_dir+'/libraryID')

code+='''
assert not os.system('mkdir QC_DGE')
os.chdir('QC_DGE')
f=open(root_dir+'/QC_DGE/generate_QC.sh','w')
f.write('perl /PUBLIC/source/RNA/QC/QC_v2/allrunQC_v2.1.pl -fq %s -se-pe se -n %s -o %s -spe %s -R %s -G %s -bed %s ' %(fq,sample,root_dir+'/QC_DGE',ss,fa,gtf,root_dir+'/QC_DGE/sorted.bed'))
'''
if argv['mapfile']:
	code+='''
f.write('-mapfile %s ' % (mapfile))
'''

if argv['ad']:
	code+='''
f.write('-ad %s\\n' % (ad))
'''

if argv['index']:
	code+='''
f.write('-m_ad y -index %s\\n' % (index))
'''

code+='''
f.write('\\n')
f.close()
assert not os.system('sh %s/generate_QC.sh' %(root_dir+'/QC_DGE'))

samples=sample.split(',')
for eachsample in samples:
	os.chdir(root_dir+'/QC_DGE/'+eachsample)
	assert not os.system('qsub -V -l vf=5g,p=8 -cwd %s' % (eachsample+'_QC_se.sh'))
os.chdir(root_dir)
'''
open('DGE_step1_QC.py','w').write(code)
#################################################################################################
#for QCreport
code='''
import os

sample='%s'
project='%s'
root_dir='%s'
''' % (sample,project,root_dir)

code+='''
assert not os.system('mkdir QC_DGE/QCreport')
os.system('sh /PUBLIC/source/RNA/QC/QC_v2/QCreport/DGE_QCreport.sh -dir %s -sample %s -title %s -results %s' %(root_dir+'/QC_DGE',sample,project,root_dir+'/QC_DGE/QCreport'))
'''
open('DGE_step1_QCreport.py','w').write(code)
#################################################################################################
#for Diff analysis and Curves
code=''
code='''
import re
import os

sample='%s'
ss='%s'
gtf='%s'
group='%s'
groupname='%s'
compare='%s'
venn='%s'
fa='%s'
n='%s'
length='%s'
root_dir='%s'
''' % (sample,ss,gtf,group,groupname,compare,venn,fa,number,length,root_dir)

if argv['genenamefile']:
	code+='''
genenamefile='%s'
''' % (genenamefile)

if set([1,2]).issubset(includes):
	code+='''
##Diff
sam=root_dir+'/QC_DGE/sam'
samples = sample.split(',')
readcount=[]
for eachsample in samples:
	temp='%s/Diff_DGE/readcount/%s.readcount' % (root_dir,eachsample)
	readcount.append(temp)
readcount=','.join(readcount)
assert not os.system('mkdir Diff_DGE')
os.chdir('Diff_DGE')
f=open(root_dir+'/Diff_DGE/generate_Diff.sh','w')
f.write('perl /PUBLIC/source/RNA/RefRNA/DGE/DiffExpression/runDiff_analysis_v2.2.pl -fa %s -sam %s -g %s -o %s -group %s -groupname %s -compare %s -venn %s -spe %s  -i %s ' % (fa,sam,gtf,root_dir+'/Diff_DGE',group,groupname,compare,venn,ss,readcount))
'''
	if argv['genenamefile']:
		code+='''
f.write(' -genenamefile %s ' % (genenamefile))
'''
	code+='''
f.write('\\n')
f.close()
assert not os.system('sh generate_Diff.sh')
assert not os.system('qsub -V -cwd -l vf=5G -l p=1 runDiff_analysis.sh')
os.chdir(root_dir)
'''

if set([1,2,3]).issubset(includes):
	code+='''
##Curve
bam=root_dir+'/QC_DGE/bam'
assert not os.system('mkdir Curve_DGE')
os.chdir('Curve_DGE')
f=open(root_dir+'/Curve_DGE/generate_Curve.sh','w')
f.write('perl /PUBLIC/source/RNA/RefRNA/DGE/Curve/runCurve_v2.pl -bam %s -sample %s -n %s -r %s -fa %s -gtf %s -o %s\\n' % (bam,sample,n,length,fa,gtf,root_dir+'/Curve_DGE'))
f.close()
assert not os.system('sh generate_Curve.sh')

assert not os.system('qsub -V -cwd -l vf=10G -l p=1 runSaturation.sh')
assert not os.system('qsub -V -cwd -l vf=5G -l p=1 rundensity.sh')

os.chdir(root_dir)
'''
open('DGE_step2_Diff_Curves.py','w').write(code)
#################################################################################################
#for Enrichment and PPI

code='''
import os

fa='%s'
goann='%s'
groupname='%s'
compare='%s'
root_dir='%s'
''' % (fa,goann,groupname,compare,root_dir)

code+='''
groupnames=groupname.split(',')
compares=compare.split(',')
'''

if set([1,2,4]).issubset(includes):
	code+='''
##GOSeq
assert not os.system('mkdir GOSeq_DGE')
os.chdir('GOSeq_DGE')
go=open(root_dir+'/GOSeq_DGE/runGOSeq.sh','w')
for each in compares:
	temp=each.split(':')
	dir=groupnames[int(temp[0])-1]+'vs'+groupnames[int(temp[1])-1]
	id=root_dir+'/Diff_DGE/Diff/'+dir+'/'+dir+'.diffgeneID'
	out=root_dir+'/GOSeq_DGE/'+dir
	length=root_dir+'/Diff_DGE/Diff/genelength'
	result=root_dir+'/GOSeq_DGE/'+dir+'/'+dir+'.GO_enrichment_result.xls'
	diff=root_dir+'/Diff_DGE/Diff/'+dir+'/'+dir+'.diffgene.xls'
	go.write('###############%s#####################\\n' % (dir))
	go.write('mkdir %s\\n' % (root_dir+'/GOSeq_DGE/'+dir))
	go.write('perl /PUBLIC/source/RNA/RefRNA/DGE/Enrichment/goseq_graph_v3.pl -i %s -goann %s -n %s -o %s -length %s\\n' % (id,goann,dir,out,length))
	go.write('perl /PUBLIC/source/RNA/RefRNA/DGE/Enrichment/changeGO_up_down.pl %s %s %s\\n' % (result,diff,root_dir+'/GOSeq_DGE/'+dir+'/'+dir+'.GO_enrichment_result_up_down.xls'))
	go.write('/PUBLIC/software/public/System/R-2.15.3/bin/Rscript /PUBLIC/source/RNA/RefRNA/DGE/Enrichment/goBar.R %s %s %s\\n' % (result,out,dir))
	go.write('/PUBLIC/software/public/System/R-2.15.3/bin/Rscript /PUBLIC/source/RNA/RefRNA/DGE/Enrichment/goBar2.R %s %s %s\\n' % (root_dir+'/GOSeq_DGE/'+dir+'/'+dir+'.GO_enrichment_result_up_down.xls',out,dir))
go.close()
assert not os.system('qsub -V -cwd -l vf=5G -l p=1 runGOSeq.sh')
os.chdir(root_dir)
'''

if argv['genenamefile'] == None:
	code+='''
assert not os.system('mkdir Blast_DGE')
os.chdir('Blast_DGE')
f=open(root_dir+'/Blast_DGE/runBlast_swissprot.sh','w')
query=root_dir+'/Diff_DGE/Diff/diffgene_union.seq'
outdir1=root_dir+'/Blast_DGE/Blast_Swissprot/'
out=root_dir+'/Blast_DGE/Blast_Swissprot/diffgene_union.seq.blastout'
f.write('echo start blastx\\ndate\\n')
f.write('mkdir %s\\n' % (outdir1))
f.write('/PUBLIC/software/public/Alignment/ncbi-blast-2.2.28+/bin/blastx -query %s -db /PUBLIC/database/Common/SwissProt/uniprot_sprot.fasta -evalue 1e-5 -outfmt 5 -max_target_seqs 1 -num_threads 10 -out %s\\n' % (query,out))
f.write('echo blastx end\\ndate\\n')
f.write('perl /PUBLIC/source/RNA/RefRNA/DGE/scriptdir/extractIDsEVxml.pl %s %s\\n' % (out,outdir1+'/diffgene_union.genenames'))
compare=root_dir+'/Diff_DGE/Diff/compare.txt'
indir=root_dir+'/Diff_DGE/Diff/'
outdir2=root_dir+'/Diff_DGE/Diff/DiffGeneList'
f.write('mkdir %s\\n' % (outdir2))
f.write('perl /PUBLIC/source/RNA/RefRNA/DGE/scriptdir/getdiffGN.up_down_v2.pl %s %s %s %s\\n' % (indir,compare,outdir1+'/diffgene_union.genenames',outdir2))
f.close()
assert not os.system('qsub -V -cwd -l vf=5G -l p=10 runBlast_swissprot.sh')
os.chdir(root_dir)
'''

if set([1,2,5]).issubset(includes):
	code+='''
species='%s'
''' % (species)
	if argv['species'] == 'kaas':
		code+='''
##KOBAS
assert not os.system('mkdir KOBAS_DGE')
os.chdir('KOBAS_DGE')
ko=open(root_dir+'/KOBAS_DGE/runKOBAS.sh','w')
path=open(root_dir+'/KOBAS_DGE/runpathway.sh','w')
query=root_dir+'/Diff_DGE/Diff/diffgene_union.seq'
blastout=root_dir+'/Blast_DGE/KOBAS_blast.xml'
ko.write('###############Run kaas#####################\\n')
ko.write('/PUBLIC/software/public/Annotation/kaas_sa2/bin/auto_annotate.pl -n -s %s\\n' % (query))
ko.write('python /PUBLIC/software/RNA/gene_annotation/scripts/convert2kobas_v1.py %s /PUBLIC/database/Common/KEGG/kos %s\\n' % (root_dir + '/Diff_DGE/Diff/diffgene_union.seq.ko',root_dir + '/KOBAS_DGE/koID.annotation'))
for each in compares:
	temp=each.split(':')
	dir=groupnames[int(temp[0])-1]+'vs'+groupnames[int(temp[1])-1]
	id=root_dir+'/Diff_DGE/Diff/'+dir+'/'+dir+'.diffgeneID'
	out=root_dir+'/KOBAS_DGE/'+dir
	ko.write('###############%s#####################\\n' % (dir))
	path.write('echo "###############%s#####################"\\n' % (dir))
	path.write('cd %s\\n' % (out))
	result=root_dir+'/KOBAS_DGE/'+dir+'/'+'add.'+dir+'.identify.xls'
	diff=root_dir+'/Diff_DGE/Diff/'+dir+'/'+dir+'.diffgene.xls'
	diffID=root_dir+'/Diff_DGE/Diff/'+dir+'/'+dir+'.diffgeneID'
	path.write('python /PUBLIC/source/RNA/RefRNA/DGE/KEGG_Enrichment/bin/pathway_annotation_flow_parallel_annotationfault_tolerant3.pyc --table %s --diff %s\\n' % (result,diff))
	path.write('mv %s %s\\n' % ('add.'+dir+'.identify.xls_rendered_html_detail.html',dir+'.html'))
	ko.write('mkdir %s\\n' % (root_dir+'/KOBAS_DGE/'+dir))
	ko.write('cd %s\\n' % (root_dir+'/KOBAS_DGE/'+dir))
	ko.write('perl /PUBLIC/source/RNA/noRef/KEGG_enrichment/runKEGG_enrich.pl -diff %s -ko %s -g %s\\n' % (diffID,root_dir + '/KOBAS_DGE/koID.annotation',dir))
path.close()
ko.close()
assert not os.system('qsub -V -cwd -l vf=5G -l p=10 runKOBAS.sh')
os.chdir(root_dir)
'''
	else:
		code+='''
##KOBAS
assert not os.system('mkdir KOBAS_DGE')
os.chdir('KOBAS_DGE')
ko=open(root_dir+'/KOBAS_DGE/runKOBAS.sh','w')
path=open(root_dir+'/KOBAS_DGE/runpathway.sh','w')
query=root_dir+'/Diff_DGE/Diff/diffgene_union.seq'
blastout=root_dir+'/Blast_DGE/KOBAS_blast.xml'
ko.write('###############Run Blast#####################\\n')
ko.write('mkdir %s\\n' % (root_dir+'/Blast_DGE/'))
ko.write('perl /PUBLIC/source/RNA/RefRNA/DGE/KEGG_Enrichment/bin/KEGG_step1_blast.pl %s %s %s %s\\n' % (query,species,blastout,root_dir+'/Blast_DGE/KOBAS_blast.sh'))
ko.write('sh %s\\n' % (root_dir+'/Blast_DGE/KOBAS_blast.sh'))
for each in compares:
	temp=each.split(':')
	dir=groupnames[int(temp[0])-1]+'vs'+groupnames[int(temp[1])-1]
	id=root_dir+'/Diff_DGE/Diff/'+dir+'/'+dir+'.diffgeneID'
	out=root_dir+'/KOBAS_DGE/'+dir
	ko.write('###############%s#####################\\n' % (dir))
	path.write('echo "###############%s#####################"\\n' % (dir))
	path.write('cd %s\\n' % (out))
	result=root_dir+'/KOBAS_DGE/'+dir+'/'+'add.'+dir+'.identify.xls'
	diff=root_dir+'/Diff_DGE/Diff/'+dir+'/'+dir+'.diffgene.xls'
	path.write('python /PUBLIC/source/RNA/RefRNA/DGE/KEGG_Enrichment/bin/pathway_annotation_flow_parallel_annotationfault_tolerant3.pyc --table %s --diff %s\\n' % (result,diff))
	path.write('mv %s %s\\n' % ('add.'+dir+'.identify.xls_rendered_html_detail.html',dir+'.html'))
	ko.write('mkdir %s\\n' % (root_dir+'/KOBAS_DGE/'+dir))
	script=root_dir+'/KOBAS_DGE/'+dir+'/run.sh'
	ko.write('perl /PUBLIC/source/RNA/RefRNA/DGE/KEGG_Enrichment/bin/KEGG_step2_enrich.pl -id %s -out-dir %s -species %s -blast-result %s -sample-names %s>%s\\n' % (id,out,species,blastout,dir,script))
	ko.write('sh %s\\n' % (script))
path.close()
ko.close()
assert not os.system('qsub -V -cwd -l vf=5G -l p=10 runKOBAS.sh')
os.chdir(root_dir)
'''

if set([1,2,6]).issubset(includes):
	code+='''
ppi_number='%s'
ppi_blast='%s'
''' % (ppi_number,ppi_blast)
	code+='''
##PPI
assert not os.system('mkdir PPI_DGE')
os.chdir('PPI_DGE')
ppi=open(root_dir+'/PPI_DGE/runPPI.sh','w')
for each in compares:
	temp=each.split(':')
	dir=groupnames[int(temp[0])-1]+'vs'+groupnames[int(temp[1])-1]
	id=root_dir+'/Diff_DGE/Diff/'+dir+'/'+dir+'.diffgeneID'
	seq=root_dir+'/Diff_DGE/Diff/Diff_Gene_Seq/'+dir+'.diffgene.seq'
	out=root_dir+'/PPI_DGE/'+dir
	ppi.write('mkdir %s\\n' % (out))
'''
	if argv['ppi_blast'] == 'y':
		code+='''
	ppi.write('python /PUBLIC/source/RNA/RefRNA/DGE/Enrichment/BLASTX_TO_PPI_v4.py --species %s --fa %s --outdir %s --name %s\\n' % (ppi_number,seq,out,dir))
'''
	else:
		code+='''
	ppi_dir='/PUBLIC/database/RNA/PPI_ref/PPI_99'
	ppi.write('python /PUBLIC/source/RNA/RefRNA/DGE/Enrichment/get_my_PPI.py -p %s -g %s -o %s\\n' % (ppi_dir+'/PPI_'+ppi_number+'.txt',id,out+'/'+dir+'.ppi.txt'))
'''
	code+='''
ppi.close()
assert not os.system('qsub -V -cwd -l vf=5G -l p=10 runPPI.sh')
os.chdir(root_dir)
'''

open('DGE_step3_Enrichment.py','w').write(code)

#################################################################################################
#for Report

code='''
import re
import os
import os.path
import glob
import sys
import gzip
import linecache
reload(sys)
sys.setdefaultencoding('utf-8')
from django.template import Template, Context, loader
from django.conf import settings
settings.configure(DEBUG=True, TEMPLATE_DEBUG=True,TEMPLATE_DIRS=('/PUBLIC/source/RNA/RefRNA/Report/template/DGE',))

fq='%s'
sample='%s'
groupname='%s'
compare='%s'
project='%s'
root_dir='%s'
ss='%s'
venn_cluster_name='%s'
flag_uniform=True
''' % (fq,sample,groupname,compare,project,root_dir,ss,venn_cluster_name)

code+='''
samples=sample.split(',')

##md5
dir=root_dir+'/QC_DGE'
cd_md5=open(root_dir+'/QC_DGE/cd_md5.sh','w')
rd_md5=open(root_dir+'/QC_DGE/rd_md5.sh','w')
os.chdir(dir)
cd_md5.write('cd %s/cleandata\\n' % (dir))
cd_md5.write('echo -e \"md5\\tfile\">cd_md5.txt\\n')
rd_md5.write('echo -e \"md5\\tfile\">rd_md5.txt\\n')
for each in samples:
	cd_md5.write('md5sum %s.clean.fq>>cd_md5.txt\\n' % (each))
	cd_md5.write('unlink %s.clean.fq\\n' %(each))
	cd_temp=root_dir+'/QC_DGE/'+each+'/clean_data/'+each+'.clean.fq'
	cd_md5.write('gzip %s\\n' % (cd_temp))
	cd_md5.write('ln -s %s.gz %s\\n' % (cd_temp,dir+'/cleandata/'))
	rd_md5.write('cd %s\\n' % (root_dir+'/QC_DGE/'+each))
	rd_md5.write('md5sum %s>>%srd_md5.txt\\n' % (each+'.fq.gz',root_dir+'/QC_DGE/'))
cd_md5.close()
rd_md5.close()
#assert not os.system('qsub -V -l vf=3g,p=1 -cwd cd_md5.sh')
#assert not os.system('qsub -V -l vf=3g,p=1 -cwd rd_md5.sh')
os.chdir(root_dir)

os.chdir(root_dir+'/QC_DGE')
rdmd5=glob.glob('rd_md5.txt')
if not rdmd5:
        assert not os.system('qsub -V -l vf=3g,p=1 -cwd rd_md5.sh')
os.chdir(root_dir+'/QC_DGE/cleandata')
cdmd5=glob.glob('cd_md5.txt')
if not cdmd5:
	os.chdir(root_dir+'//QC_DGE')
        assert not os.system('qsub -V -l vf=3g,p=1 -cwd cd_md5.sh')
os.chdir(root_dir)
'''

if set([1,2,5]).issubset(includes):
	code+='''
##path
path=root_dir+'/KOBAS_DGE/runpathway.sh'
assert not os.system('sh %s' % (path))
'''

code+='''
#results

assert not os.system('mkdir '+project+'_DGE_result')
assert not os.system('mkdir '+project+'_DGE_result/'+project+'_results')
assert not os.system('mkdir '+project+'_DGE_result/'+project+'_report')
#assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/Report/src %s' % (root_dir+'/'+project+'_DGE_result/'+project+'_report'))
os.chdir(project+'_DGE_result/'+project+'_results')

result_order=0

result_order+=1
assert not os.system('mkdir '+str(result_order)+'.OriginalData')
rd_order=str(result_order)
os.chdir(''+str(result_order)+'.OriginalData')
#assert not os.system('cp %s\\n' % (root_dir+'/QC_DGE/rd_md5.txt'))
for eachsample in samples:
	f_in=gzip.open(root_dir+'/QC_DGE/'+eachsample+'/'+eachsample+'.fq.gz','rb')
	f_out=[f_in.readline() for num in range(20)]
	open(eachsample+'.example.fq.txt','w').writelines(f_out)
	f_in.close()
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/OriginalData.README.txt %s' % (root_dir+'/'+project+'_DGE_result/'+project+'_results/'+str(result_order)+'.OriginalData/'+str(result_order)+'.1OriginalData.README.txt'))
os.chdir(root_dir+'/'+project+'_DGE_result/'+project+'_results')

result_order+=1
assert not os.system('mkdir '+str(result_order)+'.QC')
qc_order=str(result_order)
os.chdir(''+str(result_order)+'.QC')
assert not os.system('mkdir '+str(result_order)+'.1.ErrorRate')
os.chdir(''+str(result_order)+'.1.ErrorRate')
for eachsample in samples:
	assert not os.system('cp %s %s' % (root_dir+'/QC_DGE/'+eachsample+'/clean_data/'+eachsample+'.Error.png',eachsample+'.error_rate_distribution.png'))
	assert not os.system('cp %s %s' % (root_dir+'/QC_DGE/'+eachsample+'/clean_data/'+eachsample+'.Error.pdf',eachsample+'.error_rate_distribution.pdf'))
	assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/Error.README.txt %s' % (root_dir+'/'+project+'_DGE_result/'+project+'_results/'+str(result_order)+'.QC/'+str(result_order)+'.1.ErrorRate/'+str(result_order)+'.1ErrorRate.README.txt'))
os.chdir('..')

assert not os.system('mkdir '+str(result_order)+'.2.GC')
os.chdir(''+str(result_order)+'.2.GC')
for eachsample in samples:
	assert not os.system('cp %s %s' % (root_dir+'/QC_DGE/'+eachsample+'/clean_data/'+eachsample+'.GC.png',eachsample+'.GC_content_distribution.png'))
	assert not os.system('cp %s %s' % (root_dir+'/QC_DGE/'+eachsample+'/clean_data/'+eachsample+'.GC.pdf',eachsample+'.GC_content_distribution.pdf'))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/GC.README.txt %s' % (root_dir+'/'+project+'_DGE_result/'+project+'_results/'+str(result_order)+'.QC/'+str(result_order)+'.2.GC/'+str(result_order)+'.2GC.README.txt'))
os.chdir('..')

assert not os.system('mkdir '+str(result_order)+'.3.ReadsClassification')
os.chdir(''+str(result_order)+'.3.ReadsClassification')
for eachsample in samples:
	assert not os.system('cp %s %s' % (root_dir+'/QC_DGE/'+eachsample+'/files/'+eachsample+'.pie3d.png',eachsample+'.raw_reads_classification.png'))
	assert not os.system('cp %s %s' % (root_dir+'/QC_DGE/'+eachsample+'/files/'+eachsample+'.pie3d.pdf',eachsample+'.raw_reads_classification.pdf'))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/Filter.README.txt %s' % (root_dir+'/'+project+'_DGE_result/'+project+'_results/'+str(result_order)+'.QC/'+str(result_order)+'.3.ReadsClassification/'+str(result_order)+'.3ReadsClassification.README.txt'))
os.chdir('..')

assert not os.system('mkdir '+str(result_order)+'.4.DataTable')
os.chdir(''+str(result_order)+'.4.DataTable')
#for eachsample in samples:
#	assert not os.system('cp %s %s' % (root_dir+'/QC_DGE/'+eachsample+'/files/dataTable',eachsample+'datatable.xls'))
alldataTable = [ root_dir + '/QC_DGE/' + sample + '/files/dataTable' for sample in samples]
alldataTable = ' '.join(alldataTable)
assert not os.system("cat %s|cut -f 1-8 > %s " % ( alldataTable,'datatable.xls.1'))
assert not os.system('cat %s %s > %s' % ('/PUBLIC/source/RNA/RefRNA/DGE/Report/datatable','datatable.xls.1','datatable.xls'))
assert not os.system('rm %s' % ('datatable.xls.1'))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/DataTable.README.txt %s' % (root_dir+'/'+project+'_DGE_result/'+project+'_results/'+str(result_order)+'.QC/'+str(result_order)+'.4.DataTable/'+str(result_order)+'.4DataTable.README.txt'))
os.chdir(root_dir+'/'+project+'_DGE_result/'+project+'_results')

result_order+=1
assert not os.system('mkdir '+str(result_order)+'.Mapping')
mp_order=str(result_order)
os.chdir(''+str(result_order)+'.Mapping')
assert not os.system('mkdir '+str(result_order)+'.1.MapStat')
os.chdir(''+str(result_order)+'.1.MapStat')
for eachsample in samples:
	assert not os.system('cut -f 2 %s > %s' % (root_dir+'/QC_DGE/'+eachsample+'/files/'+eachsample+'.stat',eachsample+'.stat.2'))
assert not os.system('cut -f 1 %s > %s' % (root_dir+'/QC_DGE/'+samples[0]+'/files/'+samples[0]+'.stat','stat.1'))
allStat_2 = [ root_dir + '/' + project+'_DGE_result/'+project+'_results/'+str(result_order)+'.Mapping/'+str(result_order)+'.1.MapStat/'+sample+'.stat.2'  for sample in samples ]
allStat_2 = ' '.join(allStat_2)
assert not os.system('paste %s %s > %s' % ('stat.1',allStat_2,'MapStat.xls'))
assert not os.system("sed -i -e '/^Reads mapped in proper pairs/d' %s" % ('MapStat.xls'))
assert not os.system("sed -i -e '/^Proper-paired reads map to different chrom/d' %s" % ('MapStat.xls'))
assert not os.system('rm stat.1')
assert not os.system('rm *.stat.2')
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/MapStat.README.txt %s' % (root_dir+'/'+project+'_DGE_result/'+project+'_results/'+str(result_order)+'.Mapping/'+str(result_order)+'.1.MapStat/'+str(result_order)+'.1MapStat.README.txt'))
os.chdir('..')

assert not os.system('mkdir '+str(result_order)+'.2.MapReg')
mr_order=str(result_order)
os.chdir(''+str(result_order)+'.2.MapReg')
for eachsample in samples:
	assert not os.system('cp %s %s' % (root_dir+'/QC_DGE/'+eachsample+'/files/'+eachsample+'.MR.png',eachsample+'.Mapped_Region.png'))
	assert not os.system('cp %s %s' % (root_dir+'/QC_DGE/'+eachsample+'/files/'+eachsample+'.MR.pdf',eachsample+'.Mapped_Region.pdf'))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/MapReg.README.txt %s' % (root_dir+'/'+project+'_DGE_result/'+project+'_results/'+str(result_order)+'.Mapping/'+str(result_order)+'.2.MapReg/'+str(result_order)+'.2MappedRegion.README.txt'))
os.chdir('..')

assert not os.system('mkdir '+str(result_order)+'.3.ChrDen')
cd_order=str(result_order)
os.chdir(''+str(result_order)+'.3.ChrDen')
for eachsample in samples:
	assert not os.system('cp %s .' % (root_dir+'/Curve_DGE/density/'+eachsample+'/'+eachsample+'.density.png'))
	assert not os.system('cp %s .' % (root_dir+'/Curve_DGE/density/'+eachsample+'/'+eachsample+'.density.pdf'))
	assert not os.system('cp %s .' % (root_dir+'/Curve_DGE/density/'+eachsample+'/'+eachsample+'.ReadsVSlength.png'))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/ChrDen.README.txt %s' % (root_dir+'/'+project+'_DGE_result/'+project+'_results/'+str(result_order)+'.Mapping/'+str(result_order)+'.3.ChrDen/'+str(result_order)+'.3ChrDen.README.txt'))
os.chdir('..')

assert not os.system('mkdir '+str(result_order)+'.4.IGV')
igv_order=str(result_order)
os.chdir(''+str(result_order)+'.4.IGV')
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/Report/IGVQuickStart.pdf %s' % (root_dir+'/'+project+'_DGE_result/'+project+'_results/'+str(result_order)+'.Mapping/'+str(result_order)+'.4.IGV'))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/IGV.README.txt %s' % (root_dir+'/'+project+'_DGE_result/'+project+'_results/'+str(result_order)+'.Mapping/'+str(result_order)+'.4.IGV/'+str(result_order)+'.4IGV.README.txt'))
os.chdir(root_dir+'/'+project+'_DGE_result/'+project+'_results')
'''

if set([1]).issubset(includes):
	code+='''
result_order+=1
assert not os.system('mkdir '+str(result_order)+'.GeneExprQuatification')
HTSeq_order=str(result_order)
os.chdir(''+str(result_order)+'.GeneExprQuatification')
assert not os.system('cp %s .' % (root_dir+'/Diff_DGE/Diff/rpkm.xls'))
assert not os.system('cp %s .' % (root_dir+'/Diff_DGE/Diff/readcount.xls'))
assert not os.system('cp %s .' % (root_dir+'/Diff_DGE/Diff/rpkm.stat.xls'))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/Quatification.README.txt %s' % (root_dir+'/'+project+'_DGE_result/'+project+'_results/'+str(result_order)+'.GeneExprQuatification/'+str(result_order)+'.GeneExprQuatification.README.txt'))
os.chdir(root_dir+'/'+project+'_DGE_result/'+project+'_results')
'''
if set([1,3]).issubset(includes):
	code+='''
result_order+=1
assert not os.system('mkdir '+str(result_order)+'.AdvancedQC')
os.chdir(''+str(result_order)+'.AdvancedQC')
assert not os.system('mkdir '+str(result_order)+'.1.SaturationCurve')
os.chdir(''+str(result_order)+'.1.SaturationCurve')
for eachsample in samples:
	assert not os.system('cp %s %s' % (root_dir+'/Curve_DGE/SaturationCurve/'+eachsample+'/'+eachsample+'.saturation.png',eachsample+'.Saturation_curve.png'))
	assert not os.system('cp %s %s' % (root_dir+'/Curve_DGE/SaturationCurve/'+eachsample+'/'+eachsample+'.saturation.pdf',eachsample+'.Saturation_curve.pdf'))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/SatCurve.README.txt %s' % (root_dir+'/'+project+'_DGE_result/'+project+'_results/'+str(result_order)+'.AdvancedQC/'+str(result_order)+'.1.SaturationCurve/'+str(result_order)+'.1SaturationCurve.README.txt'))
os.chdir('..')

if flag_uniform:
	assert not os.system('mkdir '+str(result_order)+'.2.MeanCoverage')
	os.chdir(''+str(result_order)+'.2.MeanCoverage')
	for eachsample in samples:
		assert not os.system('cp %s %s' % (root_dir+'/QC_DGE/'+eachsample+'/files/'+eachsample+'.MC.png',eachsample+'.Mean_coverage_distribution.png'))
		assert not os.system('cp %s %s' % (root_dir+'/QC_DGE/'+eachsample+'/files/'+eachsample+'.MC.pdf',eachsample+'.Mean_coverage_distribution.pdf'))
	assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/MeanCov.README.txt %s' % (root_dir+'/'+project+'_DGE_result/'+project+'_results/'+str(result_order)+'.AdvancedQC/'+str(result_order)+'.2.MeanCoverage/'+str(result_order)+'.2MeanCoverage.README.txt'))
os.chdir(root_dir+'/'+project+'_DGE_result/'+project+'_results')

if flag_uniform:
	os.chdir(''+str(result_order)+'.AdvancedQC')
	assert not os.system('mkdir '+str(result_order)+'.3.Correlation')
	os.chdir(''+str(result_order)+'.3.Correlation')
	for png in glob.iglob(root_dir+'/Diff_DGE/Diff/corr_plot/*.png'):
		assert not os.system('cp %s .' % png)
	for pdf in glob.iglob(root_dir+'/Diff_DGE/Diff/corr_plot/*.pdf'):
		assert not os.system('cp %s .' % pdf)
	assert not os.system( 'cp %s .' % (root_dir+ '/Diff_DGE/Diff/corr_plot/cor_pearson.xls' ) )
	assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/DupCorr.README.txt %s' % (root_dir+'/'+project+'_DGE_result/'+project+'_results/'+str(result_order)+'.AdvancedQC/'+str(result_order)+'.3.Correlation/'+str(result_order)+'.3Correlation.README.txt'))
	os.chdir(root_dir+'/'+project+'_DGE_result/'+project+'_results')
else:
	os.chdir(''+str(result_order)+'.AdvancedQC')
	assert not os.system('mkdir '+str(result_order)+'.2.Correlation')
	os.chdir(''+str(result_order)+'.2.Correlation')
	for png in glob.iglob(root_dir+'/Diff_DGE/Diff/corr_plot/*.png'):
		assert not os.system('cp %s .' % png)
	for pdf in glob.iglob(root_dir+'/Diff_DGE/Diff/corr_plot/*.pdf'):
		assert not os.system('cp %s .' % pdf)
	assert not os.system( 'cp %s .' % (root_dir+ '/Diff_DGE/Diff/corr_plot/cor_pearson.xls' ) )
	assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/DupCorr.README.txt %s' % (root_dir+'/'+project+'_DGE_result/'+project+'_results/'+str(result_order)+'.AdvancedQC/'+str(result_order)+'.2.Correlation/'+str(result_order)+'.2Correlation.README.txt'))
	os.chdir(root_dir+'/'+project+'_DGE_result/'+project+'_results')
'''

if set([1,2]).issubset(includes):
	code+='''
groupnames=groupname.split(',')
compares=compare.split(',')
compare_names=[]
for each in compares:
	temp=each.split(':')
	compare_names.append(groupnames[int(temp[0])-1]+'vs'+groupnames[int(temp[1])-1])
#compare_names=','.join(compare_names)
venn_cluster_vs_names=[]
flag_venn=False
for each in venn_cluster_name.split(','):
	if each.count(':') >= 2 and each.count(':') <= 5:
		flag_venn=True
		venn_cluster_vs_names.append(each.replace(':','vs'))

result_order+=1
assert not os.system('mkdir '+str(result_order)+'.DiffExprAnalysis')
os.chdir(''+str(result_order)+'.DiffExprAnalysis')
assert not os.system('mkdir '+str(result_order)+'.1.GeneExpContrast')
os.chdir(''+str(result_order)+'.1.GeneExpContrast')
assert not os.system('cp %s .' % (root_dir+'/Diff_DGE/Diff/boxplot.pdf'))
assert not os.system('cp %s .' % (root_dir+'/Diff_DGE/Diff/boxplot.png'))
assert not os.system('cp %s .' % (root_dir+'/Diff_DGE/Diff/density.pdf'))
assert not os.system('cp %s .' % (root_dir+'/Diff_DGE/Diff/density.png'))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/ExpLev.README.txt %s' % (root_dir+'/'+project+'_DGE_result/'+project+'_results/'+str(result_order)+'.DiffExprAnalysis/'+str(result_order)+'.1.GeneExpContrast/'+str(result_order)+'.1.GeneExpContrast.README.txt'))
os.chdir('..')

assert not os.system('mkdir '+str(result_order)+'.2.DEGsList')
os.chdir(''+str(result_order)+'.2.DEGsList')
for each in compare_names:
	assert not os.system('cp %s %s' % (root_dir+'/Diff_DGE/Diff/DiffGeneList/'+each+'.diffgene.genename.xls',each+'.DEG.xls'))
	assert not os.system('cp %s %s' % (root_dir+'/Diff_DGE/Diff/DiffGeneList/'+each+'.diffgene_up.genename.xls',each+'.DEG_up.xls'))
	assert not os.system('cp %s %s' % (root_dir+'/Diff_DGE/Diff/DiffGeneList/'+each+'.diffgene_down.genename.xls',each+'.DEG_down.xls'))
	assert not os.system('cp %s %s' % (root_dir+'/Diff_DGE/Diff/Diff_Gene_Seq/'+each+'.diffgene.seq',each+'.DEG.fasta'))
	assert not os.system('cp %s . ' % (root_dir+'/Diff_DGE/Diff/'+each+'/'+each+'.Differential_analysis_results.xls'))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/DiffList.README.txt %s' % (root_dir+'/'+project+'_DGE_result/'+project+'_results/'+str(result_order)+'.DiffExprAnalysis/'+str(result_order)+'.2.DEGsList/'+str(result_order)+'.2DEGsList.README.txt'))
os.chdir('..')

assert not os.system('mkdir '+str(result_order)+'.3.DEGsFilter')
os.chdir(''+str(result_order)+'.3.DEGsFilter')
for each in compare_names:
	assert not os.system('cp %s %s' % (root_dir+'/Diff_DGE/Diff/'+each+'/'+each+'.volcano.pdf',each+'.Volcanoplot.pdf'))
	assert not os.system('cp %s %s' % (root_dir+'/Diff_DGE/Diff/'+each+'/'+each+'.volcano.png',each+'.Volcanoplot.png'))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/DiffFilt.README.txt %s' % (root_dir+'/'+project+'_DGE_result/'+project+'_results/'+str(result_order)+'.DiffExprAnalysis/'+str(result_order)+'.3.DEGsFilter/'+str(result_order)+'.3DEGsFilter.README.txt'))
os.chdir('..')

assert not os.system('mkdir '+str(result_order)+'.4.DEGcluster')
os.chdir(''+str(result_order)+'.4.DEGcluster')
assert not os.system('cp %s Union_for_cluster.xls' % ( root_dir+'/Diff_DGE/Diff/union_for_cluster' ))
assert not os.system('cp %s .' % (root_dir+'/Diff_DGE/Diff/heatCluster.pdf'))
assert not os.system('cp %s .' % (root_dir+'/Diff_DGE/Diff/heatCluster.png'))
os.system('cp %s .' % (root_dir+'/Diff_DGE/Diff/heatCluster.detail.pdf'))
assert not os.system('cp %s .' % (root_dir+'/Diff_DGE/Diff/union_for_cluster_hclust/h_show_plots.png'))
assert not os.system('cp %s .' % (root_dir+'/Diff_DGE/Diff/union_for_cluster_kmeans/kmeans_cluster_plots.pdf'))
assert not os.system('cp %s .' % (root_dir+'/Diff_DGE/Diff/union_for_cluster_som/som_cluster_plots.pdf'))
assert not os.system('cp %s %s' % (root_dir+'/Diff_DGE/Diff/union_for_cluster_hclust/h_subcluster_plots.pdf','h_cluster_plots.pdf'))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/Cluster.README.txt %s' % (root_dir+'/'+project+'_DGE_result/'+project+'_results/'+str(result_order)+'.DiffExprAnalysis/'+str(result_order)+'.4.DEGcluster/'+str(result_order)+'.4DEGcluster.README.txt'))
assert not os.system('mkdir '+'Subcluster')
assert not os.system('mkdir '+'Subcluster/kmeans')
assert not os.system('mkdir '+'Subcluster/hclust')
assert not os.system('mkdir '+'Subcluster/som')
assert not os.system('ls %s|while read a;do cut -f 1 $a>%s${a##*/}.xls;done' % (root_dir+'/Diff_DGE/Diff/union_for_cluster_kmeans/subcluster*','Subcluster/kmeans/'))
assert not os.system('ls %s|while read a;do cut -f 1 $a>%s${a##*/}.xls;done' % (root_dir+'/Diff_DGE/Diff/union_for_cluster_hclust/subcluster*','Subcluster/hclust/'))
assert not os.system('ls %s|while read a;do cut -f 1 $a>%s${a##*/}.xls;done' % (root_dir+'/Diff_DGE/Diff/union_for_cluster_som/subcluster*','Subcluster/som/'))
os.chdir('..')

if flag_venn:
	assert not os.system('mkdir '+str(result_order)+'.5.VennDiagram')
	os.chdir(''+str(result_order)+'.5.VennDiagram')
	for each in venn_cluster_vs_names:
		assert not os.system('cp -r %s .' % (root_dir+'/Diff_DGE/Diff/Venn/'+each))
	assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/Venn.README.txt %s' % (root_dir+'/'+project+'_DGE_result/'+project+'_results/'+str(result_order)+'.DiffExprAnalysis/'+str(result_order)+'.5.VennDiagram/'+str(result_order)+'.5VennDiagram.README.txt'))
os.chdir(root_dir+'/'+project+'_DGE_result/'+project+'_results')
'''

if set([1,2,4]).issubset(includes):
	code+='''
result_order+=1
assert not os.system('mkdir '+str(result_order)+'.DEG_GOEnrichment')
os.chdir(''+str(result_order)+'.DEG_GOEnrichment')
assert not os.system('mkdir '+str(result_order)+'.1.DEG_GOList')
os.chdir(''+str(result_order)+'.1.DEG_GOList')
for each in compare_names:
	assert not os.system('cp %s %s' % (root_dir+'/GOSeq_DGE/'+each+'/'+each+'.GO_enrichment_result_up_down.xls',each+'.DEG_GO_enrichment_result_up_down.xls'))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/GOList.README.txt %s' % (root_dir+'/'+project+'_DGE_result/'+project+'_results/'+str(result_order)+'.DEG_GOEnrichment/'+str(result_order)+'.1.DEG_GOList/'+str(result_order)+'.1DEG_GOList.README.txt'))
os.chdir('..')

assert not os.system('mkdir '+str(result_order)+'.2.DAG')
os.chdir(''+str(result_order)+'.2.DAG')
for each in compare_names:
	os.system('cp %s .' % (root_dir+'/GOSeq_DGE/'+each+'/'+each+'.GObp_DAG.png'))
	os.system('cp %s .' % (root_dir+'/GOSeq_DGE/'+each+'/'+each+'.GObp_DAG.pdf'))
	os.system('cp %s .' % (root_dir+'/GOSeq_DGE/'+each+'/'+each+'.GOmf_DAG.png'))
	os.system('cp %s .' % (root_dir+'/GOSeq_DGE/'+each+'/'+each+'.GOmf_DAG.pdf'))
	os.system('cp %s .' % (root_dir+'/GOSeq_DGE/'+each+'/'+each+'.GOcc_DAG.png'))
	os.system('cp %s .' % (root_dir+'/GOSeq_DGE/'+each+'/'+each+'.GOcc_DAG.pdf'))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/DAG.README.txt %s' % (root_dir+'/'+project+'_DGE_result/'+project+'_results/'+str(result_order)+'.DEG_GOEnrichment/'+str(result_order)+'.2.DAG/'+str(result_order)+'.2DAG.README.txt'))
os.chdir('..')

assert not os.system('mkdir '+str(result_order)+'.3.BAR')
os.chdir(''+str(result_order)+'.3.BAR')
for each in compare_names:
	assert not os.system('cp %s .' % (root_dir+'/GOSeq_DGE/'+each+'/'+each+'.bar_graph.png'))
	assert not os.system('cp %s .' % (root_dir+'/GOSeq_DGE/'+each+'/'+each+'.bar_graph.pdf'))
	assert not os.system('cp %s .' % (root_dir+'/GOSeq_DGE/'+each+'/'+each+'.bar_graph_updown.png'))
	assert not os.system('cp %s .' % (root_dir+'/GOSeq_DGE/'+each+'/'+each+'.bar_graph_updown.pdf'))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/BAR.README.txt %s' % (root_dir+'/'+project+'_DGE_result/'+project+'_results/'+str(result_order)+'.DEG_GOEnrichment/'+str(result_order)+'.3.BAR/'+str(result_order)+'.3BAR.README.txt'))
os.chdir(root_dir+'/'+project+'_DGE_result/'+project+'_results')
'''

if set([1,2,5]).issubset(includes):
	code+='''
result_order+=1
assert not os.system('mkdir '+str(result_order)+'.DEG_KEGGEnrichment')
os.chdir(''+str(result_order)+'.DEG_KEGGEnrichment')
assert not os.system('mkdir '+str(result_order)+'.1.DEG_KEGGList')
os.chdir(''+str(result_order)+'.1.DEG_KEGGList')
for each in compare_names:
	assert not os.system('cp %s %s' % (root_dir+'/KOBAS_DGE/'+each+'/add.'+each+'.identify.xls',each+'.DEG_KEGG_pathway_enrichment_result.xls'))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/KEGGList.README.txt %s' % (root_dir+'/'+project+'_DGE_result/'+project+'_results/'+str(result_order)+'.DEG_KEGGEnrichment/'+str(result_order)+'.1.DEG_KEGGList/'+str(result_order)+'.1DEG_KEGGList.README.txt'))
os.chdir('..')

assert not os.system('mkdir '+str(result_order)+'.2.DEG_KEGGScat')
os.chdir(''+str(result_order)+'.2.DEG_KEGGScat')
for each in compare_names:
	os.system('cp %s .' % (root_dir+'/KOBAS_DGE/'+each+'/'+each+'.DEG_enriched_KEGG_pathway_scatterplot.png'))
	os.system('cp %s .' % (root_dir+'/KOBAS_DGE/'+each+'/'+each+'.DEG_enriched_KEGG_pathway_scatterplot.pdf'))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/KEGGScat.README.txt %s' % (root_dir+'/'+project+'_DGE_result/'+project+'_results/'+str(result_order)+'.DEG_KEGGEnrichment/'+str(result_order)+'.2.DEG_KEGGScat/'+str(result_order)+'.2DEG_KEGGScat.README.txt'))
os.chdir('..')

assert not os.system('mkdir '+str(result_order)+'.3.DEG_KEGGPath')
os.chdir(''+str(result_order)+'.3.DEG_KEGGPath')
for each in compare_names:
	assert not os.system('mkdir '+each)
	assert not os.system('cp %s %s' % (root_dir+'/KOBAS_DGE/'+each+'/'+each+'.html',each))
	assert not os.system('cp -r %s %s' % (root_dir+'/KOBAS_DGE/'+each+'/src',each))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/KEGGPath.README.txt %s' % (root_dir+'/'+project+'_DGE_result/'+project+'_results/'+str(result_order)+'.DEG_KEGGEnrichment/'+str(result_order)+'.3.DEG_KEGGPath/'+str(result_order)+'.3DEG_KEGGPath.README.txt'))
os.chdir(root_dir+'/'+project+'_DGE_result/'+project+'_results')
'''

if set([1,2,6]).issubset(includes):
	code+='''
result_order+=1
assert not os.system('mkdir '+str(result_order)+'.DEG_PPI')
os.chdir(''+str(result_order)+'.DEG_PPI')
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/Report/CytoscapeQuickStart.pdf .')
for each in compare_names:
	assert not os.system('cp %s .' % (root_dir+'/PPI_DGE/'+each+'/'+each+'.ppi.txt'))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/ProteinNetwork.README.txt %s' % (root_dir+'/'+project+'_DGE_result/'+project+'_results/'+str(result_order)+'.DEG_PPI/'+str(result_order)+'.DEG_PPI.README.txt'))
os.chdir(root_dir+'/'+project+'_DGE_result/'+project+'_results')
'''

code+='''

print project+' report generating...'
os.chdir(root_dir)
os.chdir(project+'_DGE_result/'+project+'_results')
os.system('tree -d -v -L 3 -o ../'+project+'_report/DirectoryTree.html -H ../'+project+'_results')
DirectoryTree = root_dir+'/'+project+'_DGE_result/'+project+'_report/DirectoryTree.html'
NewDirectoryTree = DirectoryTree + '.new'
os.system( 'perl /PUBLIC/source/RNA/RefRNA/Report/bin/result_tree.pl %s %s' % (DirectoryTree, NewDirectoryTree) )
os.system('mv %s %s' % (NewDirectoryTree, DirectoryTree))
#os.system('perl /PUBLIC/source/RNA/RefRNA/TransRef/Report/dtree2jtree.pl %s > %s' % (root_dir+'/'+project+'_result/report/DirectoryTree.html',root_dir+'/'+project+'_result/report/src/tables/jtree'))
os.chdir(root_dir)
#generate the TransRef html report (separate the processes from above for clearness)
os.chdir(root_dir+'/'+project+'_DGE_result/'+project+'_report')
assert not os.system('cp -r /PUBLIC/source/RNA/RefRNA/Report/template/DGE/src .')
#os.system('perl /PUBLIC/source/RNA/RefRNA/TransRef/Report/dtree2jtree.pl %s > %s' % (root_dir+'/'+project+'_DGE_result/'+project+'_report/DirectoryTree.html',root_dir+'/'+project+'_DGE_result/'+project+'_report/src/tables/jtree'))

render={}
render['title']=project
render['flag_qc']=False
render['flag_tophat']=False
render['flag_htseq']=False
render['flag_rnaseqqc']=False
render['diffexp']=False
render['flag_go']=False
render['flag_kegg']=False
render['flag_ppi']=False

figure_order=0
table_order=0
html_order=0

render['flag_qc']=True
html_order+=1	#for raw data stat
html_order+=1	#for qc
render['html_order_qc']=str(html_order)
figure_order+=1
render['figure_order_qc1']=str(figure_order)
figure_order+=1
render['figure_order_qc2']=str(figure_order)
figure_order+=1
render['figure_order_qc3']=str(figure_order)
table_order+=1
render['table_order_qc']=str(table_order)

assert not os.system('mkdir ./src/pictures')
assert not os.system('mkdir ./src/pictures/2.1Error')
#assert not os.system('mkdir ./src/pictures/2.2GC')
#assert not os.system('mkdir ./src/pictures/2.3Filter')
#assert not os.system('mkdir ./src/pictures/2.4.DataTable')
render['figure2_1']=[]
render['figure2_2']=[]
render['figure2_3']=[]
render['table2_4']=[]
for eachsample in samples:
	assert not os.system('convert -resize 600 %s ./src/pictures/%s.error_rate_distribution.png' % (root_dir+'/QC_DGE/'+eachsample+'/files/'+eachsample+'.Error.png',eachsample))
	assert not os.system('convert -resize 90 %s ./src/pictures/%s.error_rate_distribution.JPEG' % (root_dir+'/QC_DGE/'+eachsample+'/files/'+eachsample+'.Error.png',eachsample))
	render['figure2_1'].append(["'"+'src/pictures/'+eachsample+'.error_rate_distribution.png'+"'","'"+'src/pictures/'+eachsample+'.error_rate_distribution.JPEG'+"'"])
	assert not os.system('convert -resize 600 %s ./src/pictures/%s.GC_content_distribution.png' % (root_dir+'/QC_DGE/'+eachsample+'/files/'+eachsample+'.GC.png',eachsample))
	assert not os.system('convert -resize 90 %s ./src/pictures/%s.GC_content_distribution.JPEG' % (root_dir+'/QC_DGE/'+eachsample+'/files/'+eachsample+'.GC.png',eachsample))
	render['figure2_2'].append(["'"+'src/pictures/'+eachsample+'.GC_content_distribution.png'+"'","'"+'src/pictures/'+eachsample+'.GC_content_distribution.JPEG'+"'"])
	assert not os.system('convert -resize 600 %s ./src/pictures/%s.raw_reads_classification.png' % (root_dir+'/QC_DGE/'+eachsample+'/files/'+eachsample+'.pie3d.png',eachsample))
	assert not os.system('convert -resize 90 %s ./src/pictures/%s.raw_reads_classification.JPEG' % (root_dir+'/QC_DGE/'+eachsample+'/files/'+eachsample+'.pie3d.png',eachsample))
	render['figure2_3'].append(["'"+'src/pictures/'+eachsample+'.raw_reads_classification.png'+"'","'"+'src/pictures/'+eachsample+'.raw_reads_classification.JPEG'+"'"])
	for eachLine in open(root_dir+'/QC_DGE/'+eachsample+'/files/dataTable'):
		if eachLine.strip() != '':
			render['table2_4'].append(eachLine.strip().split())

render['flag_tophat']=True
render['flag_flow']=True
html_order+=1
render['html_order_tophat']=str(html_order)

table_order+=1
render['table_order_tophat']=str(table_order)
render['table3_1']=[]
temp = linecache.getline(root_dir+'/'+project+'_DGE_result/'+project+'_results/3.Mapping/3.1.MapStat/MapStat.xls',1)
tophat_sample = temp.strip().split('\t')[1:]
render['tophat_head']='</th><th>'.join(tophat_sample)
cont=open(root_dir+'/'+project+'_DGE_result/'+project+'_results/3.Mapping/3.1.MapStat/MapStat.xls').readlines()
cont.pop(0)
for eachLine in cont:
	if (not eachLine.startswith('Sample')) and eachLine.strip() != '':
		temp=eachLine.split('\\t')
		tophat='</td><td>'.join(temp)
		render['table3_1'].append(tophat)

figure_order+=1
render['figure_order_tophat2']=str(figure_order)
figure_order+=1
render['figure_order_tophat3']=str(figure_order)
figure_order+=1
render['figure_order_tophat4']=str(figure_order)
render['figure3_2']=[]
render['figure3_3']=[]
assert not os.system('mkdir ./src/pictures/3.2MapReg')
assert not os.system('mkdir ./src/pictures/3.3ChrDen')
for eachsample in samples:
	assert not os.system('convert -resize 600 %s ./src/pictures/3.2MapReg/%s.Mapped_Region.png' % (root_dir+'/QC_DGE/'+eachsample+'/files/'+eachsample+'.MR.png',eachsample))
	assert not os.system('convert -resize 90 %s ./src/pictures/3.2MapReg/%s.Mapped_Region.JPEG' % (root_dir+'/QC_DGE/'+eachsample+'/files/'+eachsample+'.MR.png',eachsample))
	render['figure3_2'].append(["'"+'src/pictures/3.2MapReg/'+eachsample+'.Mapped_Region.png'+"'","'"+'src/pictures/3.2MapReg/'+eachsample+'.Mapped_Region.JPEG'+"'"])
	assert not os.system('convert -resize 600 %s ./src/pictures/3.3ChrDen/%s.density.png' % (root_dir+'/Curve_DGE/density/'+eachsample+'/'+eachsample+'.density.png',eachsample))
	assert not os.system('convert -resize 90 %s ./src/pictures/3.3ChrDen/%s.density.JPEG' % (root_dir+'/Curve_DGE/density/'+eachsample+'/'+eachsample+'.density.png',eachsample))
	render['figure3_3'].append(["'"+'src/pictures/3.3ChrDen/'+eachsample+'.density.png'+"'","'"+'src/pictures/3.3ChrDen/'+eachsample+'.density.JPEG'+"'"])
	assert not os.system('convert -resize 600 %s ./src/pictures/3.3ChrDen/%s.ReadsVSlength.png' % (root_dir+'/Curve_DGE/density/'+eachsample+'/'+eachsample+'.ReadsVSlength.png',eachsample))
	assert not os.system('convert -resize 90 %s ./src/pictures/3.3ChrDen/%s.ReadsVSlength.JPEG' % (root_dir+'/Curve_DGE/density/'+eachsample+'/'+eachsample+'.ReadsVSlength.png',eachsample))
	render['figure3_3'].append(["'"+'src/pictures/3.3ChrDen/'+eachsample+'.ReadsVSlength.png'+"'","'"+'src/pictures/3.3ChrDen/'+eachsample+'.ReadsVSlength.JPEG'+"'"])
'''

if set([1]).issubset(includes):
	code+='''
render['flag_htseq']=True
html_order+=1
render['html_order_htseq']=str(html_order)

table_order+=1
render['table_order_htseq1']=str(table_order)
render['table4_1']=[]
render['rpkm_stat_head']='</th><th>'.join(samples)
cont=open(root_dir+'/Diff_DGE/Diff/rpkm.stat.xls').readlines()
cont.pop(0)
for eachLine in cont:
	if (not eachLine.startswith('RPKM')) and eachLine.strip() !='':
		temp=eachLine.split('\\t')
		rpkm_stat='</td><td>'.join(temp)
		render['table4_1'].append(rpkm_stat)

table_order+=1
render['table_order_htseq2']=str(table_order)
render['table4_2']=[]
k=4
render['rpkm_head']='</th><th>'.join(samples)
cont=open(root_dir+'/Diff_DGE/Diff/rpkm.xls').readlines()
cont.pop(0)
for eachLine in cont:
	if (not eachLine.startswith('geneID')) and eachLine.strip() !='':
		temp=eachLine.split('\\t')
		rpkm='</td><td>'.join(temp)
		render['table4_2'].append(rpkm)
		k-=1
	if k == 0:
		break
'''

if set([1,3]).issubset(includes):
	code+='''
render['flag_rnaseqqc']=True
html_order+=1
render['html_order_rnaseqqc']=str(html_order)

figure_order+=1
render['figure_order_rnaseqqc1']=str(figure_order)
figure_order+=1
render['figure_order_rnaseqqc2']=str(figure_order)
render['figure5_1']=[]
render['figure5_2']=[]
for eachsample in samples:
	assert not os.system('convert -resize 600 %s ./src/pictures/%s.saturation.png' % (root_dir+'/Curve_DGE/SaturationCurve/'+eachsample+'/'+eachsample+'.saturation.png',eachsample))
	assert not os.system('convert -resize 90 %s ./src/pictures/%s.saturation.JPEG' % (root_dir+'/Curve_DGE/SaturationCurve/'+eachsample+'/'+eachsample+'.saturation.png',eachsample))
	render['figure5_1'].append(["'"+'src/pictures/'+eachsample+'.saturation.png'+"'","'"+'src/pictures/'+eachsample+'.saturation.JPEG'+"'"])
	assert not os.system('convert -resize 600 %s ./src/pictures/%s.Mean_coverage_distribution.png' % (root_dir+'/QC_DGE/'+eachsample+'/files/'+eachsample+'.MC.png',eachsample))
	assert not os.system('convert -resize 90 %s ./src/pictures/%s.Mean_coverage_distribution.JPEG' % (root_dir+'/QC_DGE/'+eachsample+'/files/'+eachsample+'.MC.png',eachsample))
	render['figure5_2'].append(["'"+'src/pictures/'+eachsample+'.Mean_coverage_distribution.png'+"'","'"+'src/pictures/'+eachsample+'.Mean_coverage_distribution.JPEG'+"'"])
'''

if set([1,2,3]).issubset(includes):
	code+='''
render['flag_rnaseqqc_cor']=True
figure_order+=1
render['figure_order_rnaseqqc_cor']=str(figure_order)
render['figure5_3']=[]
assert not os.system('convert -resize 600 %s %s' % (root_dir+'/Diff_DGE/Diff/corr_plot/cor_pearson.png','./src/pictures/cor_pearson.png'))
assert not os.system('convert -resize 90 %s %s' % (root_dir+'/Diff_DGE/Diff/corr_plot/cor_pearson.png','./src/pictures/cor_pearson.JPEG'))
render['figure5_3'].append(["'"+'./src/pictures/cor_pearson.png'+"'","'"+'./src/pictures/cor_pearson.JPEG'+"'"])
for png in glob.iglob(root_dir+'/Diff_DGE/Diff/corr_plot/*.png'):
	if png !=root_dir+'/Diff_DGE/Diff/corr_plot/cor_pearson.png':
		name=re.search(r'(/.*/)(.*)\.png',png).group(2)
		assert not os.system('convert -resize 600 %s %s' % (png,'./src/pictures/'+name+'.png'))
		assert not os.system('convert -resize 90 %s %s' % (png,'./src/pictures/'+name+'.JPEG'))
		render['figure5_3'].append(["'"+'src/pictures/'+os.path.basename(png)+"'","'"+'src/pictures/'+name+'.JPEG'+"'"])
'''

if set([1,2]).issubset(includes):
	code+='''
render['flag_diffexp']=True
html_order+=1
render['html_order_diffexp']=str(html_order)

figure_order+=1
render['figure_order_diffexp1']=str(figure_order)
render['figure6_1']=[]
assert not os.system('convert -resize 600 %s ./src/pictures/boxplot.png' % (root_dir+'/Diff_DGE/Diff/boxplot.png'))
assert not os.system('convert -resize 90 %s ./src/pictures/boxplot.JPEG' % (root_dir+'/Diff_DGE/Diff/boxplot.png'))
assert not os.system('convert -resize 600 %s ./src/pictures/density.png' % (root_dir+'/Diff_DGE/Diff/density.png'))
assert not os.system('convert -resize 90 %s ./src/pictures/density.JPEG' % (root_dir+'/Diff_DGE/Diff/density.png'))
render['figure6_1'].append(["'"+'src/pictures/boxplot.png'+"'","'"+'src/pictures/boxplot.JPEG'+"'"])
render['figure6_1'].append(["'"+'src/pictures/density.png'+"'","'"+'src/pictures/density.JPEG'+"'"])

render['flag_venn']=flag_venn
table_order+=1
render['table_order_diffexp']=str(table_order)
render['table6_2']=[]
com=compare_names[0].split('vs')
render['com1']=com[0]
render['com2']=com[1]
k=4
for eachLine in open(root_dir+'/Diff_DGE/Diff/'+compare_names[0]+'/'+compare_names[0]+'.diffgene.xls'):
	if (not eachLine.startswith('Gene_id')) and eachLine.strip() != '':
		temp=eachLine.split('\\t')
		render['table6_2'].append(temp)
		k-=1
	if k == 0:
		break

figure_order+=1
render['figure_order_diffexp2']=str(figure_order)
render['figure6_3']=[]
for each in compare_names:
	assert not os.system('convert -resize 600 %s ./src/pictures/%s.Volcanoplot.png' % (root_dir+'/Diff_DGE/Diff/'+each+'/'+each+'.volcano.png',each))
	assert not os.system('convert -resize 90 %s ./src/pictures/%s.Volcanoplot.JPEG' % (root_dir+'/Diff_DGE/Diff/'+each+'/'+each+'.volcano.png',each))
	render['figure6_3'].append(["'"+'src/pictures/'+each+'.Volcanoplot.png'+"'","'"+'src/pictures/'+each+'.Volcanoplot.JPEG'+"'"])

figure_order+=1
render['figure_order_diffexp3']=str(figure_order)
render['figure6_4']=[]
assert not os.system('convert -resize 600 %s ./src/pictures/heatCluster.png' % (root_dir+'/Diff_DGE/Diff/heatCluster.png'))
assert not os.system('convert -resize 90 %s ./src/pictures/heatCluster.JPEG' % (root_dir+'/Diff_DGE/Diff/heatCluster.png'))
render['figure6_4'].append(["'"+'src/pictures/heatCluster.png'+"'","'"+'src/pictures/heatCluster.JPEG'+"'"])
assert not os.system('convert -resize 600 %s ./src/pictures/h_show_plots.png' % (root_dir+'/Diff_DGE/Diff/union_for_cluster_hclust/h_show_plots.png'))
assert not os.system('convert -resize 90 %s ./src/pictures/h_show_plots.JPEG' % (root_dir+'/Diff_DGE/Diff/union_for_cluster_hclust/h_show_plots.png'))
render['figure6_4'].append(["'"+'src/pictures/h_show_plots.png'+"'","'"+'src/pictures/h_show_plots.JPEG'+"'"])

if flag_venn:
	figure_order+=1
	render['figure_order_diffexp4']=str(figure_order)
	render['figure6_5']=[]
	for each in venn_cluster_vs_names:
		assert not os.system('convert -resize 600 %s ./src/pictures/%s.DEG_Venn_diagram.png' % (root_dir+'/Diff_DGE/Diff/Venn/'+each+'/'+each+'.venn.png',each))
		assert not os.system('convert -resize 90 %s ./src/pictures/%s.DEG_Venn_diagram.JPEG' % (root_dir+'/Diff_DGE/Diff/Venn/'+each+'/'+each+'.venn.png',each))
		render['figure6_5'].append(["'"+'src/pictures/'+each+'.DEG_Venn_diagram.png'+"'","'"+'src/pictures/'+each+'.DEG_Venn_diagram.JPEG'+"'"])
'''

if set([1,2,4]).issubset(includes):
	code+='''
render['flag_go']=True
html_order+=1
render['html_order_go']=str(html_order)

table_order+=1
render['table_order_go']=str(table_order)
render['table7_1']=[]
k=4
for eachLine in open(root_dir+'/GOSeq_DGE/'+compare_names[0]+'/'+compare_names[0]+'.GO_enrichment_result_up_down.xls'):
	if (not eachLine.startswith('GO_accession')) and eachLine.strip() != '':
		temp=eachLine.split('\\t')
		render['table7_1'].append(temp[:7])
		k-=1
	if k == 0:
		break

figure_order+=1
render['figure_order_go1']=str(figure_order)
render['figure7_2']=[]
figure_order+=1
render['figure_order_go2']=str(figure_order)
render['figure7_3']=[]

for each in compare_names:
	os.system('convert -resize 600 %s ./src/pictures/%s.DEG_Enriched_GO_classification.png' % (root_dir+'/GOSeq_DGE/'+each+'/'+each+'.bar_graph.png',each))
	os.system('convert -resize 90 %s ./src/pictures/%s.DEG_Enriched_GO_classification.JPEG' % (root_dir+'/GOSeq_DGE/'+each+'/'+each+'.bar_graph.png',each))
	render['figure7_2'].append(["'"+'src/pictures/'+each+'.DEG_Enriched_GO_classification.png'+"'","'"+'src/pictures/'+each+'.DEG_Enriched_GO_classification.JPEG'+"'"])
	os.system('convert -resize 600 %s ./src/pictures/%s.DEG_Enriched_GO_updown.png' % (root_dir+'/GOSeq_DGE/'+each+'/'+each+'.bar_graph_updown.png',each))
	os.system('convert -resize 90 %s ./src/pictures/%s.DEG_Enriched_GO_updown.JPEG' % (root_dir+'/GOSeq_DGE/'+each+'/'+each+'.bar_graph_updown.png',each))
	render['figure7_2'].append(["'"+'src/pictures/'+each+'.DEG_Enriched_GO_updown.png'+"'","'"+'src/pictures/'+each+'.DEG_Enriched_GO_updown.JPEG'+"'"])
	os.system('convert -resize 600 %s ./src/pictures/%s.DEG_Enriched_GO_bp_DAG.png' % (root_dir+'/GOSeq_DGE/'+each+'/'+each+'.GObp_DAG.png',each))
	os.system('convert -resize 90 %s ./src/pictures/%s.DEG_Enriched_GO_bp_DAG.JPEG' % (root_dir+'/GOSeq_DGE/'+each+'/'+each+'.GObp_DAG.png',each))
	render['figure7_3'].append(["'"+'src/pictures/'+each+'.DEG_Enriched_GO_bp_DAG.png'+"'","'"+'src/pictures/'+each+'.DEG_Enriched_GO_bp_DAG.JPEG'+"'"])
	os.system('convert -resize 600 %s ./src/pictures/%s.DEG_Enriched_GO_cc_DAG.png' % (root_dir+'/GOSeq_DGE/'+each+'/'+each+'.GOcc_DAG.png',each))
	os.system('convert -resize 90 %s ./src/pictures/%s.DEG_Enriched_GO_cc_DAG.JPEG' % (root_dir+'/GOSeq_DGE/'+each+'/'+each+'.GOcc_DAG.png',each))
	render['figure7_3'].append(["'"+'src/pictures/'+each+'.DEG_Enriched_GO_cc_DAG.png'+"'","'"+'src/pictures/'+each+'.DEG_Enriched_GO_cc_DAG.JPEG'+"'"])
	os.system('convert -resize 600 %s ./src/pictures/%s.DEG_Enriched_GO_mf_DAG.png' % (root_dir+'/GOSeq_DGE/'+each+'/'+each+'.GOmf_DAG.png',each))
	os.system('convert -resize 90 %s ./src/pictures/%s.DEG_Enriched_GO_mf_DAG.JPEG' % (root_dir+'/GOSeq_DGE/'+each+'/'+each+'.GOmf_DAG.png',each))
	render['figure7_3'].append(["'"+'src/pictures/'+each+'.DEG_Enriched_GO_mf_DAG.png'+"'","'"+'src/pictures/'+each+'.DEG_Enriched_GO_mf_DAG.JPEG'+"'"])
'''

if set([1,2,5]).issubset(includes):
	code+='''
render['flag_kegg']=True
html_order+=1
render['html_order_kegg']=str(html_order)

table_order+=1
render['table_order_kegg']=str(table_order)
flag=0

filelength=[]
for eachlen in compare_names:

        count=len(open(root_dir+'/KOBAS_DGE/'+eachlen+'/add.'+eachlen+'.identify.xls','r').readlines())
        filelength.append(count)
maxindex = filelength.index(max(filelength))
kcount = max(filelength) - 12

if(kcount >= 4):
        k=4
else:
        k=kcount

sig_map=''
kegg_png=[]
render['table8_1']=[]
for eachLine in open(root_dir+'/KOBAS_DGE/'+compare_names[maxindex]+'/add.'+compare_names[maxindex]+'.identify.xls'):
	if eachLine.startswith('#Term'):
		flag=1
		continue
	if flag == 1:
		if eachLine.strip() != '':
			temp=eachLine.split('\\t')
			if sig_map == '':
				sig_map=temp[2].strip()
				sig_map_png=root_dir+'/KOBAS_DGE/'+compare_names[maxindex]+'/src/'+sig_map+'.png'
				if os.path.isfile(sig_map_png):
					sig_map=temp[2].strip()
					sig_map_png=root_dir+'/KOBAS_DGE/'+compare_names[maxindex]+'/src/'+sig_map+'.png'
					kegg_png.append(sig_map_png)
					sig_map = ''
				else:
					sig_map = ''
			render['table8_1'].append(temp[:7])
			k-=1
		if k == 0:
			break

figure_order+=1
render['figure_order_kegg1']=str(figure_order)
render['figure8_2']=[]
render['figure20_lable']=[]
for each in compare_names:
	os.system('convert -resize 600 %s ./src/pictures/%s.DEG_enriched_KEGG_pathway_scatterplot.png' % (root_dir+'/KOBAS_DGE/'+each+'/'+each+'.DEG_enriched_KEGG_pathway_scatterplot.png',each))
	os.system('convert -resize 90 %s ./src/pictures/%s.DEG_enriched_KEGG_pathway_scatterplot.JPEG' % (root_dir+'/KOBAS_DGE/'+each+'/'+each+'.DEG_enriched_KEGG_pathway_scatterplot.png',each))
	render['figure8_2'].append(["'"+'src/pictures/'+each+'.DEG_enriched_KEGG_pathway_scatterplot.png'+"'","'"+'src/pictures/'+each+'.DEG_enriched_KEGG_pathway_scatterplot.JPEG'+"'"])
	render['figure20_lable'].append(each)
render['figure20_lable']=','.join(render['figure20_lable'])

figure_order+=1
render['figure_order_kegg2']=str(figure_order)
render['figure8_3']=[]
for png in kegg_png:
	name=re.search(r'(/.*/)(.*)\.png',png).group(2)
	assert not os.system('convert -resize 600 %s %s' % (png,'./src/pictures/'+name+'.png'))
	assert not os.system('convert -resize 90 %s %s' % (png,'./src/pictures/'+name+'.JPEG'))
	render['figure8_3'].append(["'"+'src/pictures/'+os.path.basename(png)+"'","'"+'src/pictures/'+name+'.JPEG'+"'"])
'''

if set([1,2,6]).issubset(includes):
	code+='''
render['flag_ppi']=True
html_order+=1
render['html_order_ppi']=str(html_order)

figure_order+=1
render['figure_order_ppi']=str(figure_order)
'''

code+='''
t=loader.get_template('DGE_report.html')
c=Context(render)
html=t.render(c)
open(project+'_Report.html','w').write(html)
assert not os.system("/PUBLIC/software/RNA/wkhtmltopdf/wkhtmltopdf --page-width 350mm --page-height 495mm -n --print-media-type --footer-center '[page] / [topage]' "+ project+'_Report.html '+project+'_Report.pdf')
assert not os.system("sed -i 's/#00//g' "+project+'_Report.pdf')
'''

code+='''
os.chdir(root_dir)
assert not os.system('perl /PUBLIC/source/RNA/RefRNA/TransRef/data_sum/data_sum.pl -n %s -p %s -ss %s -dir %s')
''' % (sample,project+'_DGE',ss,root_dir)

open('DGE_step4_Report.py','w').write(code)
#################################################################################################
#for data release
code='''
import os

fq='%s'
fa='%s'
gtf='%s'
goann='%s'
sample='%s'
root_dir='%s'
project='%s'
''' % (fq,fa,gtf,goann,sample,root_dir,project)

if argv['genenamefile']:
	code+='''
genenamefile='%s'
''' % (genenamefile)
else:
	code+='''
genenamefile='%s'
''' % (root_dir+'/Blast_DGE/Blast_Swissprot/diffgene_union.genenames')

code+='''
assert not os.system('mkdir data_release')
data_release=open(root_dir+'/data_release/data_release.sh','w')
data_release.write('mkdir %s\\n' % (root_dir+'/data_release/data_give'))
dir=root_dir+'/data_release/data_give'
data_release.write('echo "###############prepare SuppFiles#####################"\\n')
data_release.write('mkdir %s\\n' % (root_dir+'/'+project+'_DGE_result/'+project+'_results/0.SuppFiles'))
SuppFiles_dir=root_dir+'/'+project+'_DGE_result/'+project+'_results/0.SuppFiles'
data_release.write('cp /PUBLIC/source/RNA/RefRNA/TransRef/data_release/SuppFiles.README.txt %s\\n' % (SuppFiles_dir))
data_release.write('cp %s %s\\n' % (genenamefile,SuppFiles_dir+'/gene.description.xls'))
data_release.write('perl /PUBLIC/source/RNA/RefRNA/TransRef/data_release/GO_class.pl %s /PUBLIC/software/RNA/GOseq/gene_ontology.1_2.obo.ab %s\\n' % (goann,SuppFiles_dir+'/gene.goannot.xls'))
data_release.write('perl /PUBLIC/source/RNA/RefRNA/TransRef/data_release/extractInfo_from_GTF.pl -i %s -o %s\\n' % (gtf,SuppFiles_dir+'/gene.info.xls'))
data_release.write('perl /PUBLIC/source/RNA/RefRNA/DGE/Enrichment/extractcDNAfromFA.pl %s %s %s\\n' % (gtf,fa,SuppFiles_dir+'/gene.fasta'))
#data_release.write('gzip %s\\n' % (SuppFiles_dir+'/gene.fasta'))
data_release.write('echo "###############copy Novofinder#####################"\\n')
data_release.write('cp /PUBLIC/source/RNA/RefRNA/TransRef/data_release/Novofinder.exe %s\\n' % (root_dir+'/'+project+'_DGE_result/Novofinder.exe'))
data_release.write('cp /PUBLIC/source/RNA/RefRNA/TransRef/data_release/Novofinder_manual.pdf %s\\n' % (root_dir+'/'+project+'_DGE_result/'))
data_release.write('echo "###############tar Report#####################"\\n')
data_release.write('cd %s\\n' % (root_dir))
data_release.write('tar -cvzf %s %s\\n' % (project + '_DGE.tar.gz',project+'_DGE_result/'))
data_release.write('mv %s %s\\n' % (project + '_DGE.tar.gz',root_dir+'/data_release/data_give/'))
data_release.write('cd %s\\n' % (root_dir))
data_release.write('echo "###############prepare clean data#####################"\\n')
data_release.write('ln -s %s %s\\n' % (root_dir+'/QC_DGE/cleandata',dir+'/cleandata'))
data_release.write('echo "###############prepare raw data#####################"\\n')
data_release.write('mkdir %s\\n' % (dir+'/rawdata'))
data_release.write('cp %s %s' % (root_dir+'/QC_DGE/rd_md5.txt',dir+'/rawdata/\\n'))

samples=sample.split(',')
for eachsample in samples:
	data_release.write('ln -s %s %s\\n' % (root_dir+'/QC_DGE/'+eachsample+'/'+eachsample+'.fq.gz',dir+'/rawdata/'))

data_release.write('echo "###############prepare IGV_data#####################"\\n')
data_release.write('mkdir %s\\n' % (dir+'/IGV_data'))
data_release.write('cp /PUBLIC/source/RNA/RefRNA/TransRef/data_release/IGV_data.README.txt %s\\n' % (dir+'/IGV_data'))
data_release.write('cp %s %s\\n' % (fa,dir+'/IGV_data/genome.fasta'))
data_release.write('gzip %s\\n' % (dir+'/IGV_data/genome.fasta'))
data_release.write('ln -s %s %s\\n' % (root_dir+'/QC_DGE/bam',dir+'/IGV_data/bam'))
data_release.close()

os.chdir('data_release')
assert not os.system('qsub -V -cwd -l vf=1G -l p=1 data_release.sh')
os.chdir(root_dir)
'''

open('DGE_step5_data_release.py','w').write(code)
#################################################################################################
#for byebye
byebye=open(root_dir+'/byebye.sh','w')
byebye.write('sh /PUBLIC/source/RNA/RefRNA/TransRef/Pipeline/byebye_v2.1.sh -dir %s -sample %s -pipline DGE\n' % (root_dir,sample))
byebye.close()
