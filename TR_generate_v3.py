import sys
import os
import os.path
import re
import argparse
root_dir=os.getcwd()

##########################################################################
##########################################################################

#parse the arguments
parser = argparse.ArgumentParser(description="TR pipline v3.0")
parser.add_argument('--project',help='project name, maybe same with the name of root dir, which will be displayed in the final report title, [REQUIRED]',required=True)
parser.add_argument('--sample',help="sample names(sample1,sample2,...)  warning: order sensitive!",required=True)
parser.add_argument('--fq',help="the original directory of the raw fastq reads, [REQUIRED]",default=None)
parser.add_argument('--mapfile',help="mapfile files (mapfile1,mapfile2.The first column of mapfile is the NH number ,the second column of mapfile is the sample name)",default=None)
parser.add_argument('--raw_dir',help="the original directory of the raw fastq reads keep order in line with mapfile files (raw_dir1,raw_dir2)",default=None)
parser.add_argument('--ad',help="the original directory of the adapter list, ",default=None)
parser.add_argument('--generate_adapter',help='whether to generate the adpepter list y or n',choices=['y','n'],default='n')
parser.add_argument('--index',help='the index number of the fastq files,if you want to generate the adapter list,the parameter will be necessary,diffrent samples split by ,warning: order sensitive!',default=None)
parser.add_argument('--ss',help="strand specific, [REQUIRED]",choices=['no','yes','reverse'],required=True)
parser.add_argument('--number',help="chromosome number, [REQUIRED for density plot]",required=True)
parser.add_argument('--length',help="the length of sequenced reads,[defalut=100]",default='100')
parser.add_argument('--fa',help="the reference FASTA file, [REQUIRED for mapping]",required=True)
parser.add_argument('--gtf',help="the annotation GTF file, [REQUIRED for mapping]",required=True)
parser.add_argument('--group',help="sample classification, e.g. sample1/sample2,sample3, [REQUIRED]",default=None)
parser.add_argument('--groupname',help="group names, [group1,group2,... ]",default=None)
parser.add_argument('--compare',help="group comparison strategy 1:2,1:3,..., [REQUIRED]",default=None)
parser.add_argument('--venn',help="venn and cluster mode, defult is all groups, 1:2_1:3,1:3_2:3,... ",default=None)
parser.add_argument('--goann',help="gene to GO annotations file, [REQUIRED]",default=None)
parser.add_argument('--species',help="abbreviation for species, for kegg pathway. (note: all species: /PUBLIC/database/RNA//kobas2.0-data-20120208/seq_pep_v2), [defalut=kaas]",default='kaas')
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

argv = vars(parser.parse_args())
project=argv['project'].strip()
sample=argv['sample'].strip()
samples=sample.split(',')
samples_tmp=list(set(samples))
for s in samples:
  checkSample(s)
assert len(samples)==len(samples_tmp)
#-----------------------------------------------------------------------
mapfiles=[each.strip() for each in argv['mapfile'].strip().split(',') if each.strip() != '']
mapfile=' '.join(mapfiles)
assert not os.system('cat %s |sort -u |sed \'/^$/d\' >%s' % (mapfile,root_dir+'/libraryID'))
if argv['fq']:
	fq = argv['fq'].strip()
	fq = os.path.abspath(fq)
else:
	assert not os.system('mkdir raw_data')
	assert not os.system('perl /PUBLIC/source/RNA/RefRNA/ln_raw_data.pl %s %s pe raw_data' %(argv['mapfile'],argv['raw_dir']))
	fq = root_dir + '/raw_data'
for each in samples:
	fq_tmp1=fq+'/'+each+'_1.fq.gz'
	fq_tmp2=fq+'/'+each+'_2.fq.gz'
	assert os.path.isfile(fq_tmp1)
	assert os.path.isfile(fq_tmp2)
if argv['ad']:
	ad = argv['ad'].strip()
	ad = os.path.abspath(ad)
	for each in samples:
		ad_tmp1=ad+'/'+each+'_1.adapter.list.gz'
		ad_tmp2=ad+'/'+each+'_2.adapter.list.gz'
		assert os.path.isfile(ad_tmp1)
		assert os.path.isfile(ad_tmp2)
if argv['generate_adapter']:
	generate_adapter = argv['generate_adapter'].strip()
else:
	generate_adapter = 'n'
if argv['index']:
	index=argv['index'].strip()
	indexes=index.split(',')
	assert len(samples)==len(indexes)
if generate_adapter == 'y':
	if argv['ad'] !=None:
		print 'Error:  the parameters --ad and --generate_adapter are not consensus!\n'
		exit()
	if argv['index']==None:
		print 'Error:  the parameters --index and --generate_adapter are not consensus!\n'
		exit()
else:
	if argv['index'] != None:
		print 'Error:  the parameters --index and --generate_adapter are not consensus!\n'
		exit()
#-----------------------------------------------------------------------
all_content=set([1,2,3,4,5,6,7,8,9])
if argv['ex'] != None:
	excludes=argv['ex'].strip().strip(',').strip().split(',')
	excludes=[int(each.strip()) for each in excludes]
	for each1 in excludes:
		assert each1 in all_content
else:
	excludes=[] #list
if argv['group'].find(':') == -1:
	excludes.append(3)
if ( len(argv['group'].split(',') ) == 1) or (5 in excludes):
  excludes.extend([7,8,9])
includes=all_content-set(excludes) #set
#-----------------------------------------------------------------------
ss = argv['ss'].strip()
number = argv['number'].strip()
if argv['length']:
        length = argv['length'].strip()
else:
        length = '100'
fa = argv['fa'].strip()
fa = os.path.abspath(fa)
suffix_fa=['.1.bt2','.2.bt2','.3.bt2','.4.bt2','.rev.1.bt2','.rev.2.bt2']
for each in suffix_fa:
	fa_tmp=fa+each
	assert os.path.isfile(fa_tmp)
gtf = argv['gtf'].strip()
gtf = os.path.abspath(gtf)
assert os.path.isfile(gtf)
#-----------------------------------------------------------------------
if (set([1]).issubset(includes)):
	if argv['group'] == None:
		groups=samples
                groups_iter=samples
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
		assert len(each1) == 2
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
	assert len(groupnames)==len(list(set(groupnames)))
	for i,each in enumerate(groupnames):
		group_tmp=groups[i]
		group_tmp=group_tmp.replace(':',',')
	compare_name=[]
	for i,each in enumerate(compares):
		compare_tmp=each
		compare_tmp=compare_tmp.split(':')
		assert len(compare_tmp)==2
		compare_name=groupnames[int(compare_tmp[0])-1]+'vs'+groupnames[int(compare_tmp[1])-1]
	for i,each in enumerate(venns):
		venn_tmp=each
		venn_tmp=venn_tmp.split('_')
		for each2 in venn_tmp:
			assert each2 in compares
			venn_tmp2=each2.split(':')
			venn_name=groupnames[int(venn_tmp2[0])-1]+'vs'+groupnames[int(venn_tmp2[1])-1]
#-----------------------------------------------------------------------
goann = argv['goann'].strip()
goann = os.path.abspath(goann)
species = ppi_number = ppi_blast = ''
assert os.path.isfile(goann)
if set([4,5,8]).issubset(includes):
	if argv['species']:
		species = argv['species'].strip()
	else:
		species = 'kaas'
if set([4,5,9]).issubset(includes):
	if argv['ppi_number']:
        	ppi_number = argv['ppi_number'].strip()
	if argv['ppi_blast']:
	        ppi_blast = argv['ppi_blast'].strip()
	if argv['ppi_blast']:
		if argv['ppi_number'] == None:
			print 'Error:  the parameters --ppi_blast and --ppi_number are not consensus!\n'
			exit()
	else:
		if argv['ppi_number']:
			print 'Error:  the parameters --ppi_blast and --ppi_number are not consensus!\n'
			exit()
if argv['genenamefile']:
        genenamefile = argv['genenamefile'].strip()
	genenamefile = os.path.abspath(genenamefile)
	assert os.path.isfile(genenamefile)

##display all parameters##
display = open(root_dir + '/' + 'TR_command.txt','w')
display.write('project: %s\n' % (project))
display.write('sample: %s\n' % (sample))
display.write('fq: %s\n' % (fq))
display.write('adapter: ')
if argv['ad']:  
	display.write('%s\n' % (argv['ad'].rstrip()))
if argv['generate_adapter']:
	generate_adapter = argv['generate_adapter']
	display.write('generate_adapter: %s\n' % (generate_adapter))
if argv['index']:
	indexes = argv['index'].rstrip().split(',')
	for i,index_tmp in enumerate(samples):
		display.write('%s:\t%s\n' % (index_tmp,indexes[i]))
		display.write('\n')
display.write('length: %s\n' % (length))
display.write('fa: %s\n' % (fa))
display.write('gtf: %s\n' % (gtf))
display.write('\ngroups:\n')
for i,each in enumerate(groupnames):
	group_tmp=groups[i]
	group_tmp=group_tmp.replace(':',',')
	display.write('%s: %s\n' % (each,group_tmp))
display.write('\ncompare:\n')
compare_name=[]
for i,each in enumerate(compares):
	compare_tmp=each
	compare_tmp=compare_tmp.split(':')
	assert len(compare_tmp)==2
	compare_name=groupnames[int(compare_tmp[0])-1]+'vs'+groupnames[int(compare_tmp[1])-1]
	display.write('%s: \t' % (each))
	display.write('%s\n' % (compare_name))
display.write('\nvenn:\n')
for i,each in enumerate(venns):
	display.write('%s: \t' % (each))
	venn_tmp=each.split('_')
	for each2 in venn_tmp:
		venn_tmp2=each2.split(':')
		venn_name=groupnames[int(venn_tmp2[0])-1]+'vs'+groupnames[int(venn_tmp2[1])-1]
		display.write('%s,' % (venn_name))
	display.write('\n')
display.write('\n')
display.write('goann: %s\n' % (goann))
if species:
	display.write('KEGG species: %s\n' % (species) )
display.write('PPI number: %s\n' % (ppi_number) )
display.write('PPI blast: %s\n' % (ppi_blast))
if argv['genenamefile']:
	display.write('genenamefile: %s\n' % (genenamefile) )
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
assert not os.system('mkdir QC_TR')
os.chdir('QC_TR')
f=open(root_dir+'/QC_TR/generate_QC.sh','w')
f.write('awk \\'{if($3==\\"exon\\"){print $0}}\\' %s > %s\\n' % (gtf,root_dir+'/QC_TR/exon.gtf'))
f.write('msort -k mf1 -k nf4 %s > %s\\n' % (root_dir+'/QC_TR/exon.gtf',root_dir+'/QC_TR/sorted.gtf'))
f.write('perl /PUBLIC/source/RNA/QC/QC_v2/gtf2bed.pl %s %s\\n' % (root_dir+'/QC_TR/sorted.gtf',root_dir+'/QC_TR/sorted.bed'))
f.write('perl /PUBLIC/source/RNA/QC/QC_v2/allrunQC_v2.1.pl -fq %s -se-pe pe -n %s -o %s -spe %s -R %s -G %s -bed %s ' %(fq,sample,root_dir+'/QC_TR',ss,fa,gtf,root_dir+'/QC_TR/sorted.bed'))
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
assert not os.system('sh %s/generate_QC.sh' %(root_dir+'/QC_TR'))

samples=sample.split(',')
for eachsample in samples:
	os.chdir(root_dir+'/QC_TR/'+eachsample)
	assert not os.system('qsub -V -l vf=5g,p=8 -cwd %s' % (eachsample+'_QC.sh'))
os.chdir(root_dir)
'''
open('TR_step1_QC.py','w').write(code)
#################################################################################################
#for QCreport
code='''
import os

sample='%s'
project='%s'
root_dir='%s'
''' % (sample,project,root_dir)

code+='''

assert not os.system('mkdir QC_TR/QCreport')
os.system('sh /PUBLIC/source/RNA/QC/QC_v2/QCreport/TR_QCreport.sh -dir %s -sample %s -title %s -results %s' %(root_dir+'/QC_TR',sample,project,root_dir+'/QC_TR/QCreport'))
'''
open('TR_step1_QCreport.py','w').write(code)

#################################################################################################
#for CAN and SNP
code=''
code='''

import os

ss='%s'
sample='%s'
fa='%s'
gtf='%s'
bam='%s'
sam='%s'
group='%s'
groupname='%s'
goann='%s'
root_dir='%s'
''' % (ss,sample,fa,gtf,root_dir+'/QC_TR/bam',root_dir+'/QC_TR/sam',group,groupname,goann,root_dir)


if set([1]).issubset(includes):
	code+='''

if ss == 'no':
	lib='fr-unstranded'
elif ss == 'yes':
	lib = 'fr-firststrand'
else:
	lib = 'fr-secondstrand'

#for CAN
assert not os.system('mkdir CAN_TR')
os.chdir('CAN_TR')
f=open(root_dir+'/CAN_TR/generate_CAN.sh','w')
f.write('perl /PUBLIC/source/RNA/RefRNA/TransRef/scripts/runCAN.pl -R %s -G %s -i %s -sam %s -o %s -lib %s -group %s -groupname %s \\n' % (fa,gtf,bam,sam,root_dir+'/CAN_TR/CAN',lib,group,groupname))
f.write('sh %s/runCAN.sh\\n' % (root_dir+'/CAN_TR'))
novel_gtf=root_dir+'/CAN_TR/CAN/NovelGene/novelGene.gtf.xls'
f.write('perl /PUBLIC/source/RNA/RefRNA/TransRef/scripts/novelhmm.pl %s %s %s %s \\n' % (novel_gtf,fa,goann,root_dir+'/CAN_TR/CAN/NovelGene/GO'))
f.write('sh %s/novelhmm.sh\\n' % (root_dir+'/CAN_TR'))
f.close()
assert not os.system('qsub -V -l vf=5g,p=8 -cwd %s/generate_CAN.sh' % (root_dir+'/CAN_TR'))
os.chdir(root_dir)
'''
if set([1,2]).issubset(includes):
	code+='''
#for SNP
assert not os.system('mkdir SNP_TR')
os.chdir('SNP_TR')
f=open(root_dir+'/SNP_TR/generate_SNP.sh','w')
f.write('perl /PUBLIC/source/RNA/RefRNA/TransRef/GATKsnp/runGATK2.1.pl -R %s -t bam -i %s -o %s -b %s -n %s -gff %s\\n' % (fa,bam,root_dir+'/SNP_TR/SNP',sample,sample,gtf))
f.write('mv %s/workflow.sh %s/workflow.sh\\n' % (root_dir+'/SNP_TR/SNP',root_dir+'/SNP_TR'))
f.close()
assert not os.system('sh %s' % (root_dir+'/SNP_TR/generate_SNP.sh'))
assert not os.system('qsub -V -l vf=10g,p=6 -cwd %s/workflow.sh' % (root_dir+'/SNP_TR'))
os.chdir(root_dir)
'''
open('TR_step2_CAN_SNP.py','w').write(code)
#################################################################################################
#for DEXSeq
code=''
if set([1,3]).issubset(includes):
	code='''
import os

sam='%s'
group='%s'
groupname='%s'
compare='%s'
ss='%s'
allsamples_gtf='%s'
allsamples_tmap='%s'
root_dir='%s'
''' % (root_dir+'/QC_TR/sam',group,groupname,compare,ss,root_dir+'/CAN_TR/CAN/merged_gtf_tmap/allsamples.gtf',root_dir+'/CAN_TR/CAN/merged_gtf_tmap/allsamples.allsamples.gtf.tmap',root_dir)

	DEXseq_groups=group.split(',')
	for i,each in enumerate(DEXseq_groups):
		DEXseq_groups[i]=DEXseq_groups[i].split(':')
	DEXseq_compare=[]
	DEXseq_compares=compare.split(',')
	for each in DEXseq_compares:
		temp=each.split(':')
		tmp1=DEXseq_groups[int(temp[0])-1]
		tmp2=DEXseq_groups[int(temp[1])-1]
		tmp=tmp1+tmp2
		if (len(DEXseq_groups[int(temp[0])-1])>1 and len(DEXseq_groups[int(temp[1])-1])>1 and len(list(set(tmp)))==len(tmp1)+len(tmp2)):
			DEXseq_compare.append(each)
	DEXseq_compare=','.join(DEXseq_compare)

	code+='''
DEXseq_compare='%s'

assert not os.system('mkdir DEXseq_TR')
os.chdir('DEXseq_TR')
f=open(root_dir+'/DEXseq_TR/generate_DEXseq.sh','w')
f.write('/PUBLIC/software/public/System/Python-2.7.6/bin/python /PUBLIC/source/RNA/RefRNA/TransRef/DEXSeq/runDEXSeq_v2.py -G %%s -i %%s -g %%s -n %%s -c %%s -o %%s -p yes -s %%s -a 10 -m %%s \\n' %% (allsamples_gtf,sam,group,groupname,DEXseq_compare,root_dir+'/DEXseq_TR/DEXseq',ss,allsamples_tmap))
f.close()

assert not os.system('sh %%s' %% (root_dir+'/DEXseq_TR/generate_DEXseq.sh'))
assert not os.system('qsub -cwd -l vf=2G,p=1 -V runDEXSeq.sh')
os.chdir(root_dir)
''' % (DEXseq_compare)

open('TR_step3_DEXseq.py','w').write(code)
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
''' % (sample,ss,root_dir+'/CAN_TR/CAN/NovelGene/novelcombined.gtf',group,groupname,compare,venn,fa,number,length,root_dir)

if argv['genenamefile']:
	code+='''
genenamefile='%s'
''' % (genenamefile)

if set([1,5]).issubset(includes):
	code+='''
##Diff
sam=root_dir+'/QC_TR/sam'
samples = sample.split(',')
readcount=[]
for eachsample in samples:
	temp='%s/Diff_TR/readcount/%s.readcount' % (root_dir,eachsample)
	readcount.append(temp)
readcount=','.join(readcount)
assert not os.system('mkdir Diff_TR')
os.chdir('Diff_TR')
f=open(root_dir+'/Diff_TR/generate_Diff.sh','w')
f.write('perl /PUBLIC/source/RNA/RefRNA/DGE/DiffExpression/runDiff_analysis_v2.2.pl -fa %s -sam %s -g %s -o %s -group %s -groupname %s -compare %s -venn %s -spe %s -i %s ' % (fa,sam,gtf,root_dir+'/Diff_TR',group,groupname,compare,venn,ss,readcount))
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

if set([1,4,6]).issubset(includes):
	code+='''
##Curve
bam=root_dir+'/QC_TR/bam'
assert not os.system('mkdir Curve_TR')
os.chdir('Curve_TR')
f=open(root_dir+'/Curve_TR/generate_Curve.sh','w')
f.write('perl /PUBLIC/source/RNA/RefRNA/DGE/Curve/runCurve_v2.pl -bam %s -sample %s -n %s -r %s -fa %s -gtf %s -o %s\\n' % (bam,sample,n,length,fa,gtf,root_dir+'/Curve_TR'))
f.close()
assert not os.system('sh generate_Curve.sh')

assert not os.system('qsub -V -cwd -l vf=10G -l p=1 runSaturation.sh')
assert not os.system('qsub -V -cwd -l vf=5G -l p=1 rundensity.sh')

os.chdir(root_dir)
'''
open('TR_step4_Diff_Curves.py','w').write(code)
#################################################################################################
#for Enrichment and PPI

code='''
import os

fa='%s'
goann='%s'
groupname='%s'
compare='%s'
root_dir='%s'
''' % (fa,root_dir+'/CAN_TR/CAN/NovelGene/GO/gene.go',groupname,compare,root_dir)

code+='''
groupnames=groupname.split(',')
compares=compare.split(',')
'''
if set([1,4,5,7]).issubset(includes):
	code+='''
##GOSeq
assert not os.system('mkdir GOSeq_TR')
os.chdir('GOSeq_TR')
go=open(root_dir+'/GOSeq_TR/runGOSeq.sh','w')
for each in compares:
	temp=each.split(':')
	dir=groupnames[int(temp[0])-1]+'vs'+groupnames[int(temp[1])-1]
	id=root_dir+'/Diff_TR/Diff/'+dir+'/'+dir+'.diffgeneID'
	out=root_dir+'/GOSeq_TR/'+dir
	length=root_dir+'/Diff_TR/Diff/genelength'
	result=root_dir+'/GOSeq_TR/'+dir+'/'+dir+'.GO_enrichment_result.xls'
	diff=root_dir+'/Diff_TR/Diff/'+dir+'/'+dir+'.diffgene.xls'
	go.write('###############%s#####################\\n' % (dir))
	go.write('mkdir %s\\n' % (root_dir+'/GOSeq_TR/'+dir))
	go.write('perl /PUBLIC/source/RNA/RefRNA/DGE/Enrichment/goseq_graph_v3.pl -i %s -goann %s -n %s -o %s -length %s\\n' % (id,goann,dir,out,length))
	go.write('perl /PUBLIC/source/RNA/RefRNA/DGE/Enrichment/changeGO_up_down.pl %s %s %s\\n' % (result,diff,root_dir+'/GOSeq_TR/'+dir+'/'+dir+'.GO_enrichment_result_up_down.xls'))
	go.write('/PUBLIC/software/public/System/R-2.15.3/bin/Rscript /PUBLIC/source/RNA/RefRNA/DGE/Enrichment/goBar.R %s %s %s\\n' % (result,out,dir))
	go.write('/PUBLIC/software/public/System/R-2.15.3/bin/Rscript /PUBLIC/source/RNA/RefRNA/DGE/Enrichment/goBar2.R %s %s %s\\n' % (root_dir+'/GOSeq_TR/'+dir+'/'+dir+'.GO_enrichment_result_up_down.xls',out,dir))
go.close()
assert not os.system('qsub -V -cwd -l vf=5G -l p=1 runGOSeq.sh')
os.chdir(root_dir)
'''

if set([1,4,5]).issubset(includes) and argv['genenamefile'] == None:
	code+='''
assert not os.system('mkdir Blast_TR')
os.chdir('Blast_TR')
f=open(root_dir+'/Blast_TR/runBlast_swissprot.sh','w')
query=root_dir+'/Diff_TR/Diff/diffgene_union.seq'
outdir1=root_dir+'/Blast_TR/Blast_Swissprot/'
out=root_dir+'/Blast_TR/Blast_Swissprot/diffgene_union.seq.blastout'
f.write('echo start blastx\\ndate\\n')
f.write('mkdir %s\\n' % (outdir1))
f.write('/PUBLIC/software/public/Alignment/ncbi-blast-2.2.28+/bin/blastx -query %s -db /PUBLIC/database/Common/SwissProt/uniprot_sprot.fasta -evalue 1e-5 -outfmt 5 -max_target_seqs 1 -num_threads 10 -out %s\\n' % (query,out))
f.write('echo blastx end\\ndate\\n')
f.write('perl /PUBLIC/source/RNA/RefRNA/DGE/scriptdir/extractIDsEVxml.pl %s %s\\n' % (out,outdir1+'/diffgene_union.genenames'))
compare=root_dir+'/Diff_TR/Diff/compare.txt'
indir=root_dir+'/Diff_TR/Diff/'
outdir2=root_dir+'/Diff_TR/Diff/DiffGeneList'
f.write('mkdir %s\\n' % (outdir2))
f.write('perl /PUBLIC/source/RNA/RefRNA/DGE/scriptdir/getdiffGN.up_down_v2.pl %s %s %s %s\\n' % (indir,compare,outdir1+'/diffgene_union.genenames',outdir2))
f.close()
assert not os.system('qsub -V -cwd -l vf=5G -l p=10 runBlast_swissprot.sh')
os.chdir(root_dir)
'''
if set([1,4,5,8]).issubset(includes):
	code+='''
species='%s'
''' % (species)
	if argv['species'] == 'kaas':
		code+='''
##KOBAS
assert not os.system('mkdir KOBAS_TR')
os.chdir('KOBAS_TR')
ko=open(root_dir+'/KOBAS_TR/runKOBAS.sh','w')
path=open(root_dir+'/KOBAS_TR/runpathway.sh','w')
query=root_dir+'/Diff_TR/Diff/diffgene_union.seq'
blastout=root_dir+'/Blast_TR/KOBAS_blast.xml'
ko.write('###############Run kaas#####################\\n')
ko.write('/PUBLIC/software/public/Annotation/kaas_sa2/bin/auto_annotate.pl -n -s %s\\n' % (query))
ko.write('python /PUBLIC/software/RNA/gene_annotation/scripts/convert2kobas_v1.py %s /PUBLIC/database/Common/KEGG/kos %s\\n' % (root_dir + '/Diff_TR/Diff/diffgene_union.seq.ko',root_dir + '/KOBAS_TR/koID.annotation'))
for each in compares:
        temp=each.split(':')
        dir=groupnames[int(temp[0])-1]+'vs'+groupnames[int(temp[1])-1]
        id=root_dir+'/Diff_TR/Diff/'+dir+'/'+dir+'.diffgeneID'
        out=root_dir+'/KOBAS_TR/'+dir
        ko.write('###############%s#####################\\n' % (dir))
        path.write('echo "###############%s#####################"\\n' % (dir))
        path.write('cd %s\\n' % (out))
        result=root_dir+'/KOBAS_TR/'+dir+'/'+'add.'+dir+'.identify.xls'
        diff=root_dir+'/Diff_TR/Diff/'+dir+'/'+dir+'.diffgene.xls'
        diffID=root_dir+'/Diff_TR/Diff/'+dir+'/'+dir+'.diffgeneID'
        path.write('python /PUBLIC/source/RNA/RefRNA/DGE/KEGG_Enrichment/bin/pathway_annotation_flow_parallel_annotationfault_tolerant3.pyc --table %s --diff %s\\n' % (result,diff))
        path.write('mv %s %s\\n' % ('add.'+dir+'.identify.xls_rendered_html_detail.html',dir+'.html'))
        ko.write('mkdir %s\\n' % (root_dir+'/KOBAS_TR/'+dir))
	ko.write('cd %s\\n' % (root_dir+'/KOBAS_TR/'+dir))
        ko.write('perl /PUBLIC/source/RNA/noRef/KEGG_enrichment/runKEGG_enrich.pl -diff %s -ko %s -g %s\\n' % (diffID,root_dir + '/KOBAS_TR/koID.annotation',dir))
path.close()
ko.close()
assert not os.system('qsub -V -cwd -l vf=5G -l p=10 runKOBAS.sh')
os.chdir(root_dir)
'''
	else:
		code+='''	
##KOBAS
assert not os.system('mkdir KOBAS_TR')
os.chdir('KOBAS_TR')
ko=open(root_dir+'/KOBAS_TR/runKOBAS.sh','w')
path=open(root_dir+'/KOBAS_TR/runpathway.sh','w')
query=root_dir+'/Diff_TR/Diff/diffgene_union.seq'
blastout=root_dir+'/Blast_TR/KOBAS_blast.xml'
ko.write('###############Run Blast#####################\\n')
ko.write('mkdir %s\\n' % (root_dir+'/Blast_TR/'))
ko.write('perl /PUBLIC/source/RNA/RefRNA/DGE/KEGG_Enrichment/bin/KEGG_step1_blast.pl %s %s %s %s\\n' % (query,species,blastout,root_dir+'/Blast_TR/KOBAS_blast.sh'))
ko.write('sh %s\\n' % (root_dir+'/Blast_TR/KOBAS_blast.sh'))
for each in compares:
	temp=each.split(':')
	dir=groupnames[int(temp[0])-1]+'vs'+groupnames[int(temp[1])-1]
	id=root_dir+'/Diff_TR/Diff/'+dir+'/'+dir+'.diffgeneID'
	out=root_dir+'/KOBAS_TR/'+dir
	ko.write('###############%s#####################\\n' % (dir))
	path.write('echo "###############%s#####################"\\n' % (dir))
	path.write('cd %s\\n' % (out))
	result=root_dir+'/KOBAS_TR/'+dir+'/'+'add.'+dir+'.identify.xls'
	diff=root_dir+'/Diff_TR/Diff/'+dir+'/'+dir+'.diffgene.xls'
	path.write('python /PUBLIC/source/RNA/RefRNA/DGE/KEGG_Enrichment/bin/pathway_annotation_flow_parallel_annotationfault_tolerant3.pyc --table %s --diff %s\\n' % (result,diff))
	path.write('mv %s %s\\n' % ('add.'+dir+'.identify.xls_rendered_html_detail.html',dir+'.html'))
	ko.write('mkdir %s\\n' % (root_dir+'/KOBAS_TR/'+dir))
	script=root_dir+'/KOBAS_TR/'+dir+'/run.sh'
	ko.write('perl /PUBLIC/source/RNA/RefRNA/DGE/KEGG_Enrichment/bin/KEGG_step2_enrich.pl -id %s -out-dir %s -species %s -blast-result %s -sample-names %s>%s\\n' % (id,out,species,blastout,dir,script))
	ko.write('sh %s\\n' % (script))
path.close()
ko.close()
assert not os.system('qsub -V -cwd -l vf=5G -l p=10 runKOBAS.sh')
os.chdir(root_dir)
'''

if set([1,4,5,9]).issubset(includes):
	code+='''
ppi_number='%s'
ppi_blast='%s'
''' % (ppi_number,ppi_blast)
	code+='''
##PPI
assert not os.system('mkdir PPI_TR')
os.chdir('PPI_TR')
ppi=open(root_dir+'/PPI_TR/runPPI.sh','w')
for each in compares:
	temp=each.split(':')
	dir=groupnames[int(temp[0])-1]+'vs'+groupnames[int(temp[1])-1]
	id=root_dir+'/Diff_TR/Diff/'+dir+'/'+dir+'.diffgeneID'
	seq=root_dir+'/Diff_TR/Diff/Diff_Gene_Seq/'+dir+'.diffgene.seq'
	out=root_dir+'/PPI_TR/'+dir
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
open('TR_step5_Enrichment.py','w').write(code)

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
settings.configure(DEBUG=True, TEMPLATE_DEBUG=True,TEMPLATE_DIRS=('/PUBLIC/source/RNA/RefRNA/Report/template/TR',))

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
dir=root_dir+'/QC_TR'
cd_md5=open(root_dir+'/QC_TR/cd_md5.sh','w')
rd_md5=open(root_dir+'/QC_TR/rd_md5.sh','w')
os.chdir(dir)
cd_md5.write('cd %s/cleandata\\n' % (dir))
cd_md5.write('echo -e \"md5\\tfile\">cd_md5.txt\\n')
rd_md5.write('echo -e \"md5\\tfile\">rd_md5.txt\\n')
for each in samples:
	cd_md5.write('md5sum %s_1.clean.fq>>cd_md5.txt\\n' % (each))
	cd_md5.write('md5sum %s_2.clean.fq>>cd_md5.txt\\n' % (each))
	cd_md5.write('unlink %s_1.clean.fq\\n' %(each))
	cd_md5.write('unlink %s_2.clean.fq\\n' %(each))
	cd_temp1=root_dir+'/QC_TR/'+each+'/clean_data/'+each+'_1.clean.fq'
	cd_temp2=root_dir+'/QC_TR/'+each+'/clean_data/'+each+'_2.clean.fq'
	cd_md5.write('gzip %s\\n' % (cd_temp1))
	cd_md5.write('gzip %s\\n' % (cd_temp2))
	cd_md5.write('ln -s %s.gz %s\\n' % (cd_temp1,dir+'/cleandata/'))
	cd_md5.write('ln -s %s.gz %s\\n' % (cd_temp2,dir+'/cleandata/'))
	rd_md5.write('cd %s\\n' % (root_dir+'/QC_TR/'+each))
	rd_md5.write('md5sum %s>>%srd_md5.txt\\n' % (each+'_1.fq.gz',root_dir+'/QC_TR/'))
	rd_md5.write('md5sum %s>>%srd_md5.txt\\n' % (each+'_2.fq.gz',root_dir+'/QC_TR/'))
cd_md5.close()
rd_md5.close()
os.chdir(root_dir)

os.chdir(root_dir+'/QC_TR')
rdmd5=glob.glob('rd_md5.txt')
if not rdmd5:
        assert not os.system('qsub -V -l vf=3g,p=1 -cwd rd_md5.sh')
os.chdir(root_dir+'/QC_TR/cleandata')
cdmd5=glob.glob('cd_md5.txt')
if not cdmd5:
        os.chdir(root_dir+'//QC_TR')
        assert not os.system('qsub -V -l vf=3g,p=1 -cwd cd_md5.sh')
os.chdir(root_dir)
'''
if set([1,4,5,8]).issubset(includes):
	code+='''
##path
path=root_dir+'/KOBAS_TR/runpathway.sh'
assert not os.system('sh %s' % (path))
'''

code+='''
#results

assert not os.system('mkdir '+project+'_TR_result')
assert not os.system('mkdir '+project+'_TR_result/'+project+'_results')
assert not os.system('mkdir '+project+'_TR_result/'+project+'_report')
os.chdir(project+'_TR_result/'+project+'_results')

result_order=0

result_order+=1
assert not os.system('mkdir '+str(result_order)+'.OriginalData')
rd_order=str(result_order)
os.chdir(''+str(result_order)+'.OriginalData')
for eachsample in samples:
	f_in=gzip.open(root_dir+'/QC_TR/'+eachsample+'/'+eachsample+'_1.fq.gz','rb')
	f_out=[f_in.readline() for num in range(20)]
	open(eachsample+'.example.fq.txt','w').writelines(f_out)
	f_in.close()
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/OriginalData.README.txt %s' % (root_dir+'/'+project+'_TR_result/'+project+'_results/'+str(result_order)+'.OriginalData/'+str(result_order)+'.OriginalData.README.txt'))
os.chdir(root_dir+'/'+project+'_TR_result/'+project+'_results')

result_order+=1
assert not os.system('mkdir '+str(result_order)+'.QC')
qc_order=str(result_order)
os.chdir(''+str(result_order)+'.QC')
assert not os.system('mkdir '+str(result_order)+'.1.ErrorRate')
os.chdir(''+str(result_order)+'.1.ErrorRate')
for eachsample in samples:
	assert not os.system('cp %s %s' % (root_dir+'/QC_TR/'+eachsample+'/clean_data/'+eachsample+'.Error.png',eachsample+'.error_rate_distribution.png'))
	assert not os.system('cp %s %s' % (root_dir+'/QC_TR/'+eachsample+'/clean_data/'+eachsample+'.Error.pdf',eachsample+'.error_rate_distribution.pdf'))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/Error.README.txt %s' % (root_dir+'/'+project+'_TR_result/'+project+'_results/'+str(result_order)+'.QC/'+str(result_order)+'.1.ErrorRate/'+str(result_order)+'.1ErrorRate.README.txt'))
os.chdir('..')

assert not os.system('mkdir '+str(result_order)+'.2.GC')
os.chdir(''+str(result_order)+'.2.GC')
for eachsample in samples:
	assert not os.system('cp %s %s' % (root_dir+'/QC_TR/'+eachsample+'/clean_data/'+eachsample+'.GC.png',eachsample+'.GC_content_distribution.png'))
	assert not os.system('cp %s %s' % (root_dir+'/QC_TR/'+eachsample+'/clean_data/'+eachsample+'.GC.pdf',eachsample+'.GC_content_distribution.pdf'))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/GC.README.txt %s' % (root_dir+'/'+project+'_TR_result/'+project+'_results/'+str(result_order)+'.QC/'+str(result_order)+'.2.GC/'+str(result_order)+'.2GC.README.txt'))
os.chdir('..')

assert not os.system('mkdir '+str(result_order)+'.3.ReadsClassification')
os.chdir(''+str(result_order)+'.3.ReadsClassification')
for eachsample in samples:
	assert not os.system('cp %s %s' % (root_dir+'/QC_TR/'+eachsample+'/files/'+eachsample+'.pie3d.png',eachsample+'.raw_reads_classification.png'))
	assert not os.system('cp %s %s' % (root_dir+'/QC_TR/'+eachsample+'/files/'+eachsample+'.pie3d.pdf',eachsample+'.raw_reads_classification.pdf'))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/Filter.README.txt %s' % (root_dir+'/'+project+'_TR_result/'+project+'_results/'+str(result_order)+'.QC/'+str(result_order)+'.3.ReadsClassification/'+str(result_order)+'.3ReadsClassification.README.txt'))
os.chdir('..')

assert not os.system('mkdir '+str(result_order)+'.4.DataTable')
os.chdir(''+str(result_order)+'.4.DataTable')
alldataTable = [ root_dir + '/QC_TR/' + sample + '/files/dataTable' for sample in samples]
alldataTable = ' '.join(alldataTable)
assert not os.system("cat %s|cut -f 1-8 > %s " % ( alldataTable ,'datatable.xls.1'))	
assert not os.system('cat %s %s > %s' % ('/PUBLIC/source/RNA/RefRNA/DGE/Report/datatable','datatable.xls.1','datatable.xls'))
assert not os.system('rm %s' % ('datatable.xls.1'))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/DataTable.README.txt %s' % (root_dir+'/'+project+'_TR_result/'+project+'_results/'+str(result_order)+'.QC/'+str(result_order)+'.4.DataTable/'+str(result_order)+'.4DataTable.README.txt'))
os.chdir(root_dir+'/'+project+'_TR_result/'+project+'_results')

result_order+=1
assert not os.system('mkdir '+str(result_order)+'.Mapping')
mp_order=str(result_order)
os.chdir(''+str(result_order)+'.Mapping')
assert not os.system('mkdir '+str(result_order)+'.1.MapStat')
os.chdir(''+str(result_order)+'.1.MapStat')
for eachsample in samples:
	assert not os.system('cut -f 2 %s > %s' % (root_dir+'/QC_TR/'+eachsample+'/files/'+eachsample+'.stat',eachsample+'.stat.2'))
assert not os.system('cut -f 1 %s > %s' % (root_dir+'/QC_TR/'+samples[0]+'/files/'+samples[0]+'.stat','stat.1'))
allStat_2 = [ root_dir + '/' + project+'_TR_result/'+project+'_results/'+str(result_order)+'.Mapping/'+str(result_order)+'.1.MapStat/'+sample+'.stat.2'  for sample in samples ]
allStat_2 = ' '.join(allStat_2)
assert not os.system('paste %s %s > %s' % ('stat.1',allStat_2,'MapStat.xls'))
assert not os.system("sed -i -e '/^Reads mapped in proper pairs/d' %s" % ('MapStat.xls'))
assert not os.system("sed -i -e '/^Proper-paired reads map to different chrom/d' %s" % ('MapStat.xls'))
assert not os.system('rm stat.1')
assert not os.system('rm *.stat.2')
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/MapStat.README.txt %s' % (root_dir+'/'+project+'_TR_result/'+project+'_results/'+str(result_order)+'.Mapping/'+str(result_order)+'.1.MapStat/'+str(result_order)+'.1MapStat.README.txt'))
os.chdir('..')

assert not os.system('mkdir '+str(result_order)+'.2.MapReg')
mr_order=str(result_order)
os.chdir(''+str(result_order)+'.2.MapReg')
for eachsample in samples:
	assert not os.system('cp %s %s' % (root_dir+'/QC_TR/'+eachsample+'/files/'+eachsample+'.MR.png',eachsample+'.Mapped_Region.png'))
	assert not os.system('cp %s %s' % (root_dir+'/QC_TR/'+eachsample+'/files/'+eachsample+'.MR.pdf',eachsample+'.Mapped_Region.pdf'))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/MapReg.README.txt %s' % (root_dir+'/'+project+'_TR_result/'+project+'_results/'+str(result_order)+'.Mapping/'+str(result_order)+'.2.MapReg/'+str(result_order)+'.2MappedRegion.README.txt'))
os.chdir('..')

assert not os.system('mkdir '+str(result_order)+'.3.ChrDen')
cd_order=str(result_order)
os.chdir(''+str(result_order)+'.3.ChrDen')
for eachsample in samples:
	assert not os.system('cp %s .' % (root_dir+'/Curve_TR/density/'+eachsample+'/'+eachsample+'.density.png'))
	assert not os.system('cp %s .' % (root_dir+'/Curve_TR/density/'+eachsample+'/'+eachsample+'.density.pdf'))
	assert not os.system('cp %s .' % (root_dir+'/Curve_TR/density/'+eachsample+'/'+eachsample+'.ReadsVSlength.png'))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/ChrDen.README.txt %s' % (root_dir+'/'+project+'_TR_result/'+project+'_results/'+str(result_order)+'.Mapping/'+str(result_order)+'.3.ChrDen/'+str(result_order)+'.3ChrDen.README.txt'))
os.chdir('..')

assert not os.system('mkdir '+str(result_order)+'.4.IGV')
igv_order=str(result_order)
os.chdir(''+str(result_order)+'.4.IGV')
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/Report/IGVQuickStart.pdf %s' % (root_dir+'/'+project+'_TR_result/'+project+'_results/'+str(result_order)+'.Mapping/'+str(result_order)+'.4.IGV'))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/IGV.README.txt %s' % (root_dir+'/'+project+'_TR_result/'+project+'_results/'+str(result_order)+'.Mapping/'+str(result_order)+'.4.IGV/'+str(result_order)+'.4IGV.README.txt'))
os.chdir(root_dir+'/'+project+'_TR_result/'+project+'_results')
'''

if set([1]).issubset(includes):
	code+='''
result_order+=1
assert not os.system('mkdir '+str(result_order)+'.AS')
as_order=str(result_order)
os.chdir(''+str(result_order)+'.AS')
assert not os.system('cp %s .' % (root_dir+'/CAN_TR/CAN/ASprofile/AS.png'))
assert not os.system('cp %s .' % (root_dir+'/CAN_TR/CAN/ASprofile/AS.pdf'))
for eachsample in samples:
	assert not os.system('cp %s .' % (root_dir+'/CAN_TR/CAN/ASprofile/'+eachsample+'/'+eachsample+'.fpkm.xls'))
assert not os.system('cp %s .' % (root_dir+'/CAN_TR/CAN/ASprofile/ASevent.stat.xls'))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/AS.README.txt %s' % (root_dir+'/'+project+'_TR_result/'+project+'_results/'+str(result_order)+'.AS/'+str(result_order)+'.AS.README.txt'))
os.chdir(root_dir+'/'+project+'_TR_result/'+project+'_results')

result_order+=1
assert not os.system('mkdir '+str(result_order)+'.NovelGene')
ng_order=str(result_order)
os.chdir(''+str(result_order)+'.NovelGene')
assert not os.system('cp %s .' % (root_dir+'/CAN_TR/CAN/NovelGene/novelGene.gtf.xls'))
assert not os.system('cp %s .' % (root_dir+'/CAN_TR/CAN/NovelGene/novelIsoform.gtf.xls'))
assert not os.system('cp %s .' % (root_dir+'/CAN_TR/CAN/NovelGene/geneStructOpt.xls'))
assert not os.system('cp %s .' % (root_dir+'/CAN_TR/CAN/NovelGene/GO/Novelgene.fa'))
assert not os.system('cp %s .' % (root_dir+'/CAN_TR/CAN/NovelGene/GO/Novelgene.pfam.txt'))
assert not os.system('cp %s .' % (root_dir+'/CAN_TR/CAN/NovelGene/GO/Novelgene.hmm_go.txt'))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/NovelGene.README.txt %s' % (root_dir+'/'+project+'_TR_result/'+project+'_results/'+str(result_order)+'.NovelGene/'+str(result_order)+'.NovelGene.README.txt'))
os.chdir(root_dir+'/'+project+'_TR_result/'+project+'_results')
'''

if set([1,2]).issubset(includes):
	code+='''
result_order+=1
assert not os.system('mkdir '+str(result_order)+'.SNP')
snp_order=str(result_order)
os.chdir(''+str(result_order)+'.SNP')
assert not os.system('cp %s .' % (root_dir+'/SNP_TR/SNP/ResultsQ30/InDels.xls'))
assert not os.system('cp %s .' % (root_dir+'/SNP_TR/SNP/ResultsQ30/SNPs.xls'))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/SNP.README.txt %s' % (root_dir+'/'+project+'_TR_result/'+project+'_results/'+str(result_order)+'.SNP/'+str(result_order)+'.SNP.README.txt'))
os.chdir(root_dir+'/'+project+'_TR_result/'+project+'_results')
'''

if set([1,4,5]).issubset(includes):
	code+='''
result_order+=1
assert not os.system('mkdir '+str(result_order)+'.GeneExprQuatification')
HTSeq_order=str(result_order)
os.chdir(''+str(result_order)+'.GeneExprQuatification')
assert not os.system('cp %s .' % (root_dir+'/Diff_TR/Diff/rpkm.xls'))
assert not os.system('cp %s .' % (root_dir+'/Diff_TR/Diff/readcount.xls'))
assert not os.system('cp %s .' % (root_dir+'/Diff_TR/Diff/rpkm.stat.xls'))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/Quatification.README.txt %s' % (root_dir+'/'+project+'_TR_result/'+project+'_results/'+str(result_order)+'.GeneExprQuatification/'+str(result_order)+'GeneExprQuatification.README.txt'))
os.chdir(root_dir+'/'+project+'_TR_result/'+project+'_results')
'''

if set([1,4,6]).issubset(includes):
	code+='''
result_order+=1
assert not os.system('mkdir '+str(result_order)+'.AdvancedQC')
os.chdir(''+str(result_order)+'.AdvancedQC')
assert not os.system('mkdir '+str(result_order)+'.1.SaturationCurve')
os.chdir(''+str(result_order)+'.1.SaturationCurve')
for eachsample in samples:
	assert not os.system('cp %s .' % (root_dir+'/Curve_TR/SaturationCurve/'+eachsample+'/'+eachsample+'.saturation.png'))
	assert not os.system('cp %s .' % (root_dir+'/Curve_TR/SaturationCurve/'+eachsample+'/'+eachsample+'.saturation.pdf'))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/SatCurve.README.txt %s' % (root_dir+'/'+project+'_TR_result/'+project+'_results/'+str(result_order)+'.AdvancedQC/'+str(result_order)+'.1.SaturationCurve/'+str(result_order)+'.1SaturationCurve.README.txt'))
os.chdir('..')

if flag_uniform:
	assert not os.system('mkdir '+str(result_order)+'.2.MeanCoverage')
	os.chdir(''+str(result_order)+'.2.MeanCoverage')
	for eachsample in samples:
		assert not os.system('cp %s %s' % (root_dir+'/QC_TR/'+eachsample+'/files/'+eachsample+'.MC.png',eachsample+'.Mean_coverage_distribution.png'))
		assert not os.system('cp %s %s' % (root_dir+'/QC_TR/'+eachsample+'/files/'+eachsample+'.MC.pdf',eachsample+'.Mean_coverage_distribution.pdf'))
	assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/MeanCov.README.txt %s' % (root_dir+'/'+project+'_TR_result/'+project+'_results/'+str(result_order)+'.AdvancedQC/'+str(result_order)+'.2.MeanCoverage/'+str(result_order)+'.2MeanCoverage.README.txt'))
os.chdir(root_dir+'/'+project+'_TR_result/'+project+'_results')
'''
if set([1,4,5]).issubset(includes):
	code +='''
if flag_uniform:
	os.chdir(''+str(result_order)+'.AdvancedQC')
	assert not os.system('mkdir '+str(result_order)+'.3.Correlation')
	os.chdir(''+str(result_order)+'.3.Correlation')
	for png in glob.iglob(root_dir+'/Diff_TR/Diff/corr_plot/*.png'):
		assert not os.system('cp %s .' % png)
	for pdf in glob.iglob(root_dir+'/Diff_TR/Diff/corr_plot/*.pdf'):
		assert not os.system('cp %s .' % pdf)
	assert not os.system( 'cp %s .' % (root_dir+'/Diff_TR/Diff/corr_plot/cor_pearson.xls') )
	assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/DupCorr.README.txt %s' % (root_dir+'/'+project+'_TR_result/'+project+'_results/'+str(result_order)+'.AdvancedQC/'+str(result_order)+'.3.Correlation/'+str(result_order)+'.3Correlation.README.txt'))
	os.chdir(root_dir+'/'+project+'_TR_result/'+project+'_results')
else:
	os.chdir(''+str(result_order)+'.AdvancedQC')
	assert not os.system('mkdir '+str(result_order)+'.2.Correlation')
	os.chdir(''+str(result_order)+'.2.Correlation')
	for png in glob.iglob(root_dir+'/Diff_TR/Diff/corr_plot/*.png'):
		assert not os.system('cp %s .' % png)
	for pdf in glob.iglob(root_dir+'/Diff_TR/Diff/corr_plot/*.pdf'):
		assert not os.system('cp %s .' % pdf)
	assert not os.system( 'cp %s .' % (root_dir+'/Diff_TR/Diff/corr_plot/cor_pearson.xls') )
	assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/DupCorr.README.txt %s' % (root_dir+'/'+project+'_TR_result/'+project+'_results/'+str(result_order)+'.AdvancedQC/'+str(result_order)+'.2.Correlation/'+str(result_order)+'.2Correlation.README.txt'))
	os.chdir(root_dir+'/'+project+'_TR_result/'+project+'_results')
'''

if set([1,4,5]).issubset(includes):
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
assert not os.system('cp %s .' % (root_dir+'/Diff_TR/Diff/boxplot.pdf'))
assert not os.system('cp %s .' % (root_dir+'/Diff_TR/Diff/boxplot.png'))
assert not os.system('cp %s .' % (root_dir+'/Diff_TR/Diff/density.pdf'))
assert not os.system('cp %s .' % (root_dir+'/Diff_TR/Diff/density.png'))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/ExpLev.README.txt %s' % (root_dir+'/'+project+'_TR_result/'+project+'_results/'+str(result_order)+'.DiffExprAnalysis/'+str(result_order)+'.1.GeneExpContrast/'+str(result_order)+'.1GeneExpContrast.README.txt'))
os.chdir('..')

assert not os.system('mkdir '+str(result_order)+'.2.DEGsList')
os.chdir(''+str(result_order)+'.2.DEGsList')
for each in compare_names:
	assert not os.system('cp %s %s' % (root_dir+'/Diff_TR/Diff/DiffGeneList/'+each+'.diffgene.genename.xls',each+'.DEG.xls'))
	assert not os.system('cp %s %s' % (root_dir+'/Diff_TR/Diff/DiffGeneList/'+each+'.diffgene_up.genename.xls',each+'.DEG_up.xls'))
	assert not os.system('cp %s %s' % (root_dir+'/Diff_TR/Diff/DiffGeneList/'+each+'.diffgene_down.genename.xls',each+'.DEG_down.xls'))
	assert not os.system('cp %s %s' % (root_dir+'/Diff_TR/Diff/Diff_Gene_Seq/'+each+'.diffgene.seq',each+'.DEG.fasta'))
	assert not os.system('cp %s .' % (root_dir+'/Diff_TR/Diff/'+each+'/'+each+'.Differential_analysis_results.xls'))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/DiffList.README.txt %s' % (root_dir+'/'+project+'_TR_result/'+project+'_results/'+str(result_order)+'.DiffExprAnalysis/'+str(result_order)+'.2.DEGsList/'+str(result_order)+'.2DEGsList.README.txt'))
os.chdir('..')

assert not os.system('mkdir '+str(result_order)+'.3.DEGsFilter')
os.chdir(''+str(result_order)+'.3.DEGsFilter')
for each in compare_names:
	assert not os.system('cp %s %s' % (root_dir+'/Diff_TR/Diff/'+each+'/'+each+'.volcano.pdf',each+'.Volcanoplot.pdf'))
	assert not os.system('cp %s %s' % (root_dir+'/Diff_TR/Diff/'+each+'/'+each+'.volcano.png',each+'.Volcanoplot.png'))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/DiffFilt.README.txt %s' % (root_dir+'/'+project+'_TR_result/'+project+'_results/'+str(result_order)+'.DiffExprAnalysis/'+str(result_order)+'.3.DEGsFilter/'+str(result_order)+'.3DEGsFilter.README.txt'))
os.chdir('..')

assert not os.system('mkdir '+str(result_order)+'.4.DEGcluster')
os.chdir(''+str(result_order)+'.4.DEGcluster')
assert not os.system('cp %s Union_for_cluster.xls' % ( root_dir+'/Diff_TR/Diff/union_for_cluster' ))
assert not os.system('cp %s .' % (root_dir+'/Diff_TR/Diff/heatCluster.pdf'))
os.system('cp %s .' % (root_dir+'/Diff_TR/Diff/heatCluster.png'))
assert not os.system('cp %s .' % (root_dir+'/Diff_TR/Diff/heatCluster.detail.pdf'))
assert not os.system('cp %s .' % (root_dir+'/Diff_TR/Diff/union_for_cluster_hclust/h_show_plots.png'))
assert not os.system('cp %s .' % (root_dir+'/Diff_TR/Diff/union_for_cluster_kmeans/kmeans_cluster_plots.pdf'))
assert not os.system('cp %s .' % (root_dir+'/Diff_TR/Diff/union_for_cluster_som/som_cluster_plots.pdf'))
assert not os.system('cp %s %s' % (root_dir+'/Diff_TR/Diff/union_for_cluster_hclust/h_subcluster_plots.pdf','h_cluster_plots.pdf'))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/Cluster.README.txt %s' % (root_dir+'/'+project+'_TR_result/'+project+'_results/'+str(result_order)+'.DiffExprAnalysis/'+str(result_order)+'.4.DEGcluster/'+str(result_order)+'.4DEGcluster.README.txt'))
assert not os.system('mkdir '+'Subcluster')
assert not os.system('mkdir '+'Subcluster/kmeans')
assert not os.system('mkdir '+'Subcluster/hclust')
assert not os.system('mkdir '+'Subcluster/som')
assert not os.system('ls %s|while read a;do cut -f 1 $a>%s${a##*/}.xls;done' % (root_dir+'/Diff_TR/Diff/union_for_cluster_kmeans/subcluster*','Subcluster/kmeans/'))
assert not os.system('ls %s|while read a;do cut -f 1 $a>%s${a##*/}.xls;done' % (root_dir+'/Diff_TR/Diff/union_for_cluster_hclust/subcluster*','Subcluster/hclust/'))
assert not os.system('ls %s|while read a;do cut -f 1 $a>%s${a##*/}.xls;done' % (root_dir+'/Diff_TR/Diff/union_for_cluster_som/subcluster*','Subcluster/som/'))
os.chdir('..')

if flag_venn:
	assert not os.system('mkdir '+str(result_order)+'.5.VennDiagram')
	os.chdir(''+str(result_order)+'.5.VennDiagram')
	for each in venn_cluster_vs_names:
		assert not os.system('cp -r %s .' % (root_dir+'/Diff_TR/Diff/Venn/'+each))
	assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/Venn.README.txt %s' % (root_dir+'/'+project+'_TR_result/'+project+'_results/'+str(result_order)+'.DiffExprAnalysis/'+str(result_order)+'.5.VennDiagram/'+str(result_order)+'.5VennDiagram.README.txt'))
os.chdir(root_dir+'/'+project+'_TR_result/'+project+'_results')
'''

if set([1,4,5,7]).issubset(includes):
	code+='''
result_order+=1
assert not os.system('mkdir '+str(result_order)+'.DEG_GOEnrichment')
os.chdir(''+str(result_order)+'.DEG_GOEnrichment')
assert not os.system('mkdir '+str(result_order)+'.1.DEG_GOList')
os.chdir(''+str(result_order)+'.1.DEG_GOList')
for each in compare_names:
	assert not os.system('cp %s .' % (root_dir+'/GOSeq_TR/'+each+'/'+each+'.GO_enrichment_result_up_down.xls'))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/GOList.README.txt %s' % (root_dir+'/'+project+'_TR_result/'+project+'_results/'+str(result_order)+'.DEG_GOEnrichment/'+str(result_order)+'.1.DEG_GOList/'+str(result_order)+'.1DEG_GOList.README.txt'))
os.chdir('..')

assert not os.system('mkdir '+str(result_order)+'.2.DAG')
os.chdir(''+str(result_order)+'.2.DAG')
for each in compare_names:
	os.system('cp %s .' % (root_dir+'/GOSeq_TR/'+each+'/'+each+'.GObp_DAG.png'))
	os.system('cp %s .' % (root_dir+'/GOSeq_TR/'+each+'/'+each+'.GObp_DAG.pdf'))
	os.system('cp %s .' % (root_dir+'/GOSeq_TR/'+each+'/'+each+'.GOmf_DAG.png'))
	os.system('cp %s .' % (root_dir+'/GOSeq_TR/'+each+'/'+each+'.GOmf_DAG.pdf'))
	os.system('cp %s .' % (root_dir+'/GOSeq_TR/'+each+'/'+each+'.GOcc_DAG.png'))
	os.system('cp %s .' % (root_dir+'/GOSeq_TR/'+each+'/'+each+'.GOcc_DAG.pdf'))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/DAG.README.txt %s' % (root_dir+'/'+project+'_TR_result/'+project+'_results/'+str(result_order)+'.DEG_GOEnrichment/'+str(result_order)+'.2.DAG/'+str(result_order)+'.2DAG.README.txt'))
os.chdir('..')

assert not os.system('mkdir '+str(result_order)+'.3.BAR')
os.chdir(''+str(result_order)+'.3.BAR')
for each in compare_names:
	assert not os.system('cp %s .' % (root_dir+'/GOSeq_TR/'+each+'/'+each+'.bar_graph.png'))
	assert not os.system('cp %s .' % (root_dir+'/GOSeq_TR/'+each+'/'+each+'.bar_graph.pdf'))
	assert not os.system('cp %s .' % (root_dir+'/GOSeq_TR/'+each+'/'+each+'.bar_graph_updown.png'))
	assert not os.system('cp %s .' % (root_dir+'/GOSeq_TR/'+each+'/'+each+'.bar_graph_updown.pdf'))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/BAR.README.txt %s' % (root_dir+'/'+project+'_TR_result/'+project+'_results/'+str(result_order)+'.DEG_GOEnrichment/'+str(result_order)+'.3.BAR/'+str(result_order)+'.3BAR.README.txt'))
os.chdir(root_dir+'/'+project+'_TR_result/'+project+'_results')
'''

if set([1,4,5,8]).issubset(includes):
	code+='''
result_order+=1
assert not os.system('mkdir '+str(result_order)+'.DEG_KEGGEnrichment')
os.chdir(''+str(result_order)+'.DEG_KEGGEnrichment')
assert not os.system('mkdir '+str(result_order)+'.1.DEG_KEGGList')
os.chdir(''+str(result_order)+'.1.DEG_KEGGList')
for each in compare_names:
	assert not os.system('cp %s %s' % (root_dir+'/KOBAS_TR/'+each+'/add.'+each+'.identify.xls',each+'.DEG_KEGG_pathway_enrichment_result.xls'))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/KEGGList.README.txt %s' % (root_dir+'/'+project+'_TR_result/'+project+'_results/'+str(result_order)+'.DEG_KEGGEnrichment/'+str(result_order)+'.1.DEG_KEGGList/'+str(result_order)+'.1DEG_KEGGList.README.txt'))
os.chdir('..')

assert not os.system('mkdir '+str(result_order)+'.2.DEG_KEGGScat')
os.chdir(''+str(result_order)+'.2.DEG_KEGGScat')
for each in compare_names:
	assert not os.system('cp %s .' % (root_dir+'/KOBAS_TR/'+each+'/'+each+'.DEG_enriched_KEGG_pathway_scatterplot.png'))
	assert not os.system('cp %s .' % (root_dir+'/KOBAS_TR/'+each+'/'+each+'.DEG_enriched_KEGG_pathway_scatterplot.pdf'))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/KEGGScat.README.txt %s' % (root_dir+'/'+project+'_TR_result/'+project+'_results/'+str(result_order)+'.DEG_KEGGEnrichment/'+str(result_order)+'.2.DEG_KEGGScat/'+str(result_order)+'.2DEG_KEGGScat.README.txt'))
os.chdir('..')

assert not os.system('mkdir '+str(result_order)+'.3.DEG_KEGGPath')
os.chdir(''+str(result_order)+'.3.DEG_KEGGPath')
for each in compare_names:
	assert not os.system('mkdir '+each)
	assert not os.system('cp %s %s' % (root_dir+'/KOBAS_TR/'+each+'/'+each+'.html',each))
	assert not os.system('cp -r %s %s' % (root_dir+'/KOBAS_TR/'+each+'/src',each))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/KEGGPath.README.txt %s' % (root_dir+'/'+project+'_TR_result/'+project+'_results/'+str(result_order)+'.DEG_KEGGEnrichment/'+str(result_order)+'.3.DEG_KEGGPath/'+str(result_order)+'.3DEG_KEGGPath.README.txt'))
os.chdir(root_dir+'/'+project+'_TR_result/'+project+'_results')
'''

if set([1,2,4,5,9]).issubset(includes):
	code+='''
result_order+=1
assert not os.system('mkdir '+str(result_order)+'.DEG_PPI')
os.chdir(''+str(result_order)+'.DEG_PPI')
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/Report/CytoscapeQuickStart.pdf .')
for each in compare_names:
	assert not os.system('cp %s .' % (root_dir+'/PPI_TR/'+each+'/'+each+'.ppi.txt'))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/ProteinNetwork.README.txt %s' % (root_dir+'/'+project+'_TR_result/'+project+'_results/'+str(result_order)+'.DEG_PPI/'+str(result_order)+'.ProteinNetwork.README.txt'))
os.chdir(root_dir+'/'+project+'_TR_result/'+project+'_results')
'''

if set([1,3]).issubset(includes):
	code+='''
DEXseq_compare = '%s'
''' % (DEXseq_compare)

	code+='''
DEXseq_compare_names=[]
if DEXseq_compare !='':
        DEXseq_compares=DEXseq_compare.split(',')
        for each in DEXseq_compares:
                tmp=each.split(':')
                DEXseq_compare_names.append(groupnames[int(tmp[0])-1]+'vs'+groupnames[int(tmp[1])-1])
#        DEXseq_compare_names=','.join(DEXseq_compare_names)
#DEXseq_compare='%s'
#if DEXseq_compare !='':
#        DEXseq_compares=DEXseq_compare.split(',')
#        for each in DEXseq_compares:
#                tmp=each.split(':')
#                DEXseq_compare_names.append(groupnames[int(tmp[0])-1]+'vs'+groupnames[int(tmp[1])-1])
#        DEXseq_compare_names=','.join(DEXseq_compare_names)

result_order+=1
assert not os.system('mkdir '+str(result_order)+'.DEU')
os.chdir(''+str(result_order)+'.DEU')
for each in DEXseq_compare_names:
	assert not os.system('cp -R %s .' % (root_dir+'/DEXseq_TR/DEXseq/'+each+'/'+each+'.DEU'))
	assert not os.system('cp %s .' % (root_dir+'/DEXseq_TR/DEXseq/'+each+'/'+each+'.DEUdiff.xls'))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/DEU.README.txt %s' % (root_dir+'/'+project+'_TR_result/'+project+'_results/'+str(result_order)+'.DEU/'+str(result_order)+'.DEU.README.txt'))
os.chdir(root_dir+'/'+project+'_TR_result/'+project+'_results')
'''

if ss == 'yes' or ss == 'reverse':
	code+='''
result_order+=1
assert not os.system('mkdir '+str(result_order)+'.AntiTrans')
os.chdir(''+str(result_order)+'.AntiTrans')
for eachsample in samples:
	assert not os.system('cp %s .' % (root_dir+'/CAN_TR/CAN/stranded_specific/cisNAT/'+eachsample+'.nat.xls'))
	assert not os.system('cp %s .' % (root_dir+'/CAN_TR/CAN/stranded_specific/cisNAT/'+eachsample+'.nat.stat.xls'))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/AntiTrans.README.txt %s' % (root_dir+'/'+project+'_TR_result/'+project+'_results/'+str(result_order)+'.AntiTrans/'+str(result_order)+'.AntiTrans.README.txt'))
os.chdir(root_dir+'/'+project+'_TR_result/'+project+'_results')

result_order+=1
assert not os.system('mkdir '+str(result_order)+'.TSS_TTS')
os.chdir(''+str(result_order)+'.TSS_TTS')
assert not os.system('cp %s %s' % (root_dir+'/CAN_TR/CAN/stranded_specific/TSS.xls','TSS_TTS.xls'))
assert not os.system('cp /PUBLIC/source/RNA/RefRNA/TransRef/README_v2/TSS_TTS.README.txt %s' % (root_dir+'/'+project+'_TR_result/'+project+'_results/'+str(result_order)+'.TSS_TTS/'+str(result_order)+'.TSS_TTS.README.txt'))
os.chdir(root_dir+'/'+project+'_TR_result/'+project+'_results')
'''

code+='''

print project+' report generating...'
os.chdir(root_dir)
os.chdir(project+'_TR_result/'+project+'_results')
os.system('tree -d -v -L 3 -o ../'+project+'_report/DirectoryTree.html -H ../'+project+'_results')
# add chinese desciption for DirectoryTree
DirectoryTree = root_dir+'/'+project+'_TR_result/'+project+'_report/DirectoryTree.html'
NewDirectoryTree = DirectoryTree + '.new'
os.system( 'perl /PUBLIC/source/RNA/RefRNA/Report/bin/result_tree.pl %s %s' % (DirectoryTree, NewDirectoryTree) )
os.system('mv %s %s' % (NewDirectoryTree, DirectoryTree))
os.chdir(root_dir)
#generate the TransRef html report (separate the processes from above for clearness)
os.chdir(root_dir+'/'+project+'_TR_result/'+project+'_report')
assert not os.system('cp -r /PUBLIC/source/RNA/RefRNA/Report/template/TR/src   .')

render={}
if (ss == 'no'):
	os.system('cp /PUBLIC/source/RNA/RefRNA/Report/template/TR/Methods.TransRef.pdf src/Methods.TransRef.pdf')
	render['flag_strand']=False
else:
	os.system('cp /PUBLIC/source/RNA/RefRNA/Report/template/TR/Methods.TransRef.strand.pdf src/Methods.TransRef.pdf')
	render['flag_strand']=True

render['title']=project
render['flag_qc']=False
render['flag_tophat']=False
render['flag_as']=False
render['flag_novelgene']=False
render['flag_snp']=False
render['flag_htseq']=False
render['flag_rnaseqqc']=False
render['flag_diffexp']=False
render['flag_go']=False
render['flag_kegg']=False
render['flag_ppi']=False
render['flag_deu']=False

figure_order=0
table_order=0
html_order=0

render['flag_qc']=True
html_order+=1  #for raw data stat
html_order+=1  #for qc
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
render['figure2_1']=[]
render['figure2_2']=[]
render['figure2_3']=[]
render['table2_4']=[]
for eachsample in samples:
	assert not os.system('convert -resize 600 %s ./src/pictures/%s.error_rate_distribution.png' % (root_dir+'/QC_TR/'+eachsample+'/files/'+eachsample+'.Error.png',eachsample))
	assert not os.system('convert -resize 90 %s ./src/pictures/%s.error_rate_distribution.JPEG' % (root_dir+'/QC_TR/'+eachsample+'/files/'+eachsample+'.Error.png',eachsample))
	render['figure2_1'].append(["'"+'src/pictures/'+eachsample+'.error_rate_distribution.png'+"'","'"+'src/pictures/'+eachsample+'.error_rate_distribution.JPEG'+"'"])
	assert not os.system('convert -resize 600 %s ./src/pictures/%s.GC_content_distribution.png' % (root_dir+'/QC_TR/'+eachsample+'/files/'+eachsample+'.GC.png',eachsample))
	assert not os.system('convert -resize 90 %s ./src/pictures/%s.GC_content_distribution.JPEG' % (root_dir+'/QC_TR/'+eachsample+'/files/'+eachsample+'.GC.png',eachsample))
	render['figure2_2'].append(["'"+'src/pictures/'+eachsample+'.GC_content_distribution.png'+"'","'"+'src/pictures/'+eachsample+'.GC_content_distribution.JPEG'+"'"])
	assert not os.system('convert -resize 600 %s ./src/pictures/%s.raw_reads_classification.png' % (root_dir+'/QC_TR/'+eachsample+'/files/'+eachsample+'.pie3d.png',eachsample))
	assert not os.system('convert -resize 90 %s ./src/pictures/%s.raw_reads_classification.JPEG' % (root_dir+'/QC_TR/'+eachsample+'/files/'+eachsample+'.pie3d.png',eachsample))
	render['figure2_3'].append(["'"+'src/pictures/'+eachsample+'.raw_reads_classification.png'+"'","'"+'src/pictures/'+eachsample+'.raw_reads_classification.JPEG'+"'"])
	for eachLine in open(root_dir+'/QC_TR/'+eachsample+'/files/dataTable'):
		if eachLine.strip() != '':
			render['table2_4'].append(eachLine.strip().split())

render['flag_tophat']=True
render['flag_flow']=True
html_order+=1
render['html_order_tophat']=str(html_order)

table_order+=1
render['table_order_tophat']=str(table_order)
render['table3_1']=[]
temp = linecache.getline(root_dir+'/'+project+'_TR_result/'+project+'_results/3.Mapping/3.1.MapStat/MapStat.xls',1)
tophat_sample = temp.strip().split('\\t')[1:]
render['tophat_head']='</th><th>'.join(tophat_sample)
cont=open(root_dir+'/'+project+'_TR_result/'+project+'_results/3.Mapping/3.1.MapStat/MapStat.xls').readlines()
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
for eachsample in samples:
	assert not os.system('convert -resize 600 %s ./src/pictures/%s.Mapped_Region.png' % (root_dir+'/QC_TR/'+eachsample+'/files/'+eachsample+'.MR.png',eachsample))
	assert not os.system('convert -resize 90 %s ./src/pictures/%s.Mapped_Region.JPEG' % (root_dir+'/QC_TR/'+eachsample+'/files/'+eachsample+'.MR.png',eachsample))
	render['figure3_2'].append(["'"+'src/pictures/'+eachsample+'.Mapped_Region.png'+"'","'"+'src/pictures/'+eachsample+'.Mapped_Region.JPEG'+"'"])
	assert not os.system('convert -resize 600 %s ./src/pictures/%s.density.png' % (root_dir+'/Curve_TR/density/'+eachsample+'/'+eachsample+'.density.png',eachsample))
	assert not os.system('convert -resize 90 %s ./src/pictures/%s.density.JPEG' % (root_dir+'/Curve_TR/density/'+eachsample+'/'+eachsample+'.density.png',eachsample))
	render['figure3_3'].append(["'"+'src/pictures/'+eachsample+'.density.png'+"'","'"+'src/pictures/'+eachsample+'.density.JPEG'+"'"])
	assert not os.system('convert -resize 600 %s ./src/pictures/%s.ReadsVSlength.png' % (root_dir+'/Curve_TR/density/'+eachsample+'/'+eachsample+'.ReadsVSlength.png',eachsample))
	assert not os.system('convert -resize 90 %s ./src/pictures/%s.ReadsVSlength.JPEG' % (root_dir+'/Curve_TR/density/'+eachsample+'/'+eachsample+'.ReadsVSlength.png',eachsample))
	render['figure3_3'].append(["'"+'src/pictures/'+eachsample+'.ReadsVSlength.png'+"'","'"+'src/pictures/'+eachsample+'.ReadsVSlength.JPEG'+"'"])
'''

if set([1]).issubset(includes):
	code+='''
render['flag_as']=True
html_order+=1
render['html_order_as']=str(html_order)

figure_order+=1
render['figure_order_as']=str(figure_order)
render['figure4_1']=[]

assert not os.system('cp %s %s' % (root_dir+'/CAN_TR/CAN/ASprofile/AS.png',root_dir+'/'+project+'_TR_result/'+project+'_report/src/images/pic4.1.png'))

table_order+=1
render['table_order_as']=str(table_order)
render['table4_2']=[]
k=4
cont=open(root_dir+'/CAN_TR/CAN/ASprofile/'+samples[0]+'/'+samples[0]+'.fpkm.xls').readlines()
cont.pop(0)
for eachLine in cont:
	if (not eachLine.startswith('event_id')) and eachLine.strip() != '':
		temp=eachLine.split('\\t')
		render['table4_2'].append(temp)
		k-=1
	if k == 0:
		break

render['flag_novelgene']=True
html_order+=1
render['html_order_novogene']=str(html_order)

table_order==1
render['table_order_novel1']=str(table_order)
render['table5_1']=[]
k=4
cont=open(root_dir+'/CAN_TR/CAN/NovelGene/novelGene.gtf.xls').readlines()
cont.pop(0)
for eachLine in cont:
	temp=eachLine.split('\\t')
	render['table5_1'].append(temp)
	k-=1
	if k == 0:
		break

table_order+=1
render['table_order_novel2']=str(table_order)
render['table5_2']=[]
k=4
render['geneOpt_head']='</th><th>'.join(samples)
cont=open(root_dir+'/CAN_TR/CAN/NovelGene/geneStructOpt.xls').readlines()
cont.pop(0)
for eachLine in cont:
        if (not eachLine.startswith('Gene_id')) and eachLine.strip() !='':
                temp=eachLine.split('\\t')
                render['table5_2'].append(temp)
                k-=1
        if k == 0:
                break

'''

if set([1,2]).issubset(includes):
	code+='''
render['flag_snp']=True
html_order+=1
render['html_order_snp']=str(html_order)
table_order+=1
render['table_order_snp']=str(table_order)
render['table6_1']=[]
k=4
render['SNP_head']='</th><th>'.join(samples)
cont=open(root_dir+'/SNP_TR/SNP/ResultsQ30/SNPs.xls').readlines()
cont.pop(0)
for eachLine in cont:
	if (not eachLine.startswith('#')) and eachLine.strip() != '':
		temp=eachLine.split('\\t')
		snp='</td><td>'.join(temp)
		render['table6_1'].append(snp)
		k-=1
	if k == 0:
		break

'''

if set([1,4,5]).issubset(includes):
	code+='''
render['flag_htseq']=True
html_order+=1
render['html_order_htseq']=str(html_order)

table_order+=1
render['table_order_htseq1']=str(table_order)
render['table7_1']=[]
render['rpkm_stat_head']='</th><th>'.join(samples)
cont=open(root_dir+'/Diff_TR/Diff/rpkm.stat.xls').readlines()
cont.pop(0)
for eachLine in cont:
	if (not eachLine.startswith('RPKM')) and eachLine.strip() !='':
		temp=eachLine.split('\\t')
		rpkm_stat='</td><td>'.join(temp)
		render['table7_1'].append(rpkm_stat)

table_order+=1
render['table_order_htseq2']=str(table_order)
render['table7_2']=[]
k=4
render['rpkm_head']='</th><th>'.join(samples)
cont=open(root_dir+'/Diff_TR/Diff/rpkm.xls').readlines()
cont.pop(0)
for eachLine in cont:
	if (not eachLine.startswith('geneID')) and eachLine.strip() !='':
		temp=eachLine.split('\\t')
		rpkm='</td><td>'.join(temp)
		render['table7_2'].append(rpkm)
		k-=1
	if k == 0:
		break
'''

if set([1,4,6]).issubset(includes):
	code+='''
render['flag_rnaseqqc']=True
html_order+=1
render['html_order_rnaseqqc']=str(html_order)

figure_order+=1
render['figure_order_rnaseqqc1']=str(figure_order)
figure_order+=1
render['figure_order_rnaseqqc2']=str(figure_order)
render['figure8_1']=[]
render['figure8_2']=[]
for eachsample in samples:
	assert not os.system('convert -resize 600 %s ./src/pictures/%s.saturation.png' % (root_dir+'/Curve_TR/SaturationCurve/'+eachsample+'/'+eachsample+'.saturation.png',eachsample))
	assert not os.system('convert -resize 90 %s ./src/pictures/%s.saturation.JPEG' % (root_dir+'/Curve_TR/SaturationCurve/'+eachsample+'/'+eachsample+'.saturation.png',eachsample))
	render['figure8_1'].append(["'"+'src/pictures/'+eachsample+'.saturation.png'+"'","'"+'src/pictures/'+eachsample+'.saturation.JPEG'+"'"])
	assert not os.system('convert -resize 600 %s ./src/pictures/%s.Mean_coverage_distribution.png' % (root_dir+'/QC_TR/'+eachsample+'/files/'+eachsample+'.MC.png',eachsample))
	assert not os.system('convert -resize 90 %s ./src/pictures/%s.Mean_coverage_distribution.JPEG' % (root_dir+'/QC_TR/'+eachsample+'/files/'+eachsample+'.MC.png',eachsample))
	render['figure8_2'].append(["'"+'src/pictures/'+eachsample+'.Mean_coverage_distribution.png'+"'","'"+'src/pictures/'+eachsample+'.Mean_coverage_distribution.JPEG'+"'"])
'''

if set([1,4,5,6]).issubset(includes):
	code+='''
render['flag_rnaseqqc_cor']=True
figure_order+=1
render['figure_order_rnaseqqc_cor']=str(figure_order)
render['figure8_3']=[]
assert not os.system('convert -resize 600 %s %s' % (root_dir+'/Diff_TR/Diff/corr_plot/cor_pearson.png','./src/pictures/cor_pearson.png'))
assert not os.system('convert -resize 90 %s %s' % (root_dir+'/Diff_TR/Diff/corr_plot/cor_pearson.png','./src/pictures/cor_pearson.JPEG'))
render['figure8_3'].append(["'"+'./src/pictures/cor_pearson.png'+"'","'"+'./src/pictures/cor_pearson.JPEG'+"'"])
for png in glob.iglob(root_dir+'/Diff_TR/Diff/corr_plot/*.png'):
	if png !=root_dir+'/Diff_TR/Diff/corr_plot/cor_pearson.png':
		name=re.search(r'(/.*/)(.*)\.png',png).group(2)
		assert not os.system('convert -resize 600 %s %s' % (png,'./src/pictures/'+name+'.png'))
		assert not os.system('convert -resize 90 %s %s' % (png,'./src/pictures/'+name+'.JPEG'))
		render['figure8_3'].append(["'"+'src/pictures/'+os.path.basename(png)+"'","'"+'src/pictures/'+name+'.JPEG'+"'"])
'''

if set([1,4,5]).issubset(includes):
	code+='''
render['flag_diffexp']=True
html_order+=1
render['html_order_diffexp']=str(html_order)

figure_order+=1
render['figure_order_diffexp1']=str(figure_order)
render['figure9_1']=[]
assert not os.system('convert -resize 600 %s ./src/pictures/boxplot.png' % (root_dir+'/Diff_TR/Diff/boxplot.png'))
assert not os.system('convert -resize 90 %s ./src/pictures/boxplot.JPEG' % (root_dir+'/Diff_TR/Diff/boxplot.png'))
assert not os.system('convert -resize 600 %s ./src/pictures/density.png' % (root_dir+'/Diff_TR/Diff/density.png'))
assert not os.system('convert -resize 90 %s ./src/pictures/density.JPEG' % (root_dir+'/Diff_TR/Diff/density.png'))
render['figure9_1'].append(["'"+'src/pictures/boxplot.png'+"'","'"+'src/pictures/boxplot.JPEG'+"'"])
render['figure9_1'].append(["'"+'src/pictures/density.png'+"'","'"+'src/pictures/density.JPEG'+"'"])

render['flag_venn']=flag_venn
table_order+=1
render['table_order_diffexp']=str(table_order)
render['table9_2']=[]
com=compare_names[0].split('vs')
render['com1']=com[0]
render['com2']=com[1]
k=4
for eachLine in open(root_dir+'/Diff_TR/Diff/'+compare_names[0]+'/'+compare_names[0]+'.diffgene.xls'):
	if (not eachLine.startswith('Gene_id')) and eachLine.strip() != '':
		temp=eachLine.split('\\t')
		render['table9_2'].append(temp)
		k-=1
	if k == 0:
		break

figure_order+=1
render['figure_order_diffexp2']=str(figure_order)
render['figure9_3']=[]
for each in compare_names:
	assert not os.system('convert -resize 600 %s ./src/pictures/%s.Volcanoplot.png' % (root_dir+'/Diff_TR/Diff/'+each+'/'+each+'.volcano.png',each))
	assert not os.system('convert -resize 90 %s ./src/pictures/%s.Volcanoplot.JPEG' % (root_dir+'/Diff_TR/Diff/'+each+'/'+each+'.volcano.png',each))
	render['figure9_3'].append(["'"+'src/pictures/'+each+'.Volcanoplot.png'+"'","'"+'src/pictures/'+each+'.Volcanoplot.JPEG'+"'"])

figure_order+=1
render['figure_order_diffexp3']=str(figure_order)
render['figure9_4']=[]
assert not os.system('convert -resize 600 %s ./src/pictures/heatCluster.png' % (root_dir+'/Diff_TR/Diff/heatCluster.png'))
assert not os.system('convert -resize 90 %s ./src/pictures/heatCluster.JPEG' % (root_dir+'/Diff_TR/Diff/heatCluster.png'))
render['figure9_4'].append(["'"+'src/pictures/heatCluster.png'+"'","'"+'src/pictures/heatCluster.JPEG'+"'"])
assert not os.system('convert -resize 600 %s ./src/pictures/h_show_plots.png' % (root_dir+'/Diff_TR/Diff/union_for_cluster_hclust/h_show_plots.png'))
assert not os.system('convert -resize 90 %s ./src/pictures/h_show_plots.JPEG' % (root_dir+'/Diff_TR/Diff/union_for_cluster_hclust/h_show_plots.png'))
render['figure9_4'].append(["'"+'src/pictures/h_show_plots.png'+"'","'"+'src/pictures/h_show_plots.JPEG'+"'"])

if flag_venn:
	figure_order+=1
	render['figure_order_diffexp4']=str(figure_order)
	render['figure9_5']=[]
	for each in venn_cluster_vs_names:
		assert not os.system('convert -resize 600 %s ./src/pictures/%s.DEG_Venn_diagram.png' % (root_dir+'/Diff_TR/Diff/Venn/'+each+'/'+each+'.venn.png',each))
		assert not os.system('convert -resize 90 %s ./src/pictures/%s.DEG_Venn_diagram.JPEG' % (root_dir+'/Diff_TR/Diff/Venn/'+each+'/'+each+'.venn.png',each))
		render['figure9_5'].append(["'"+'src/pictures/'+each+'.DEG_Venn_diagram.png'+"'","'"+'src/pictures/'+each+'.DEG_Venn_diagram.JPEG'+"'"])
'''

if set([1,4,5,7]).issubset(includes):
	code+='''
render['flag_go']=True
html_order+=1
render['html_order_go']=str(html_order)

table_order+=1
render['table_order_go']=str(table_order)
render['table10_1']=[]
k=4
for eachLine in open(root_dir+'/GOSeq_TR/'+compare_names[0]+'/'+compare_names[0]+'.GO_enrichment_result_up_down.xls'):
	if (not eachLine.startswith('GO_accession')) and eachLine.strip() != '':
		temp=eachLine.split('\\t')
		render['table10_1'].append(temp[:7])
		k-=1
	if k == 0:
		break

figure_order+=1
render['figure_order_go1']=str(figure_order)
render['figure10_2']=[]
figure_order+=1
render['figure_order_go2']=str(figure_order)
render['figure10_3']=[]

for each in compare_names:
	os.system('convert -resize 600 %s ./src/pictures/%s.DEG_Enriched_GO_classification.png' % (root_dir+'/GOSeq_TR/'+each+'/'+each+'.bar_graph.png',each))
	os.system('convert -resize 90 %s ./src/pictures/%s.DEG_Enriched_GO_classification.JPEG' % (root_dir+'/GOSeq_TR/'+each+'/'+each+'.bar_graph.png',each))
	render['figure10_2'].append(["'"+'src/pictures/'+each+'.DEG_Enriched_GO_classification.png'+"'","'"+'src/pictures/'+each+'.DEG_Enriched_GO_classification.JPEG'+"'"])
	os.system('convert -resize 600 %s ./src/pictures/%s.DEG_Enriched_GO_updown.png' % (root_dir+'/GOSeq_TR/'+each+'/'+each+'.bar_graph_updown.png',each))
	os.system('convert -resize 90 %s ./src/pictures/%s.DEG_Enriched_GO_updown.JPEG' % (root_dir+'/GOSeq_TR/'+each+'/'+each+'.bar_graph_updown.png',each))
	render['figure10_2'].append(["'"+'src/pictures/'+each+'.DEG_Enriched_GO_updown.png'+"'","'"+'src/pictures/'+each+'.DEG_Enriched_GO_updown.JPEG'+"'"])
	os.system('convert -resize 600 %s ./src/pictures/%s.DEG_Enriched_GO_bp_DAG.png' % (root_dir+'/GOSeq_TR/'+each+'/'+each+'.GObp_DAG.png',each))
	os.system('convert -resize 90 %s ./src/pictures/%s.DEG_Enriched_GO_bp_DAG.JPEG' % (root_dir+'/GOSeq_TR/'+each+'/'+each+'.GObp_DAG.png',each))
	render['figure10_3'].append(["'"+'src/pictures/'+each+'.DEG_Enriched_GO_bp_DAG.png'+"'","'"+'src/pictures/'+each+'.DEG_Enriched_GO_bp_DAG.JPEG'+"'"])
	os.system('convert -resize 600 %s ./src/pictures/%s.DEG_Enriched_GO_cc_DAG.png' % (root_dir+'/GOSeq_TR/'+each+'/'+each+'.GOcc_DAG.png',each))
	os.system('convert -resize 90 %s ./src/pictures/%s.DEG_Enriched_GO_cc_DAG.JPEG' % (root_dir+'/GOSeq_TR/'+each+'/'+each+'.GOcc_DAG.png',each))
	render['figure10_3'].append(["'"+'src/pictures/'+each+'.DEG_Enriched_GO_cc_DAG.png'+"'","'"+'src/pictures/'+each+'.DEG_Enriched_GO_cc_DAG.JPEG'+"'"])
	os.system('convert -resize 600 %s ./src/pictures/%s.DEG_Enriched_GO_mf_DAG.png' % (root_dir+'/GOSeq_TR/'+each+'/'+each+'.GOmf_DAG.png',each))
	os.system('convert -resize 90 %s ./src/pictures/%s.DEG_Enriched_GO_mf_DAG.JPEG' % (root_dir+'/GOSeq_TR/'+each+'/'+each+'.GOmf_DAG.png',each))
	render['figure10_3'].append(["'"+'src/pictures/'+each+'.DEG_Enriched_GO_mf_DAG.png'+"'","'"+'src/pictures/'+each+'.DEG_Enriched_GO_mf_DAG.JPEG'+"'"])

'''

if set([1,4,5,8]).issubset(includes):
	code+='''
render['flag_kegg']=True
html_order+=1
render['html_order_kegg']=str(html_order)

table_order+=1
render['table_order_kegg']=str(table_order)
flag=0
k=4
sig_map=''
kegg_png=[]
render['table11_1']=[]
for eachLine in open(root_dir+'/KOBAS_TR/'+compare_names[0]+'/add.'+compare_names[0]+'.identify.xls'):
	if eachLine.startswith('#Term'):
		flag=1
		continue
	if flag == 1:
		if eachLine.strip() != '':
			temp=eachLine.split('\\t')
			if sig_map == '':
				sig_map=temp[2].strip()
				sig_map_png=root_dir+'/KOBAS_TR/'+compare_names[0]+'/src/'+sig_map+'.png'
				if os.path.isfile(sig_map_png):
					sig_map=temp[2].strip()
					sig_map_png=root_dir+'/KOBAS_TR/'+compare_names[0]+'/src/'+sig_map+'.png'
					kegg_png.append(sig_map_png)
					sig_map = ''
				else:
					sig_map = ''
			render['table11_1'].append(temp[:7])
			k-=1
		if k == 0:
			break
figure_order+=1
render['figure_order_kegg1']=str(figure_order)
render['figure11_2']=[]
render['figure20_lable']=[]
for each in compare_names:
	os.system('convert -resize 600 %s ./src/pictures/%s.DEG_enriched_KEGG_pathway_scatterplot.png' % (root_dir+'/KOBAS_TR/'+each+'/'+each+'.DEG_enriched_KEGG_pathway_scatterplot.png',each))
	os.system('convert -resize 90 %s ./src/pictures/%s.DEG_enriched_KEGG_pathway_scatterplot.JPEG' % (root_dir+'/KOBAS_TR/'+each+'/'+each+'.DEG_enriched_KEGG_pathway_scatterplot.png',each))
	render['figure11_2'].append(["'"+'src/pictures/'+each+'.DEG_enriched_KEGG_pathway_scatterplot.png'+"'","'"+'src/pictures/'+each+'.DEG_enriched_KEGG_pathway_scatterplot.JPEG'+"'"])
	render['figure20_lable'].append(each)
render['figure20_lable']=','.join(render['figure20_lable'])

figure_order+=1
render['figure_order_kegg2']=str(figure_order)
render['figure11_3']=[]
for png in kegg_png:
	name=re.search(r'(/.*/)(.*)\.png',png).group(2)
	assert not os.system('convert -resize 600 %s %s' % (png,'./src/pictures/'+name+'.png'))
	assert not os.system('convert -resize 90 %s %s' % (png,'./src/pictures/'+name+'.JPEG'))
	render['figure11_3'].append(["'"+'src/pictures/'+os.path.basename(png)+"'","'"+'src/pictures/'+name+'.JPEG'+"'"])
'''

if set([1,4,5,9]).issubset(includes):
	code+='''
render['flag_ppi']=True
html_order+=1
render['html_order_ppi']=str(html_order)

figure_order+=1
render['figure_order_ppi']=str(figure_order)
'''
if set([1,3]).issubset(includes):
	code+='''
import re
render['flag_deu']=True
html_order+=1
render['html_order_deu']=str(html_order)

table_order+=1
render['table_order_deu']=str(table_order)
render['table13']=[]
k=4
for eachLine in open(root_dir+'/DEXseq_TR/DEXseq/'+compare_names[0]+'/'+compare_names[0]+'.DEUdiff.xls'):
	if (not eachLine.startswith('geneID')) and eachLine.strip() != '':
		temp=eachLine.split('\\t')
		render['table13'].append(temp)
		k-=1
	if k == 0:
		break

assert not os.system('ls %s |while read a;do echo "$a">>%s;exit;done' % (root_dir+'/DEXseq_TR/DEXseq/'+compare_names[0]+'/'+compare_names[0]+'.DEU/files/*expression.svg',root_dir+'/'+project+'_TR_result/'+project+'_report/src/pic.log'))
assert not os.system('ls %s |while read a;do echo "$a">>%s;exit;done' % (root_dir+'/DEXseq_TR/DEXseq/'+compare_names[0]+'/'+compare_names[0]+'.DEU/files/*counts.svg',root_dir+'/'+project+'_TR_result/'+project+'_report/src/pic.log'))
assert not os.system('ls %s |while read a;do echo "$a">>%s;exit;done' % (root_dir+'/DEXseq_TR/DEXseq/'+compare_names[0]+'/'+compare_names[0]+'.DEU/files/*splicing.svg',root_dir+'/'+project+'_TR_result/'+project+'_report/src/pic.log'))
assert not os.system('ls %s |while read a;do echo "$a">>%s;exit;done' % (root_dir+'/DEXseq_TR/DEXseq/'+compare_names[0]+'/'+compare_names[0]+'.DEU/files/*transcripts.svg',root_dir+'/'+project+'_TR_result/'+project+'_report/src/pic.log'))

figure_order+=1
render['figure_order_deu']=str(figure_order)
render['figure14']=[]
for each_svg in open(root_dir+'/'+project+'_TR_result/'+project+'_report/src/pic.log'):
	each_svg = each_svg.rstrip()
	svg = os.path.basename(each_svg)
	if "(" in svg:
		newSvgName = svg.replace("(",".").replace(")",".").rstrip('.svg') #wash out "()" and suffix name
	else: newSvgName = svg.rstrip('.svg')
	if "(" in svg:
		assert not os.system("convert -resize 600 '%s' %s" % (each_svg,'./src/pictures/'+ newSvgName + '.png'))
		assert not os.system("convert -resize 90 '%s' %s" % (each_svg,'./src/pictures/'+ newSvgName  +'.JPEG'))
		render['figure14'].append(["'"+'./src/pictures/'+ newSvgName +'.png'+"'","'"+'./src/pictures/'+ newSvgName +'.JPEG'+"'"])
	else:
		assert not os.system("convert -resize 600 %s %s" % (each_svg,'./src/pictures/'+ newSvgName + '.png'))
		assert not os.system("convert -resize 90 %s %s" % (each_svg,'./src/pictures/'+ newSvgName  +'.JPEG'))
		render['figure14'].append(["'"+'./src/pictures/'+ newSvgName +'.png'+"'","'"+'./src/pictures/'+ newSvgName +'.JPEG'+"'"])
'''
#render['flag_strand']=flag_strand
if set([1]).issubset(includes):
  code+='''
flag_strand=False
if ss == 'yes' or ss == 'reverse':
	flag_strand = True
if flag_strand:
	html_order+=1
	render['html_order_anti']=str(html_order)
	
	table_order+=1
	render['table_order_anti1']=str(table_order)
	render['table14_1']=[]
	k=4
	for eachLine in open(root_dir+'/CAN_TR/CAN/stranded_specific/cisNAT/'+samples[0]+'.nat.xls'):
		if (not eachLine.startswith('chromosome')) and eachLine.strip() != '':
			temp=eachLine.split('\\t')
			render['table14_1'].append(temp[1:15])
			k-=1
		if k == 0:
			break

	table_order+=1
	render['table_order_anti2']=str(table_order)
	render['table14_2']=[]
	k=4
	for eachLine in open(root_dir+'/CAN_TR/CAN/stranded_specific/cisNAT/'+samples[0]+'.nat.stat.xls'):
		if (not eachLine.startswith('chromosome')) and eachLine.strip() != '':
			temp=eachLine.split('\\t')
			render['table14_2'].append(temp)
			k-=1
		if k == 0:
			break

	html_order+=1
	render['html_order_tss']=str(html_order)
	
	table_order+=1
	render['table_order_tss']=str(table_order)
	render['table15']=[]
	k=4
	for eachLine in open(root_dir+'/CAN_TR/CAN/stranded_specific/TSS.xls'):
		if (not eachLine.startswith('id')) and eachLine.strip() != '':
			temp=eachLine.split('\\t')
			render['table15'].append(temp[2:])
			k-=1
		if k == 0:
			break
'''

code+='''
t=loader.get_template('TR_report.html')
c=Context(render)
html=t.render(c)
open(project+'_Report.html','w').write(html)
assert not os.system("/PUBLIC/software/RNA/wkhtmltopdf/wkhtmltopdf --page-width 350mm --page-height 495mm -n --print-media-type --footer-center '[page] / [topage]' "+ project+'_Report.html '+project+'_Report.pdf')
assert not os.system("sed -i 's/#00//g' "+project+'_Report.pdf')
'''

code+='''
os.chdir(root_dir)
assert not os.system('perl /PUBLIC/source/RNA/RefRNA/TransRef/data_sum/data_sum.pl -n %s -p %s -ss %s -dir %s')
''' % (sample,project+'_TR',ss,root_dir)

open('TR_step6_Report.py','w').write(code)
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
''' % (fq,fa,root_dir+'/CAN_TR/CAN/NovelGene/novelcombined.gtf',root_dir+'/CAN_TR/CAN/NovelGene/GO/gene.go',sample,root_dir,project)

if argv['genenamefile']:
	code+='''
genenamefile='%s'
''' % (genenamefile)
else:
	code+='''
genenamefile='%s'
''' % (root_dir+'/Blast_TR/Blast_Swissprot/diffgene_union.genenames')

code+='''
assert not os.system('mkdir data_release')
data_release=open(root_dir+'/data_release/data_release.sh','w')
data_release.write('mkdir %s\\n' % (root_dir+'/data_release/data_give'))
dir=root_dir+'/data_release/data_give'
data_release.write('echo "###############prepare SuppFiles#####################"\\n')
data_release.write('mkdir %s\\n' % (root_dir+'/'+project+'_TR_result/'+project+'_results/0.SuppFiles'))
SuppFiles_dir=root_dir+'/'+project+'_TR_result/'+project+'_results/0.SuppFiles'
data_release.write('cp /PUBLIC/source/RNA/RefRNA/TransRef/data_release/SuppFiles.README.txt %s\\n' % (SuppFiles_dir))
data_release.write('cp %s %s\\n' % (genenamefile,SuppFiles_dir+'/gene.description.xls'))
data_release.write('perl /PUBLIC/source/RNA/RefRNA/TransRef/data_release/GO_class.pl %s /PUBLIC/software/RNA/GOseq/gene_ontology.1_2.obo.ab %s\\n' % (goann,SuppFiles_dir+'/gene.goannot.xls'))
data_release.write('perl /PUBLIC/source/RNA/RefRNA/TransRef/data_release/extractInfo_from_GTF.pl -i %s -o %s\\n' % (gtf,SuppFiles_dir+'/gene.info.xls'))
data_release.write('perl /PUBLIC/source/RNA/RefRNA/DGE/Enrichment/extractcDNAfromFA.pl %s %s %s\\n' % (gtf,fa,SuppFiles_dir+'/gene.fasta'))
#data_release.write('gzip %s\\n' % (SuppFiles_dir+'/gene.fasta'))
data_release.write('echo "###############copy Novofinder#####################"\\n')
data_release.write('cp /PUBLIC/source/RNA/RefRNA/TransRef/data_release/Novofinder.exe   %s\\n' % (root_dir+'/'+project+'_TR_result/Novofinder.exe'))
data_release.write('cp /PUBLIC/source/RNA/RefRNA/TransRef/data_release/Novofinder_manual.pdf %s\\n' % (root_dir+'/'+project+'_TR_result'))
data_release.write('echo "###############tar Report#####################"\\n')
data_release.write('cd %s\\n' % (root_dir))
data_release.write("tar -cvzf %s %s\\n" % (project + '_TR.tar.gz',project+'_TR_result/'))
data_release.write('mv %s %s\\n' % (project + '_TR.tar.gz',root_dir+'/data_release/data_give/'))
data_release.write('cd %s\\n' % (root_dir))
data_release.write('echo "###############prepare clean data#####################"\\n')
data_release.write('ln -s %s %s\\n' % (root_dir+'/QC_TR/cleandata',dir+'/cleandata'))
data_release.write('echo "###############prepare raw data#####################"\\n')
data_release.write('mkdir %s\\n' % (dir+'/rawdata'))
data_release.write('cp %s %s\\n' % (root_dir+'/QC_TR/rd_md5.txt',dir+'/rawdata/'))

samples=sample.split(',')
for eachsample in samples:
	data_release.write('ln -s %s %s\\n' % (root_dir+'/QC_TR/'+eachsample+'/'+eachsample+'_1.fq.gz',dir+'/rawdata/'))
	data_release.write('ln -s %s %s\\n' % (root_dir+'/QC_TR/'+eachsample+'/'+eachsample+'_2.fq.gz',dir+'/rawdata/'))

data_release.write('echo "###############prepare IGV_data#####################"\\n')
data_release.write('mkdir %s\\n' % (dir+'/IGV_data'))
data_release.write('cp /PUBLIC/source/RNA/RefRNA/TransRef/data_release/IGV_data.README.txt %s\\n' % (dir+'/IGV_data'))
data_release.write('cp %s %s\\n' % (fa,dir+'/IGV_data/genome.fasta'))
data_release.write('gzip %s\\n' % (dir+'/IGV_data/genome.fasta'))
data_release.write('ln -s %s %s\\n' % (root_dir+'/QC_TR/bam',dir+'/IGV_data/bam'))
data_release.close()

os.chdir('data_release')
assert not os.system('qsub -V -cwd -l vf=1G -l p=1 data_release.sh')
os.chdir(root_dir)
'''

open('TR_step7_data_release.py','w').write(code)
#################################################################################################
#for byebye
byebye=open(root_dir+'/byebye.sh','w')
byebye.write('sh /PUBLIC/source/RNA/RefRNA/TransRef/Pipeline/byebye_v2.1.sh -dir %s -sample %s -pipline TR\n' % (root_dir,sample))
byebye.close()
