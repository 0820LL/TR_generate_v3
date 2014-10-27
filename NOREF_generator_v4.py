#noref pipline
#2013.01.06
#modified 2013.01.31
#modified 2013.5.24 by hanrui
#modified 2014.8.21 by yangying
#guoyang
import sys
import os
import os.path
import re
import argparse
root_dir=os.getcwd()

##################################################################################
#parse the arguments
parser = argparse.ArgumentParser(description="noref pipeline v1.5")		
parser.add_argument('--project',help='project name, maybe same with the name of root dir, which will be displayed in the final report title',required=True)
parser.add_argument('--sample',help="sample names(sample1,sample2,...)  warning: order sensitive!",required=True)
parser.add_argument('--seq',help="the original directory of the raw fastq reads for QC, [REQUIRED]",default=None)
parser.add_argument('--ad',help="the original directory of the adapter list for QC",default=None)
parser.add_argument('--mapfile',help="mapfile files (mapfile1,mapfile2.The first column of mapfile is the NH number ,the second column of mapfile is the sample name)",required=True)
parser.add_argument('--raw_dir',help="the original directory of the raw fastq reads keep order in line with mapfile files (raw_dir1,raw_dir2)",default=None)
#parser.add_argument('--len',help="length of reads(n,n,...) for SNP, keep order in line with '--sample'",default=None)
parser.add_argument('--ex',help="the steps you do not wanna perform",default=None)
parser.add_argument('--group',help="sample classification, e.g. sample1:sample2,sample3 for differential expression analysis",default=None)
parser.add_argument('--groupname',help="group names, (default: group1,group2,... with 'sample repeat' or sample names directly without 'sample repeat'),for differential expression analysis",default=None)
parser.add_argument('--compare',help="group comparison strategy(1:2,1:3,...) warning: number only,for differential expression analysis",default=None)
parser.add_argument('--matrix',help="the score matrices for estscan program.",default='Hs')
parser.add_argument('--codon',help="the id indicates the genetic codon table.",default='1')
parser.add_argument('--ss',help='Probability of generating a read from the forward strand of a transcript.',choices=['0.0','0.5','1.0'],default='0.5')
parser.add_argument('--jc',help='(--jaccard_clip)set 1 if you have paired reads and you expect high gene density with UTR overlap (use FASTQ input file format for reads)(Default : 0)',choices=['0','1'],default='0')
parser.add_argument('--minkmercov',help='min count for K-mers to be assembled by Inchworm(Default : 2)',default='2')
parser.add_argument('--minglue',help='min number of reads needed to glue two inchworm contigs together(Default : 2)',default='2')
parser.add_argument('--db',help='Database (NR) used; choose "nr" for the whole database or "eu" for euKaryote only.',choices=['nr','eu'],default='eu')
parser.add_argument('--venn_cluster',help='compares groups for venn separately',default=None)
parser.add_argument('--generate_adapter',help='whether to generate the adpepter list y or n',default='n')
parser.add_argument('--index',help='the index number of the fastq files,if you want to generate the adapter list,the parameter will be necessary,diffrent samples split by ,',default=None)
parser.add_argument('--ppi_number',help='the number of reference for PPI',default=None)

#extract, wash and check the parameter
argv=vars(parser.parse_args())
project=argv['project'].strip()
samples=[each.strip() for each in argv['sample'].strip().split(',') if each.strip() != '']
sample=','.join(samples)
assert not os.system('perl /PUBLIC/source/RNA/noRef/Pipeline_noRef/check.sample_name.pl -sample %s ' % (sample))
assert len(set(samples)) == len(samples)
generate_adapter=argv['generate_adapter']
assert argv['seq']!=None or argv['raw_dir'] != None
mapfiles=[each.strip() for each in argv['mapfile'].strip().split(',') if each.strip() != '']
mapfile=' '.join(mapfiles)
assert not os.system('cat %s |sort -u |sed \'/^$/d\' >%s' % (mapfile,root_dir+'/libraryID'))
if argv['seq']!=None:
	sequences=[]
	seq=argv['seq'].strip()
	for eachsample in samples:
		sequences.append([seq+'/'+eachsample+'_1.fq.gz',seq+'/'+eachsample+'_2.fq.gz'])
	if generate_adapter =='n':
		adaptors=[]
		for eachsample in samples:
			adaptors.append([seq+'/'+eachsample+'_1.adapter.list.gz',seq+'/'+eachsample+'_2.adapter.list.gz'])
else:
	assert not os.system('mkdir raw_data')
	assert not os.system('perl /PUBLIC/source/RNA/noRef/QC/bin/ln_raw_data.pl %s %s pe raw_data' %(argv['mapfile'],argv['raw_dir']))
	sequences=[]
	for eachsample in samples:
		sequences.append([root_dir+'/raw_data/'+eachsample+'_1.fq.gz',root_dir+'/raw_data/'+eachsample+'_2.fq.gz'])
	if generate_adapter =='n':
		adaptors=[]
		for eachsample in samples:
			adaptors.append([root_dir+'/raw_data/'+eachsample+'_1.adapter.list.gz',root_dir+'/raw_data/'+eachsample+'_2.adapter.list.gz'])


assert len(samples) == len(sequences)
for each1 in sequences:
	for each2 in each1:
		assert os.path.isfile(each2)
sequence=','.join([':'.join(each) for each in sequences])

if argv['ad'] != None :
	ad=argv['ad'].strip()
	adaptors=[]
	for eachsample in samples:
		adaptors.append([ad+'/'+eachsample+'_1.adapter.list.gz',ad+'/'+eachsample+'_2.adapter.list.gz'])
	assert len(samples) == len(adaptors)
	for each1 in adaptors:
		assert len(each1) == 2
		if (each1[0].strip() != '') or (each1[1].strip() != ''):
			for each2 in each1:
				assert os.path.isfile(each2)
	adaptor=','.join([':'.join(each) for each in adaptors])

else:
	if argv['raw_dir'] !=None:
		if generate_adapter =='n':
			assert len(samples) == len(adaptors)
			for each1 in adaptors:
				assert len(each1) == 2
				if (each1[0].strip() != '') or (each1[1].strip() != ''):
					for each2 in each1:
						assert os.path.isfile(each2)
				adaptor=','.join([':'.join(each) for each in adaptors])

		else:
			assert argv['index']
			indexs=[each.strip() for each in argv['index'].strip().split(',')]	
			assert len(samples) == len(indexs)
			for each1 in indexs:
				index=','.join(indexs)


			adaptor=','.join([':' for i in range(len(samples))])

	else:	
		if generate_adapter !='n':
			assert argv['index']
			indexs=[each.strip() for each in argv['index'].strip().split(',')]	
			assert len(samples) == len(indexs)
			for each1 in indexs:
				index=','.join(indexs)
		adaptor=','.join([':' for i in range(len(samples))])
	
if argv['seq']!=None:
	if generate_adapter =='n':
		if argv['ad'] == None:
			print 'Error:  the parameters --ad and --generate_adapter are not consensus!\n'
			exit()

################################
all_content=set([1,2,3,4,5,6,7,8,9,10])
if argv['ex'] != None:
	excludes=argv['ex'].strip().strip(',').strip().split(',')
	excludes=[int(each.strip()) for each in excludes]
	for each1 in excludes:
		assert each1 in all_content
else:
	excludes=[] #list
includes=all_content-set(excludes) #set
################################
if set([1]).issubset(includes):
	ss=argv['ss'].strip()
	jc=argv['jc'].strip()
	minkmercov=argv['minkmercov'].strip()
	minglue=argv['minglue'].strip()

################################
'''
if set([1,3]).issubset(includes):
	assert argv['len']
	lengths=argv['len'].strip().strip(',').strip().split(',')
	lengths=[each.strip() for each in lengths]
	if len(lengths) == 1:
		lengths=lengths*len(samples)
	for each1 in lengths:
		assert each1.isdigit()
	length=','.join(lengths)
	assert len(lengths) == len(samples)
'''
################################	

if set([1,5,6]).issubset(includes):
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
		#assert len(groups_iter) == len(set(groups_iter))
		assert set(groups_iter).issubset(samples)
	group_iter=','.join(groups_iter)
	
	if argv['groupname'] == None:
		if flag_repeat == False:
			groupnames=groups
		else:
			groupnames=['group'+str(k+1) for k in range(len(groups))]		
	else:
		groupnames=[each.strip() for each in argv['groupname'].split(',') if each.strip() != '']
		assert len(groupnames) == len(groups)
	groupname=','.join(groupnames)
	#############################
	assert argv['compare']
	compares=[each.strip().split(':') for each in argv['compare'].strip().split(',') if each.strip() != '']
	M=[]
	for each1 in compares:
		assert len(each1) == 2
		for each2 in each1:
			assert each2.isdigit()
			M.append(int(each2))
	assert max(M) <= len(groupnames)
	assert min(M) > 0	
	compare=','.join([':'.join(each) for each in compares])
	temp2=[]
	for each1 in compares:
		temp1=[]
		for each2 in each1:
			temp1.append(groupnames[int(each2)-1])
		temp2.append(':'.join(temp1))
	compare_name=','.join(temp2)
	#############################################
	if argv['venn_cluster'] != None:
		venn_cluster=argv['venn_cluster'].strip()
		com_pairs=compare.split(',')
		venn_clusters=[each.split('_') for each in venn_cluster.split(',')]
		temp1=[]
		for each1 in venn_clusters:
			temp2=[]
			for each2 in each1:
				assert each2 in com_pairs
				temp3=each2.split(':')
				assert len(temp3) == 2
				temp2.append(groupnames[int(temp3[0])-1]+':'+groupnames[int(temp3[1])-1])
			temp1.append('_'.join(temp2))	
		venn_cluster_name=','.join(temp1)		
	else:
		venn_cluster=compare.replace(',','_')
		venn_cluster_name=compare_name.replace(',','_')
##########################################################################
if set([1,5,6,10]).issubset(includes):
	assert argv['ppi_number'] !=None
	ppi_number=argv['ppi_number']
########################################################################################
if set([1,2]).issubset(includes):
	matrix=argv['matrix'].strip()
	if (matrix not in ['At','Hs','Mm','Rn','Dm','Dr','Os','Zm']) and (not os.path.isfile(matrix)):
		sys.exit('matrix fail...')
	codon=argv['codon'].strip()
	if codon not in ['1','2','3','4','5','6','9','10','11','12','13','14','15','16','21','22','23']:
		sys.exit('codon fail...')
	db=argv['db']
###########################################################################	
#parameters display
display=open('parameters_display.txt','w')
display.write('project: '+project+'\n')
for i,eachsample in enumerate(samples):
	display.write('sample: '+eachsample+'\n')
	display.write('\traw data1: '+sequences[i][0]+'\n')
	display.write('\traw data2: '+sequences[i][1]+'\n')
	if argv['ad'] != None:
		display.write('\tadaptor list1: '+adaptors[i][0]+'\n')
		display.write('\tadaptor list2: '+adaptors[i][1]+'\n')
	if argv['generate_adapter'] !='n':
		display.write('\tindex list: '+indexs[i]+'\n')

	#if set([1,3]).issubset(includes):
		#display.write('\tlength: '+lengths[i]+'\n')
display.write('\n\n')

if set([1]).issubset(includes):
	display.write('trinity parameters:\n')
	display.write('\tStrand-specific RNA-Seq read orientation: %s\n' % ss)
	display.write('\tjaccard_clip: %s\n' % jc)
	display.write('\tmin count for K-mers to be assembled by Inchworm: %s\n' % minkmercov)
	display.write('\tmin number of reads needed to glue two inchworm contigs together: %s\n\n\n' % minglue)


if set([1,5,6]).issubset(includes):
	display.write('group:samples\n')
	for i,eachgroup in enumerate(groupnames):
		display.write(eachgroup+': '+str(groups[i])+'\n')
	display.write('\n')
	display.write('group comparison strategy:\n')
	display.write(compare+'\n')
	display.write(compare_name+'\n')
	display.write('\n\n')
	display.write('venn cluster:\n')
	display.write(venn_cluster.replace('_',' ').replace(',','\n'))
	display.write('venn cluster name:\n')
	display.write(venn_cluster_name.replace('_',' ').replace(',','\n'))
	display.write('\n\n\n')
	
	
if set([1,2]).issubset(includes):
	display.write('matrix: '+matrix+'\n')
	display.write('codon: '+codon+'\n')
	db_detail={'nr':'for the whole database','eu':'for eukaryote only'}
	display.write('database: '+db+' ('+db_detail[db]+')\n')
	display.write('E-value: 1e-5\n\n\n')

contents_detail={0:'QC (default)',1:'trinity assembly',2:'function annotation & CDS prediction',3:'SNP',4:'SSR',5:'RSEM',6:'DIFF_EXP',7:'RNAseq QC',8:'KEGG enrichment',9:'GO enrichment',10:'PPI'}
display.write('content(s) included:\n')
display.write('\tcontent0: QC (default)\n')
for eachone in includes:
	display.write('\tcontent%s: %s\n' % (eachone,contents_detail[eachone]))

display.write('\n\n\n\n\n\n\n\n##########################################\n')
display.write(sample+'\n')
display.write(str(samples)+'\n')
display.write(sequence+'\n')
display.write(str(sequences)+'\n')
if argv['ad'] != None:
	display.write(str(adaptors)+'\n')
display.write(adaptor+'\n')
if set([1]).issubset(includes):
	display.write(ss+'\n')
	display.write(jc+'\n')
	display.write(minkmercov+'\n')
	display.write(minglue+'\n')
'''
if set([1,3]).issubset(includes):
	display.write(length+'\n')
	display.write(str(lengths)+'\n')
'''
if set([1,5,6]).issubset(includes):
	display.write(group+'\n')
	display.write(str(groups)+'\n')
	display.write(groupname+'\n')
	display.write(str(groupnames)+'\n')
	display.write(compare+'\n')
	display.write(str(compares)+'\n')
	display.write(compare_name+'\n')
	display.write(venn_cluster+'\n')
	display.write(venn_cluster_name+'\n')

if argv['ppi_number'] != None:
	display.write('ppi_number: %s\n\n' % (ppi_number))
	
	
if set([1,2]).issubset(includes):
	display.write(matrix+'\n')
	display.write(codon+'\n')
	display.write(db+'\n')
	
display.write(str(excludes)+'\n')
display.write(str(includes)+'\n')
display.close()

###########################################################################
#for QC (run)
code='''
import os

sample='%s'
samples=sample.split(',')
sequence='%s'
sequences=[each.split(':') for each in sequence.split(',')]
adaptor='%s'
adaptors=[each.split(':') for each in adaptor.split(',')]
root_dir='%s'
generate_adapter='%s'
''' % (sample,sequence,adaptor,root_dir,generate_adapter)

if generate_adapter != 'n':
	code+='''
index='%s'
indexs=index.split(',')
	''' % (index)

code+='''
assert not os.system('mkdir QC')
os.chdir('QC')
for i,eachsample in enumerate(samples):
	assert not os.system('mkdir '+eachsample)
	os.chdir(eachsample)
	assert not os.system('ln -s %s %s' % (sequences[i][0],eachsample+'_1.fq.gz'))
	assert not os.system('ln -s %s %s' % (sequences[i][1],eachsample+'_2.fq.gz'))
	par_map=root_dir+'/libraryID'
	par_a=eachsample+'_1.fq.gz'
	par_b=eachsample+'_2.fq.gz'
	par_n=samples[i]
	if (adaptors[i][0].strip() != '') and (adaptors[i][1].strip() != ''):
		par_ad1=adaptors[i][0]
		par_ad2=adaptors[i][1]
		f=open('qc_'+eachsample+'.sh','w')
		f.write('perl /PUBLIC/source/RNA/noRef/QC/runQC_noref_v2.2.pl -a %s -b %s -ad1 %s -ad2 %s -n %s -mapfile %s' % (par_a,par_b,par_ad1,par_ad2,par_n,par_map))
		f.close()
		assert not os.system('sh qc_%s.sh' % (eachsample))
	else:
		f=open('qc_'+eachsample+'.sh','w')
		f.write('perl /PUBLIC/source/RNA/noRef/QC/runQC_noref_v2.2.pl -a %s -b %s -n %s -mapfile %s' % (par_a,par_b,par_n,par_map))
		if generate_adapter !='n':
			par_index=indexs[i]
			f.write(' -m_ad %s -index %s' % (generate_adapter,par_index))
		f.close()
		assert not os.system('sh qc_%s.sh' % (eachsample))
	assert not os.system('qsub -V -cwd -l vf=2G -l p=3 runQC_noref.sh')
	os.chdir('..')
os.chdir(root_dir)
'''

open('NOREF_step1.py','w').write(code)
##############################################################################
#for QC_report
code='''
import os

sample='%s'
project='%s'
os.chdir('QC')
f=open('QC_report.sh','w')
f.write('perl /PUBLIC/source/RNA/noRef/QC/noRef_QC_report_library.pl -n %%s -p %%s -dir .\\n' %% (sample, project))
f.write('perl /PUBLIC/source/RNA/noRef/QC/bin/clean_data_sum.pl -n %%s\\n' %% (sample))
f.close()
assert not os.system('sh QC_report.sh')
''' % (sample, project)

open('NOREF_step1_QCreport.py','w').write(code)



####################################################################
#for trinity (sh->py noly, not run) root_dir
code=''
if set([1]).issubset(includes):
	code+='''
import os
import os.path
import sys
import argparse
import re

parser = argparse.ArgumentParser(description='run trinity and resource list required')
parser.add_argument('--vf',required=True,help='momery required for qsub')
parser.add_argument('--h',required=True,help="the fat node for qsub")
argv=vars(parser.parse_args())
vf=argv['vf']
h=argv['h']

root_dir='%s'
sample='%s'
samples=sample.split(',')
ss='%s'
jc='%s'
minkmercov='%s'
minglue='%s'
project='%s'

clean_left=[]
clean_right=[]
for eachsample in samples:
	each_clean_left=root_dir+'/QC/'+eachsample+'/clean_data/'+eachsample+'_1.clean.fq'
	assert os.path.isfile(each_clean_left)
	clean_left.append(each_clean_left)
	each_clean_right=root_dir+'/QC/'+eachsample+'/clean_data/'+eachsample+'_2.clean.fq'
	assert os.path.isfile(each_clean_right)
	clean_right.append(each_clean_right)

clean_left=','.join(clean_left)
clean_right=','.join(clean_right)
rRNA=root_dir+'/QC/'+project+'.QC_results/results/1dataTable/dataTable'
assert os.path.isfile(rRNA)
cluster_num=re.search(r'([1-9])',h).group(1)
P_mem='mem'+str(cluster_num)
q_mem=P_mem+'.q'
if cluster_num==str(3) or cluster_num==str(4):
	P_mem='mem3'
	q_mem='mem3.q'
assert not os.system('mkdir TRINITY')
os.chdir('TRINITY')
f=open('trinity.sh','w')
f.write('perl /PUBLIC/source/RNA/noRef/Trinity_CDS/runTrinity_cds_v2.pl -left %%s -right %%s -ss %%s -jc %%s -minkmercov %%s -minglue %%s -JM %%s -rRNA %%s' %% (clean_left,clean_right,ss,jc,minkmercov,minglue,vf,rRNA))
f.close()
assert not os.system('sh trinity.sh')
assert not os.system('/opt/gridengine/bin/lx26-amd64/qsub -cwd -l vf=%%s -l p=10 -P %%s -q %%s -V run_Trinity_cds.sh' %% (vf,P_mem,q_mem))

''' % (root_dir,sample,ss,jc,minkmercov,minglue,project)

open('NOREF_step2.py','w').write(code)



##########################################################################
#for RSEM and SSR and the first step of SNP and ANNOTATION(sh->py only, not run)  root_dir
code='''
import os
import os.path
import sys
'''

if set([1,2]).issubset(includes):
	code+='''
import argparse
parser = argparse.ArgumentParser(description='some in-time parameters for annotation')
parser.add_argument('--new',help='where to qsub the annotation script,(1 on the New cluster,0 on the old DELL )',choices=['0','1'],default='1')
parser.add_argument('--n',help="split sequence file into given number, only for blast NR",default='0')
argv=vars(parser.parse_args())
ibm=argv['new']
n=argv['n'].strip()
assert n.isdigit()
'''

code+='''
sample='%s'
samples=sample.split(',')
root_dir='%s'
project='%s'

''' % (sample,root_dir,project)

if set([1,2]).issubset(includes):
        code+='''
##ANNOTATION
matrix='%s'
codon='%s'
db='%s'

assert not os.system('mkdir ANNOTATION_ALL')
os.chdir('ANNOTATION_ALL')
assert os.path.isfile(root_dir+'/TRINITY/assembly_INFO/unigene.fasta')
par_seq=root_dir+'/TRINITY/assembly_INFO/unigene.fasta'
num = n if int(n) else int(round(int(os.popen('grep -c \\">\\" '+par_seq).read().rstrip())/20000.0))
par_out=root_dir+'/ANNOTATION_ALL'
f=open('annot.sh','w')
f.write('perl /PUBLIC/source/RNA/noRef/ANNOTATION/annotation_v4.pl -seq %%s -out %%s -p %%s -matrix %%s -codon %%s -n %%s -db %%s -new %%s' %% (par_seq,par_out,project,matrix,codon,num,db,ibm))
f.close()
assert not os.system('sh annot.sh')
os.chdir(root_dir)
''' % (matrix,codon,db)

if set([1,5]).issubset(includes):
	code+='''
##RSEM
ss='%s'
assert not os.system('mkdir RSEM')
os.chdir('RSEM')
assert os.path.isfile(root_dir+'/TRINITY/trinity.out/Trinity.fasta')
for eachsample in samples:
	assert not os.system('mkdir %%s' %% eachsample)
	os.chdir(eachsample)
	assert os.path.isfile(root_dir+'/QC/'+eachsample+'/clean_data/'+eachsample+'_1.clean.fq')
	assert os.path.isfile(root_dir+'/QC/'+eachsample+'/clean_data/'+eachsample+'_2.clean.fq')
	par_a=root_dir+'/TRINITY/trinity.out/Trinity.fasta'
	par_r=root_dir+'/QC/'+eachsample+'/clean_data/'+eachsample+'_1.clean.fq'+','+root_dir+'/QC/'+eachsample+'/clean_data/'+eachsample+'_2.clean.fq'
	runRSEM=open('RSEM_'+eachsample+'.sh','w')
	runRSEM.write('perl /PUBLIC/source/RNA/noRef/RSEM/runRSEM_v1.5.pl -a %%s -r %%s -n %%s -s %%s -ss %%s' %% (par_a,par_r,project,eachsample,ss))
	runRSEM.close()
	assert not os.system('qsub -V -cwd -l vf=4G -l p=2 RSEM_%%s.sh' %% eachsample)
	os.chdir('..')
os.chdir(root_dir)
''' % (ss)

if set([1,3]).issubset(includes):
        code+='''
##the first step of SNP
#length='%s'
#lengths=length.split(',')
assert not os.system('mkdir SNP')
os.chdir('SNP')
assert os.path.isfile(root_dir+'/TRINITY/assembly_INFO/unigene.fasta')
assert not os.system('mkdir bam')
par_bamdir=root_dir+'/SNP/bam'
par_D=root_dir+'/TRINITY/assembly_INFO/unigene.fasta'
par_dir=root_dir+'/SNP'
f=open('snp_step1.sh','w')
f.write('/PUBLIC/software/RNA/bowtie2-2.2.3/bowtie2-build %s %s\\n' % (par_D,par_D))

for i,eachsample in enumerate(samples):
        #assert not os.system('mkdir %s' % eachsample)
        #os.chdir(eachsample)

        par_a=root_dir+'/QC/'+eachsample+'/clean_data/'+eachsample+'_1.clean.fq'
        par_b=root_dir+'/QC/'+eachsample+'/clean_data/'+eachsample+'_2.clean.fq'
        par_n=eachsample
        #par_L=lengths[i]
        par_sam=eachsample+'.sam'
        par_bam='bam/'+eachsample+'.bam'
        #par_i=root_dir+'/TRINITY/assembly_INFO/geneINFO'
        #assert os.path.isfile(par_i)

        f.write('/PUBLIC/software/RNA/bowtie2-2.2.3/bowtie2 -p 8 -x %s -1 %s -2 %s | samtools view -Sb - >%s\\n' % (par_D,par_a,par_b,par_bam))
        #os.chdir('..')
f.write('perl /PUBLIC/source/RNA/noRef/GATKsnp/runGATK2_v1.pl -ref %s -t bam -i %s \\n' % (par_D,par_bamdir))
f.write('sh workflow.sh\\n')
f.close()
assert not os.system('qsub -V -cwd -l vf=12G -l p=8 snp_step1.sh')
os.chdir(root_dir)
'''

if set([1,4]).issubset(includes):
	code+='''
##SSR
assert not os.system('mkdir SSR')
os.chdir('SSR')
assert os.path.isfile(root_dir+'/TRINITY/assembly_INFO/unigene.fasta')
f=open('ssr.sh','w')
f.write('perl /PUBLIC/source/RNA/noRef/SSR_Primer/runSSR_primer_v2.pl -a %s ' % (root_dir+'/TRINITY/assembly_INFO/unigene.fasta'))
f.close()
assert not os.system('sh ssr.sh')
assert not os.system('qsub -V -cwd -l vf=1G -l p=1 runSSR_primer.sh')
os.chdir(root_dir)
'''
open('NOREF_step3.py','w').write(code)

#########################################################
#for different expression (sh->py only, not run) root_dir
code=''
if set([1,5,6]).issubset(includes):
	code+='''
import os
import os.path
import sys

assert not os.system('mkdir DIFF_EXP')
os.chdir('DIFF_EXP')
sample='%s'
samples=sample.split(',')
root_dir='%s'
group='%s'
groupname='%s'
compare='%s'
venn_cluster='%s'

readcount=[]
for eachsample in samples:
	temp='%%s/RSEM/%%s/%%s.RSEM.out/%%s.Readcount_FPKM.xls' %% (root_dir,eachsample,eachsample,eachsample)
	assert os.path.isfile(temp)
	readcount.append(temp)
readcount=','.join(readcount)
geneinfo=root_dir+'/TRINITY/assembly_INFO/geneINFO'
assert os.path.isfile(geneinfo)

f=open('diff_exp.sh','w')
f.write('perl /PUBLIC/source/RNA/noRef/Diff_analysis/DE_analysis/runDE_analysis_v2.pl -i %%s -c %%s -group %%s -groupname %%s -g %%s -venn_cluster %%s' %% (readcount,geneinfo,group,groupname,compare,venn_cluster))
f.close()
assert not os.system('sh diff_exp.sh')
assert not os.system('qsub -V -cwd -l vf=2G -l p=2 Diff_Analysis.sh')
''' % (sample,root_dir,group,groupname,compare,venn_cluster)

open('NOREF_step4.py','w').write(code)

###########################################################
#for runAnnotationTable.sh,CDS prediction from ANNOTATION_ALL 
code='''
import os
import os.path

sample='%s'
samples=sample.split(',')
root_dir='%s'
''' % (sample,root_dir)

if set([1,2]).issubset(includes):
	code+='''
##runAnnotationTable.sh
os.chdir('ANNOTATION_ALL')
assert os.path.isfile('runAnnotationTable.sh')
assert not os.system('qsub -V -cwd -l vf=1G -l p=2 runAnnotationTable.sh')
'''

#if set([1,3]).issubset(includes):
	#code+='''
#SNPanalysis_step2.sh
#os.chdir(root_dir)
#for eachsample in samples:
	#os.chdir(root_dir+'/SNP/'+eachsample)
	#assert not os.system('qsub -V -cwd -l vf=1G -l p=1 SNPanalysis_step2.sh')
#''' 
open('NOREF_step5.py','w').write(code)

##########################################################
#for the second step of SNP, KEGG, GO and RNA-seqQC
code='''
import os
import os.path

root_dir='%s'
project='%s'
sample='%s'
samples=sample.split(',')
''' % (root_dir,project,sample)

if set([1,2,5,7]).issubset(includes):
	code+='''
##RNA-seqQC
assert not os.system('mkdir RNA_SEQ_QC')
os.chdir('RNA_SEQ_QC')
par_b=[]
for eachsample in samples:
	assert os.path.isfile(root_dir+'/RSEM/'+eachsample+'/'+eachsample+'.RSEM.out/'+eachsample+'.transcript.bam')
	par_b.append(root_dir+'/RSEM/'+eachsample+'/'+eachsample+'.RSEM.out/'+eachsample+'.transcript.bam')
par_b=','.join(par_b)
par_s=sample
par_a=root_dir+'/ANNOTATION_ALL/'+project+'.swissprot.xml'
assert os.path.isfile(par_a)
par_gi=root_dir+'/TRINITY/assembly_INFO/geneINFO'
assert os.path.isfile(par_gi)
par_ti=root_dir+'/TRINITY/assembly_INFO/isoformINFO'
assert os.path.isfile(par_ti)

f=open('rnaseq_qc.sh','w')
f.write('perl /PUBLIC/source/RNA/noRef/Uni_saturation/uni_saturation_v1.5.pl -b %s -s %s -a %s -gi %s -ti %s' % (par_b,par_s,par_a,par_gi,par_ti))
f.close()
assert not os.system('sh rnaseq_qc.sh')
assert not os.system('qsub -V -cwd -l vf=1G -l p=1 step5_Uni_saturation.sh')
os.chdir(root_dir)
'''

if set ([1,2,3]).issubset(includes):
	code+='''
##SNP step2 
os.chdir('SNP')
f=open('run.snp_step2.sh','w')
f.write('perl /PUBLIC/source/RNA/noRef/GATKsnp/get_gtf.pl '+root_dir+'/ANNOTATION_ALL/CDSprediction/'+project+'.blast.cds.fasta >'+root_dir+'/SNP/ResultsQ30/blast.cds.gtf\\n')
f.write('perl /PUBLIC/source/RNA/noRef/GATKsnp/annotationSNP.pl --outdir '+root_dir+'/SNP/ResultsQ30/ --name SNPs_annotation --gff '+root_dir+'/SNP/ResultsQ30/blast.cds.gtf --snp '+root_dir+'/SNP/ResultsQ30/SNPs.xls --fa '+root_dir+'/TRINITY/assembly_INFO/unigene.fasta\\n')
f.write('perl /PUBLIC/source/RNA/noRef/GATKsnp/snp_annot.pl '+root_dir+'/SNP/ResultsQ30/SNPs.xls '+root_dir+'/SNP/ResultsQ30/SNPs_annotation.info.txt '+root_dir+'/SNP/ResultsQ30/SNPs_non_synonymous.xls '+root_dir+'/SNP/ResultsQ30/SNPs_non_synonymous.stat.xls\\n')
f.write('rm -rf '+root_dir+'/SNP/dedup\\n')
f.close()
assert not os.system('qsub -V -cwd -l vf=0.5G run.snp_step2.sh')
os.chdir(root_dir)
'''
if set([1,2,4]).issubset(includes):
	code+='''
##SSR step2
os.chdir(root_dir+'/SSR/')
assert not os.system('python /PUBLIC/source/RNA/noRef/SSR_Primer/bin/SSR_position.py --project %s --root_dir %s'%(project,root_dir))
assert not os.system('qsub -V -cwd -l vf=1G -l p=1 SSR_utr.sh')
os.chdir(root_dir)
'''
if set([1,2,6]).issubset(includes):
        code+='''
###check DIFFgene number
assert not os.system('perl /PUBLIC/source/RNA/noRef/Diff_analysis/DE_analysis/diffreport.pl -diffdir %s ' % (root_dir+'/DIFF_EXP/'))
'''

if set([1,2,5,6,9]).issubset(includes):
	code+='''
##GO
compare_name='%s'
compare_names=[each.split(':') for each in compare_name.split(',')]
assert not os.system('mkdir GO_ENRICHMENT')
os.chdir('GO_ENRICHMENT')

assert not os.system('mkdir ALL')
os.chdir('ALL')
for each in compare_names:
	temp='vs'.join(each)
	par_i=root_dir+'/DIFF_EXP/Diff_analysis.out/'+temp+'/'+temp+'.DEGlist.txt'
	assert os.path.isfile(par_i)
	par_goann=root_dir+'/ANNOTATION_ALL/'+project+'.go.txt'
	assert os.path.isfile(par_goann)
	par_length=root_dir+'/TRINITY/assembly_INFO/geneINFO'
	assert os.path.isfile(par_length)
	assert not os.system('mkdir %%s' %% temp)
	os.chdir(temp)
	par_o=os.getcwd()
	py=open('run_goseq_'+temp+'.sh','w')
	py.write('perl /PUBLIC/software/RNA/GOseq/goseq_graph_v3.pl -i %%s -goann %%s -o %%s -length %%s -p %%s' %% (par_i,par_goann,par_o,par_length,temp))
	py.close()
	assert not os.system('qsub -V -cwd -l vf=2G -l p=1 run_goseq_%%s.sh' %% temp)
	os.chdir('..')
	
os.chdir('..')
assert not os.system('mkdir UP')
os.chdir('UP')
for each in compare_names:
	temp='vs'.join(each)
	par_i=root_dir+'/DIFF_EXP/Diff_analysis.out/'+temp+'/'+temp+'.DEGlist_up.txt'
	assert os.path.isfile(par_i)
	par_goann=root_dir+'/ANNOTATION_ALL/'+project+'.go.txt'
	assert os.path.isfile(par_goann)
	par_length=root_dir+'/TRINITY/assembly_INFO/geneINFO'
	assert os.path.isfile(par_length)
	assert not os.system('mkdir %%s' %% temp)
	os.chdir(temp)
	par_o=os.getcwd()
	py=open('run_goseq_up_'+temp+'.sh','w')
	py.write('perl /PUBLIC/software/RNA/GOseq/goseq_graph_v3.pl -i %%s -goann %%s -o %%s -length %%s -p %%s' %% (par_i,par_goann,par_o,par_length,temp+'_up'))
	py.close()
	assert not os.system('qsub -V -cwd -l vf=2G -l p=1 run_goseq_up_%%s.sh' %% temp)
	os.chdir('..')

os.chdir('..')
assert not os.system('mkdir DOWN')
os.chdir('DOWN')
for each in compare_names:
	temp='vs'.join(each)
	par_i=root_dir+'/DIFF_EXP/Diff_analysis.out/'+temp+'/'+temp+'.DEGlist_down.txt'
	assert os.path.isfile(par_i)
	par_goann=root_dir+'/ANNOTATION_ALL/'+project+'.go.txt'
	assert os.path.isfile(par_goann)
	par_length=root_dir+'/TRINITY/assembly_INFO/geneINFO'
	assert os.path.isfile(par_length)
	assert not os.system('mkdir %%s' %% temp)
	os.chdir(temp)
	par_o=os.getcwd()
	py=open('run_goseq_down_'+temp+'.sh','w')
	py.write('perl /PUBLIC/software/RNA/GOseq/goseq_graph_v3.pl -i %%s -goann %%s -o %%s -length %%s -p %%s' %% (par_i,par_goann,par_o,par_length,temp+'_down'))
	py.close()
	assert not os.system('qsub -V -cwd -l vf=2G -l p=1 run_goseq_down_%%s.sh' %% temp)
	os.chdir('..')

os.chdir(root_dir)
''' % (compare_name)

if set([1,2,5,6,8]).issubset(includes):
	code+='''
##KEGG
compare_name='%s'
compare_names=[each.split(':') for each in compare_name.split(',')]
assert not os.system('mkdir KEGG')
os.chdir('KEGG')

assert not os.system('mkdir ALL')
os.chdir('ALL')
for each in compare_names:
	temp='vs'.join(each)
	assert not os.system('mkdir %%s' %% temp)
	os.chdir(temp)
	par_diff=root_dir+'/DIFF_EXP/Diff_analysis.out/'+temp+'/'+temp+'.DEGlist.txt'
	par_ko=root_dir+'/ANNOTATION_ALL/koID.annotation'
	assert os.path.isfile(par_ko)
	py=open('runKEGG_enrich_'+temp+'.sh','w')
	py.write('perl /PUBLIC/source/RNA/noRef/KEGG_enrichment/runKEGG_enrich.pl -diff %%s -ko %%s -g %%s' %% (par_diff,par_ko,temp))
	py.close()
	assert not os.system('qsub -V -cwd -l vf=1G -l p=1 runKEGG_enrich_%%s.sh' %% temp)
	os.chdir('..')
	
os.chdir('..')	
assert not os.system('mkdir UP')
os.chdir('UP')
for each in compare_names:
	temp='vs'.join(each)
	assert not os.system('mkdir %%s' %% temp)
	os.chdir(temp)
	par_diff=root_dir+'/DIFF_EXP/Diff_analysis.out/'+temp+'/'+temp+'.DEGlist_up.txt'
	par_ko=root_dir+'/ANNOTATION_ALL/koID.annotation'
	assert os.path.isfile(par_ko)
	py=open('runKEGG_enrich_up_'+temp+'.sh','w')
	py.write('perl /PUBLIC/source/RNA/noRef/KEGG_enrichment/runKEGG_enrich.pl -diff %%s -ko %%s -g %%s' %% (par_diff,par_ko,temp+'_up'))
	py.close()
	assert not os.system('qsub -V -cwd -l vf=1G -l p=1 runKEGG_enrich_up_%%s.sh' %% temp)
	os.chdir('..')

os.chdir('..')	
assert not os.system('mkdir DOWN')
os.chdir('DOWN')
for each in compare_names:
	temp='vs'.join(each)
	assert not os.system('mkdir %%s' %% temp)
	os.chdir(temp)
	par_diff=root_dir+'/DIFF_EXP/Diff_analysis.out/'+temp+'/'+temp+'.DEGlist_down.txt'
	par_ko=root_dir+'/ANNOTATION_ALL/koID.annotation'
	assert os.path.isfile(par_ko)
	py=open('runKEGG_enrich_down_'+temp+'.sh','w')
	py.write('perl /PUBLIC/source/RNA/noRef/KEGG_enrichment/runKEGG_enrich.pl -diff %%s -ko %%s -g %%s' %% (par_diff,par_ko,temp+'_down'))
	py.close()
	assert not os.system('qsub -V -cwd -l vf=1G -l p=1 runKEGG_enrich_down_%%s.sh' %% temp)
	os.chdir('..')
	
os.chdir(root_dir)
''' % (compare_name)

if set([1,2,5,6,10]).issubset(includes):
        code+='''
##PPI
compare_name='%s'
compare_names=[each.split(':') for each in compare_name.split(',')]
ppi_number='%s'
assert not os.system('mkdir PPI')
os.chdir('PPI')
assert not os.system('mkdir ALL')
os.chdir('ALL')
for each in compare_names:
        temp='vs'.join(each)
        assert not os.system('mkdir %%s' %% temp)
        os.chdir(temp)
        par_diff=root_dir+'/DIFF_EXP/Diff_analysis.out/'+temp+'/'+temp+'.DEGlist.txt'
        par_a=root_dir+'/TRINITY/assembly_INFO/unigene.fasta'
        py=open('runPPI_'+temp+'.sh','w')
        py.write('python /PUBLIC/source/RNA/guoyang_tools/BLASTX_TO_PPI_v4.py --entries %%s --species %%s --fa %%s --name %%s --output . \\n' %% (par_diff,ppi_number,par_a,temp))
        py.write('cp %%s %%s' %% (root_dir+'/PPI/ALL/'+temp+'/'+temp+'.ppi.txt',root_dir+'/PPI/ALL/'+temp+'/'+temp+'.PPI.txt'))
        py.close()
        assert not os.system('qsub -V -cwd -l vf=1G -l p=4 runPPI_%%s.sh' %% temp)
        os.chdir('..')

os.chdir('..')
assert not os.system('mkdir UP')
os.chdir('UP')
for each in compare_names:
        temp='vs'.join(each)
        assert not os.system('mkdir %%s' %% temp)
	name=temp+'_up'
        os.chdir(temp)
        par_diff=root_dir+'/DIFF_EXP/Diff_analysis.out/'+temp+'/'+temp+'.DEGlist_up.txt'
        par_a=root_dir+'/TRINITY/assembly_INFO/unigene.fasta'
        py=open('runPPI_up_'+temp+'.sh','w')
        py.write('python /PUBLIC/source/RNA/guoyang_tools/BLASTX_TO_PPI_v4.py --entries %%s --species %%s --fa %%s --name %%s --output . \\n' %% (par_diff,ppi_number,par_a,name))
        py.write('cp %%s %%s' %% (root_dir+'/PPI/UP/'+temp+'/'+name+'.ppi.txt',root_dir+'/PPI/UP/'+temp+'/'+temp+'_up.PPI.txt'))
        py.close()
        assert not os.system('qsub -V -cwd -l vf=1G -l p=4 runPPI_up_%%s.sh' %% temp)
        os.chdir('..')

os.chdir('..')
assert not os.system('mkdir DOWN')
os.chdir('DOWN')
for each in compare_names:
        temp='vs'.join(each)
        assert not os.system('mkdir %%s' %% temp)
	name=temp+'_down'
        os.chdir(temp)
        par_diff=root_dir+'/DIFF_EXP/Diff_analysis.out/'+temp+'/'+temp+'.DEGlist_down.txt'
        par_a=root_dir+'/TRINITY/assembly_INFO/unigene.fasta'
        py=open('runPPI_down_'+temp+'.sh','w')
        py.write('python /PUBLIC/source/RNA/guoyang_tools/BLASTX_TO_PPI_v4.py --entries %%s --species %%s --fa %%s --name %%s --output . \\n' %% (par_diff,ppi_number,par_a,name))
        py.write('cp %%s %%s' %% (root_dir+'/PPI/DOWN/'+temp+'/'+name+'.ppi.txt',root_dir+'/PPI/DOWN/'+temp+'/'+temp+'_down.PPI.txt'))
        py.close()
        assert not os.system('qsub -V -cwd -l vf=1G -l p=4 runPPI_down_%%s.sh' %% temp)
        os.chdir('..')

os.chdir(root_dir)
''' % (compare_name,ppi_number)

open('NOREF_step6.py','w').write(code)

####################################################
#for web Visualization for KEGG      due to the network problem, can not qsub this step 
code='''
import os
import os.path
root_dir='%s'
''' % (root_dir)

if set([1,2]).issubset(includes):
	code+='''
os.chdir(root_dir+'/ANNOTATION_ALL')
par_ko=root_dir+'/ANNOTATION_ALL/Annotation_KO.xls'
assert os.path.isfile(par_ko)
assert not os.system('python /PUBLIC/source/RNA/guoyang_tools/pathway_ko_annotation_parallel_tolerant.pyc %s' % par_ko)
os.chdir(root_dir)
'''


if set([1,5,6,8]).issubset(includes):
	code+='''
compare_name='%s'
compare_names=[each.split(':') for each in compare_name.split(',')]
for each in compare_names:
	temp='vs'.join(each)
	os.chdir(root_dir+'/KEGG/ALL/'+temp)
	par_table=temp+'.DEG_KEGG_pathway_enrichment_add.xls' 
	par_diff=root_dir+'/DIFF_EXP/Diff_analysis.out/'+temp+'/'+temp+'.DEG.xls'
	if os.path.isfile(par_table) and os.path.isfile(par_diff):
		print '###'+temp+'###'
		assert not os.system('python /PUBLIC/source/RNA/guoyang_tools/pathway_annotation_flow_parallel_annotationfault_tolerant.pyc --table %%s --diff %%s' %% (par_table,par_diff))
	else:
		print '!!!'+temp+'!!! file(s) unavailable...'
		
	os.chdir(root_dir+'/KEGG/UP/'+temp)
	par_table=temp+'_up.DEG_KEGG_pathway_enrichment_add.xls'
	par_diff=root_dir+'/DIFF_EXP/Diff_analysis.out/'+temp+'/'+temp+'.DEG_up.xls'
	if os.path.isfile(par_table) and os.path.isfile(par_diff):
		print '###'+temp+' up###'
		assert not os.system('python /PUBLIC/source/RNA/guoyang_tools/pathway_annotation_flow_parallel_annotationfault_tolerant.pyc --table %%s --diff %%s' %% (par_table,par_diff))
	else:
		print '!!!'+temp+' up!!! file(s) unavailable...'
	
	os.chdir(root_dir+'/KEGG/DOWN/'+temp)
	par_table=temp+'_down.DEG_KEGG_pathway_enrichment_add.xls'
	par_diff=root_dir+'/DIFF_EXP/Diff_analysis.out/'+temp+'/'+temp+'.DEG_down.xls'
	if os.path.isfile(par_table) and os.path.isfile(par_diff):
		print '###'+temp+' down###'
		assert not os.system('python /PUBLIC/source/RNA/guoyang_tools/pathway_annotation_flow_parallel_annotationfault_tolerant.pyc --table %%s --diff %%s' %% (par_table,par_diff))
	else:
		print '!!!'+temp+' down!!! file(s) unavailable...'
''' % (compare_name)
open('NOREF_step7.py','w').write(code)
####################################################
# for summary
code='''
import os
import os.path
import glob
import gzip
import sys
import re
reload(sys)
sys.setdefaultencoding('utf-8')
from django.template import Template, Context, loader
from django.conf import settings
settings.configure(DEBUG=True, TEMPLATE_DEBUG=True,TEMPLATE_DIRS=('/PUBLIC/source/RNA/noRef/Pipeline_noRef/report_v2',))

def fasta_trim(file,k,typ):
	fa=[]
	for eachLine in open(file):
		if eachLine.startswith('>'):
			k-=1
		if k == 0:
			break
		if typ == 1:
			fa.append(eachLine)
		else:
			fa.append(eachLine.strip()+'<br />')
	return fa


def len_stat(file):
	len1=0 #200-500bp
	len2=0 #500-1kbp
	len3=0 #1k-2k
	len4=0 #>2k
	for eachLine in open(file):
		if eachLine.strip() != '':
			temp=int(eachLine.split()[-1].strip())
			if 200 <= temp < 500:
				len1+=1
				continue
			if 500 <= temp < 1000:
				len2+=1
				continue
			if 1000 <= temp < 2000:
				len3+=1
				continue
			if temp >= 2000:
				len4+=1
	return '200-500bp: '+str(len1)+'\\n500-1Kbp: '+str(len2)+'\\n1K-2Kbp: '+str(len3)+'\\n>2Kbp: '+str(len4)+'\\n'

root_dir='%s'
sample='%s'
samples=sample.split(',')
sequence='%s'
project='%s'
ss='%s'
flag_uniform=True
flag_gc=True
flag_GC_order=0
sequences=[each.strip().split(':') for each in sequence.split(',')]

#results
	
assert not os.system('mkdir '+project+'_result')
assert not os.system('mkdir '+project+'_result/'+project+'_results')
os.chdir(project+'_result/'+project+'_results')

result_order=0

result_order+=1
assert not os.system('mkdir '+str(result_order)+'.RawData')
rd_order=str(result_order)
os.chdir(''+str(result_order)+'.RawData')
f=open(root_dir+'/QC/rd_md5.sh','w')
for i,eachsample in enumerate(samples):
	f_in=gzip.open(sequences[i][0],'rb')
	f_out=[f_in.readline() for num in range(20)]
	open(eachsample+'.sampled_reads.raw.fastq','w').writelines(f_out)
	f_in.close()
	f.write('cd %%s\\n' %% (root_dir+'/QC/'+eachsample))
	f.write('md5sum %%s\\n' %% (eachsample+'_1.fq.gz'))
	f.write('md5sum %%s\\n' %% (eachsample+'_2.fq.gz'))
#f.write('cd %%s\\n' %% (os.getcwd()))
#f.write('mv rd_md5.sh.o* RawData_MD5.txt\\n')
#f.write('rm rd_md5.*\\n')
f.close()
#assert not os.system('qsub -V -cwd -l vf=0.5G -l p=1 rd_md5.sh')
assert not os.system('cp /PUBLIC/source/RNA/noRef/readmefiles_v3/1.RawData.README.txt %%s' %% (root_dir+'/'+project+'_result/'+project+'_results/'+str(result_order)+'.RawData'))
os.chdir(root_dir+'/'+project+'_result/'+project+'_results')

result_order+=1
assert not os.system('mkdir '+str(result_order)+'.QC')
qc_order=str(result_order)
os.chdir(''+str(result_order)+'.QC')
flag_GC_order+=1
assert not os.system('mkdir '+str(result_order)+'.1.ErrorRate')
os.chdir(''+str(result_order)+'.1.ErrorRate')
for eachsample in samples:
	assert not os.system('cp %%s %%s' %% (root_dir+'/QC/'+eachsample+'/files/'+eachsample+'.Error.png',eachsample+'.error_rate_distribution.png'))
assert not os.system('cp /PUBLIC/source/RNA/noRef/readmefiles_v3/2.1.ErrorRate.README.txt %%s' %% (root_dir+'/'+project+'_result/'+project+'_results/'+str(result_order)+'.QC/'+str(result_order)+'.1.ErrorRate'))
os.chdir('..')
if flag_gc:
	flag_GC_order+=1
	assert not os.system('mkdir '+str(result_order)+'.2.GC')
	os.chdir(''+str(result_order)+'.2.GC')
	for eachsample in samples:
		assert not os.system('cp %%s %%s' %% (root_dir+'/QC/'+eachsample+'/files/'+eachsample+'.GC.png',eachsample+'.GC_content_distribution.png'))
	assert not os.system('cp /PUBLIC/source/RNA/noRef/readmefiles_v3/2.2.GC.README.txt %%s' %% (root_dir+'/'+project+'_result/'+project+'_results/'+str(result_order)+'.QC/'+str(result_order)+'.2.GC'))
	os.chdir('..')

flag_GC_order+=1
assert not os.system('mkdir '+str(result_order)+'.'+str(flag_GC_order)+'.ReadsClassification')
os.chdir(''+str(result_order)+'.'+str(flag_GC_order)+'.ReadsClassification')
for eachsample in samples:
	assert not os.system('cp %%s %%s' %% (root_dir+'/QC/'+eachsample+'/files/'+eachsample+'.pie3d.png',eachsample+'.raw_reads_classification.png'))
assert not os.system('cp /PUBLIC/source/RNA/noRef/readmefiles_v3/2.3.ReadsClassification.README.txt %%s' %% (root_dir+'/'+project+'_result/'+project+'_results/'+str(result_order)+'.QC/'+str(result_order)+'.'+str(flag_GC_order)+'.ReadsClassification/'+str(result_order)+'.'+str(flag_GC_order)+'.ReadsClassification.README.txt'))
os.chdir('..')

flag_GC_order+=1
assert not os.system('mkdir '+str(result_order)+'.'+str(flag_GC_order)+'.CleanData_QCsummary')
os.chdir(''+str(result_order)+'.'+str(flag_GC_order)+'.CleanData_QCsummary')
f=open(root_dir+'/QC/cd_md5.sh','w')
for eachsample in samples:
	assert not os.system('head -20 %%s > %%s' %% (root_dir+'/QC/'+eachsample+'/clean_data/'+eachsample+'_1.clean.fq',eachsample+'.sampled_reads.clean.fastq'))
	f.write('cd %%s\\n' %% (root_dir+'/QC/'+eachsample+'/clean_data/'))
	f.write('gzip -c %%s > %%s\\n' %% (eachsample+'_1.clean.fq',eachsample+'_1.clean.fq.gz'))
	f.write('md5sum %%s\\n' %% (eachsample+'_1.clean.fq.gz'))
	f.write('gzip -c %%s > %%s\\n' %% (eachsample+'_2.clean.fq',eachsample+'_2.clean.fq.gz'))
	f.write('md5sum %%s\\n' %% (eachsample+'_2.clean.fq.gz'))
f.write('cd %%s\\n' %% (root_dir+'/QC/'))
#f.write('perl /PUBLIC/source/RNA/noRef/QC/bin/clean_data_sum.pl -n %%s\\n' %% (sample))
assert not os.system('cp %%s %%s\\n' %% (root_dir+'/QC/CleanData_QCsummary.xls',root_dir+'/'+project+'_result/'+project+'_results/'+str(result_order)+'.QC/'+str(result_order)+'.'+str(flag_GC_order)+'.CleanData_QCsummary'))
#f.write('cd %%s\\n' %% (os.getcwd()))
#f.write('mv cd_md5.sh.o* CleanData_MD5.txt\\n')
#f.write('rm cd_md5.*\\n')
f.close()
#assert not os.system('qsub -V -cwd -l vf=1G -l p=1 cd_md5.sh')
assert not os.system('cp /PUBLIC/source/RNA/noRef/readmefiles_v3/2.4.CleanData_QCsummary.README.txt %%s' %% (root_dir+'/'+project+'_result/'+project+'_results/'+str(result_order)+'.QC/'+str(result_order)+'.'+str(flag_GC_order)+'.CleanData_QCsummary/'+str(result_order)+'.'+str(flag_GC_order)+'.CleanData_QCsummary.README.txt'))
os.chdir(root_dir+'/'+project+'_result/'+project+'_results')
''' % (root_dir,sample,sequence,project,ss)

if set([1]).issubset(includes):
	code+='''
result_order+=1
minkmercov='%s'
assert not os.system('mkdir '+str(result_order)+'.TranscriptomeAssembly')
os.chdir(''+str(result_order)+'.TranscriptomeAssembly')
assert not os.system('mkdir '+str(result_order)+'.1.AssembledTranscriptome')
os.chdir(''+str(result_order)+'.1.AssembledTranscriptome')
assert not os.system('cp %%s .' %% (root_dir+'/TRINITY/trinity.out/Trinity.fasta'))
assert not os.system('cp %%s .' %% (root_dir+'/TRINITY/assembly_INFO/unigene.fasta'))
open('Sampled_trinity.fasta','w').writelines(fasta_trim('Trinity.fasta',5,1))
open('Sampled_unigene.fasta','w').writelines(fasta_trim('unigene.fasta',5,1))
assert not os.system('cp /PUBLIC/source/RNA/noRef/readmefiles_v3/3.1.AssembledTranscriptome.README.txt %%s' %% (root_dir+'/'+project+'_result/'+project+'_results/'+str(result_order)+'.TranscriptomeAssembly/'+str(result_order)+'.1.AssembledTranscriptome'))
os.chdir('..')

assert not os.system('mkdir '+str(result_order)+'.2.AssemblyINFO')
os.chdir(''+str(result_order)+'.2.AssemblyINFO')
assert not os.system('cp %%s .' %% (root_dir+'/TRINITY/assembly_INFO/Transcript_length_distribution.png'))
assert not os.system('cp %%s .' %% (root_dir+'/TRINITY/assembly_INFO/Unigene_length_distribution.png'))
#open('Gene_length_stat.txt','w').write(len_stat(root_dir+'/TRINITY/assembly_INFO/geneINFO'))
#open('Transcript_length_stat.txt','w').write(len_stat(root_dir+'/TRINITY/assembly_INFO/isoformINFO'))
# os.system('cp %%s %%s' %% (root_dir+'/TRINITY/assembly_INFO/geneINFO','Gene_length_stat.txt'))
# os.system('cp %%s %%s' %% (root_dir+'/TRINITY/assembly_INFO/isoformINFO','Transcript_length_stat.txt'))
assert not os.system('cp %%s .' %% (root_dir+'/TRINITY/assembly_INFO/Unigene_length_stat.txt'))
assert not os.system('cp %%s .' %% (root_dir+'/TRINITY/assembly_INFO/Transcript_length_stat.txt'))
assert not os.system('cp %%s .' %% (root_dir+'/TRINITY/assembly_INFO/Transcript_bar.png'))
assert not os.system('cp %%s .' %% (root_dir+'/TRINITY/assembly_INFO/Transcript_bar.pdf'))
assert not os.system('cp %%s .' %% (root_dir+'/TRINITY/assembly_INFO/Unigene_bar.png'))
assert not os.system('cp %%s .' %% (root_dir+'/TRINITY/assembly_INFO/Unigene_bar.pdf'))
assert not os.system('cp %%s .' %% (root_dir+'/TRINITY/assembly_INFO/Transcript_n50.txt'))
assert not os.system('cp %%s .' %% (root_dir+'/TRINITY/assembly_INFO/Unigene_n50.txt'))
assert not os.system('cp %%s .' %% (root_dir+'/TRINITY/assembly_INFO/Transcript_length.txt'))
assert not os.system('cp %%s .' %% (root_dir+'/TRINITY/assembly_INFO/Unigene_length.txt'))
assert not os.system('cp %%s .' %% (root_dir+'/TRINITY/assembly_INFO/Unigene_Transcript.length_distribution.png'))
assert not os.system('cp %%s .' %% (root_dir+'/TRINITY/assembly_INFO/Unigene_Transcript.length_distribution.pdf'))
assert not os.system('cp /PUBLIC/source/RNA/noRef/readmefiles_v3/3.2.AssemblyINFO.README.txt %%s' %% (root_dir+'/'+project+'_result/'+project+'_results/'+str(result_order)+'.TranscriptomeAssembly/'+str(result_order)+'.2.AssemblyINFO'))
os.chdir(root_dir+'/'+project+'_result/'+project+'_results')
''' % (minkmercov)

if set([1,2]).issubset(includes):
	code+='''
result_order+=1
assert not os.system('mkdir '+str(result_order)+'.GeneFunctionalAnnotation')
os.chdir(''+str(result_order)+'.GeneFunctionalAnnotation')
assert not os.system('mkdir '+str(result_order)+'.1.GeneFunctionalAnnotation')
os.chdir(''+str(result_order)+'.1.GeneFunctionalAnnotation')
assert not os.system('cp %s .' % (root_dir+'/ANNOTATION_ALL/Annotation_Summary.xls'))
assert not os.system('cp %s .' % (root_dir+'/ANNOTATION_ALL/Annotation_GO.xls'))
assert not os.system('cp %s .' % (root_dir+'/ANNOTATION_ALL/Annotation_KOG.xls'))
assert not os.system('cp %s .' % (root_dir+'/ANNOTATION_ALL/Annotation_Pfam.xls'))
assert not os.system('cp %s .' % (root_dir+'/ANNOTATION_ALL/Annotation_KO.xls'))
assert not os.system('cp %s .' % (root_dir+'/ANNOTATION_ALL/Blast_NR.xls'))
assert not os.system('cp %s .' % (root_dir+'/ANNOTATION_ALL/Blast_NT.xls'))
assert not os.system('cp %s .' % (root_dir+'/ANNOTATION_ALL/Blast_Swissprot.xls'))
assert not os.system('cp %s .' % (root_dir+'/ANNOTATION_ALL/Evalue_distribution.txt'))
assert not os.system('cp %s .' % (root_dir+'/ANNOTATION_ALL/Evalue_Distribution.3dpie.pdf'))
assert not os.system('cp %s .' % (root_dir+'/ANNOTATION_ALL/Evalue_Distribution.3dpie.png'))
assert not os.system('cp %s .' % (root_dir+'/ANNOTATION_ALL/Similarity_distribution.txt'))
assert not os.system('cp %s .' % (root_dir+'/ANNOTATION_ALL/Similarity_Distribution.3dpie.png'))
assert not os.system('cp %s .' % (root_dir+'/ANNOTATION_ALL/Similarity_Distribution.3dpie.pdf'))
assert not os.system('cp %s .' % (root_dir+'/ANNOTATION_ALL/Species_classification.txt'))
assert not os.system('cp %s .' % (root_dir+'/ANNOTATION_ALL/Species_classification.3dpie.pdf'))
assert not os.system('cp %s .' % (root_dir+'/ANNOTATION_ALL/Species_classification.3dpie.png'))
#assert not os.system('cp %s .' % (root_dir+'/ANNOTATION_ALL/Evalue_distribution.png'))
#assert not os.system('cp %s .' % (root_dir+'/ANNOTATION_ALL/Evalue_distribution.pdf'))
#assert not os.system('cp %s .' % (root_dir+'/ANNOTATION_ALL/Evalue_distribution.txt'))
#assert not os.system('cp %s .' % (root_dir+'/ANNOTATION_ALL/Similarity_distribution.txt'))
#assert not os.system('cp %s .' % (root_dir+'/ANNOTATION_ALL/Similarity_distribution.png'))
#assert not os.system('cp %s .' % (root_dir+'/ANNOTATION_ALL/Similarity_distribution.pdf'))
#assert not os.system('cp %s .' % (root_dir+'/ANNOTATION_ALL/go2gene.dat'))
assert not os.system('cp %s .' % (root_dir+'/ANNOTATION_ALL/complete_go.txt'))
assert not os.system('cp %s ./Annotation.statistics.xls' % (root_dir+'/ANNOTATION_ALL/Annotation.statistics'))

assert not os.system('mkdir Annotation_KO')
assert not os.system('cp %s Annotation_KO' % (root_dir+'/ANNOTATION_ALL/Annotation_KO.html'))
assert not os.system('cp -r %s Annotation_KO' % (root_dir+'/ANNOTATION_ALL/src'))
assert not os.system('cp /PUBLIC/source/RNA/noRef/readmefiles_v3/4.1.GeneFunctionalAnnotation.README.txt %s' % (root_dir+'/'+project+'_result/'+project+'_results/'+str(result_order)+'.GeneFunctionalAnnotation/'+str(result_order)+'.1.GeneFunctionalAnnotation/'+str(result_order)+'.1.GeneFunctionalAnnotation.README.txt'))
os.chdir('..')

assert not os.system('mkdir '+str(result_order)+'.2.GOclassification')
os.chdir(''+str(result_order)+'.2.GOclassification')
assert not os.system('cp %s .' % (root_dir+'/ANNOTATION_ALL/GO_classification_bar.png'))
assert not os.system('cp %s .' % (root_dir+'/ANNOTATION_ALL/GO_classification_bar.pdf'))
assert not os.system('cp %s .' % (root_dir+'/ANNOTATION_ALL/GO_classification.xls'))
assert not os.system('cp %s .' % (root_dir+'/ANNOTATION_ALL/GO_classification_count.txt'))
assert not os.system('cp %s .' % (root_dir+'/ANNOTATION_ALL/GO_classification.dat'))
assert not os.system('cp /PUBLIC/source/RNA/noRef/readmefiles_v3/4.2.GOclassification.README.txt %s' % (root_dir+'/'+project+'_result/'+project+'_results/'+str(result_order)+'.GeneFunctionalAnnotation/'+str(result_order)+'.2.GOclassification/'+str(result_order)+'.2.GOclassification.README.txt'))
os.chdir('..')

assert not os.system('mkdir '+str(result_order)+'.3.KOGclassification')
os.chdir(''+str(result_order)+'.3.KOGclassification')
assert not os.system('cp %s .' % (root_dir+'/ANNOTATION_ALL/KOG_classification.png'))
assert not os.system('cp %s .' % (root_dir+'/ANNOTATION_ALL/KOG_classification.pdf'))
assert not os.system('cp %s .' % (root_dir+'/ANNOTATION_ALL/KOG_classification.xls'))
assert not os.system('cp %s .' % (root_dir+'/ANNOTATION_ALL/KOG_classification_count.txt'))
assert not os.system('cp /PUBLIC/source/RNA/noRef/readmefiles_v3/4.3.KOGclassification.README.txt %s' % (root_dir+'/'+project+'_result/'+project+'_results/'+str(result_order)+'.GeneFunctionalAnnotation/'+str(result_order)+'.3.KOGclassification/'+str(result_order)+'.3.KOGclassification.README.txt'))
os.chdir('..')

assert not os.system('mkdir '+str(result_order)+'.4.KEGGclassification')
os.chdir(''+str(result_order)+'.4.KEGGclassification')
assert not os.system('cp %s .' % (root_dir+'/ANNOTATION_ALL/KEGG_classification.png'))
assert not os.system('cp %s .' % (root_dir+'/ANNOTATION_ALL/KEGG_classification.pdf')) 
assert not os.system('cp %s .' % (root_dir+'/ANNOTATION_ALL/KEGG_classification.xls'))
assert not os.system('cp %s .' % (root_dir+'/ANNOTATION_ALL/KEGG_classification_count.txt'))
assert not os.system('cp /PUBLIC/source/RNA/noRef/readmefiles_v3/4.4.KEGGclassification.README.txt %s' % (root_dir+'/'+project+'_result/'+project+'_results/'+str(result_order)+'.GeneFunctionalAnnotation/'+str(result_order)+'.4.KEGGclassification'))
os.chdir(root_dir+'/'+project+'_result/'+project+'_results')

result_order+=1
assert not os.system('mkdir '+str(result_order)+'.CDSprediction')
os.chdir(''+str(result_order)+'.CDSprediction')
assert not os.system('cp %s %s' % (root_dir+'/ANNOTATION_ALL/CDSprediction/'+project+'.blast.cds.fasta','unigene.blast.cds.fasta'))
assert not os.system('cp %s %s' % (root_dir+'/ANNOTATION_ALL/CDSprediction/'+project+'.blast.cds.length.txt','unigene.blast.cds.length.txt'))
assert not os.system('cp %s %s' % (root_dir+'/ANNOTATION_ALL/CDSprediction/'+project+'.blast.cds.pdf','unigene.blast.cds.pdf'))
assert not os.system('cp %s %s' % (root_dir+'/ANNOTATION_ALL/CDSprediction/'+project+'.blast.cds.png','unigene.blast.cds.png'))
assert not os.system('cp %s %s' % (root_dir+'/ANNOTATION_ALL/CDSprediction/'+project+'.blast.pep.fasta','unigene.blast.pep.fasta'))
assert not os.system('cp %s %s' % (root_dir+'/ANNOTATION_ALL/CDSprediction/'+project+'.blast.pep.length.txt','unigene.blast.pep.length.txt'))
assert not os.system('cp %s %s' % (root_dir+'/ANNOTATION_ALL/CDSprediction/'+project+'.blast.pep.pdf','unigene.blast.pep.pdf'))
assert not os.system('cp %s %s' % (root_dir+'/ANNOTATION_ALL/CDSprediction/'+project+'.blast.pep.png','unigene.blast.pep.png'))
assert not os.system('cp %s %s' % (root_dir+'/ANNOTATION_ALL/CDSprediction/'+project+'.estscan.cds.fasta','unigene.estscan.cds.fasta'))
assert not os.system('cp %s %s' % (root_dir+'/ANNOTATION_ALL/CDSprediction/'+project+'.estscan.cds.length.txt','unigene.estscan.cds.length.txt'))
assert not os.system('cp %s %s' % (root_dir+'/ANNOTATION_ALL/CDSprediction/'+project+'.estscan.cds.pdf','unigene.estscan.cds.pdf'))
assert not os.system('cp %s %s' % (root_dir+'/ANNOTATION_ALL/CDSprediction/'+project+'.estscan.cds.png','unigene.estscan.cds.png'))
assert not os.system('cp %s %s' % (root_dir+'/ANNOTATION_ALL/CDSprediction/'+project+'.estscan.pep.fasta','unigene.estscan.pep.fasta'))
assert not os.system('cp %s %s' % (root_dir+'/ANNOTATION_ALL/CDSprediction/'+project+'.estscan.pep.length.txt','unigene.estscan.pep.length.txt'))
assert not os.system('cp %s %s' % (root_dir+'/ANNOTATION_ALL/CDSprediction/'+project+'.estscan.pep.pdf','unigene.estscan.pep.pdf'))
assert not os.system('cp %s %s' % (root_dir+'/ANNOTATION_ALL/CDSprediction/'+project+'.estscan.pep.png','unigene.estscan.pep.png'))
assert not os.system('cp /PUBLIC/source/RNA/noRef/readmefiles_v3/5.CDSprediction.README.txt %s' % (root_dir+'/'+project+'_result/'+project+'_results/'+str(result_order)+'.CDSprediction/'+str(result_order)+'.CDSprediction.README.txt'))
os.chdir(root_dir+'/'+project+'_result/'+project+'_results')
'''

if set([1,3]).issubset(includes):
	code+='''
result_order+=1
assert not os.system('mkdir '+str(result_order)+'.SNPcalling')
os.chdir(''+str(result_order)+'.SNPcalling')
assert not os.system('cp %s .' % (root_dir+'/SNP/ResultsQ30/SNPs.xls'))
assert not os.system('cp %s .' % (root_dir+'/SNP/ResultsQ30/InDels.xls'))
assert not os.system('cp %s .' % (root_dir+'/SNP/ResultsQ30/SNPs_non_synonymous.xls'))
assert not os.system('cp %s .'% (root_dir+'/SNP/ResultsQ30/SNPs_non_synonymous.stat.xls'))
assert not os.system('cp /PUBLIC/source/RNA/noRef/readmefiles_v3/6.SNPcalling.README.txt %s' % (root_dir+'/'+project+'_result/'+project+'_results/'+str(result_order)+'.SNPcalling/'+str(result_order)+'.SNPcalling.README.txt'))
os.chdir(root_dir+'/'+project+'_result/'+project+'_results')
'''

if set([1,4]).issubset(includes):
	code+='''
result_order+=1
assert not os.system('mkdir '+str(result_order)+'.SSRdetection')
os.chdir(''+str(result_order)+'.SSRdetection')
assert not os.system('cp %s %s' % (root_dir+'/SSR/SSR_Primer/'+project+'.SSR.pos.misa','unigene.fasta.misa'))
assert not os.system('cp %s .' % (root_dir+'/SSR/SSR_Primer/unigene.fasta.SSR_primer'))
assert not os.system('cp %s .' % (root_dir+'/SSR/SSR_Primer/unigene.fasta.SSR_motifs_distribution.png'))
assert not os.system('cp %s .' % (root_dir+'/SSR/SSR_Primer/unigene.fasta.SSR_motifs_distribution.pdf'))
assert not os.system('cp %s %s' % (root_dir+'/SSR/SSR_Primer/unigene.fasta.ssr_classification.txt','unigene.fasta.SSR_classification.txt '))
assert not os.system('cp %s %s' % (root_dir+'/SSR/SSR_Primer/unigene.fasta.ssr_summary.txt','unigene.fasta.SSR_summary.txt'))
assert not os.system('cp %s %s' % (root_dir+'/SSR/SSR_Primer/unigene.fasta.ssr_motif_frequency.txt','unigene.fasta.SSR_motif_frequency.txt'))
assert not os.system('cp %s .' % (root_dir+'/SSR/SSR_Primer/unigene.fasta.repeat_type_frequency.txt'))
assert not os.system('cp /PUBLIC/source/RNA/noRef/readmefiles_v3/7.newSSRdetection.README.txt %s' % (root_dir+'/'+project+'_result/'+project+'_results/'+str(result_order)+'.SSRdetection/'+str(result_order)+'.SSRdetection.README.txt'))
os.chdir(root_dir+'/'+project+'_result/'+project+'_results')
'''

if set([1,5]).issubset(includes):
	code+='''
result_order+=1
assert not os.system('mkdir '+str(result_order)+'.GeneExprQuantification')
RSEM_order=result_order
os.chdir(''+str(result_order)+'.GeneExprQuantification')
for eachsample in samples:
	assert not os.system('cp %s .' % (root_dir+'/RSEM/'+eachsample+'/'+eachsample+'.RSEM.out/'+eachsample+'.Readcount_FPKM.xls'))
	assert not os.system('cp %s .' % (root_dir+'/RSEM/'+eachsample+'/'+eachsample+'.RSEM.out/'+eachsample+'.FPKM_density_distribution.pdf'))
	assert not os.system('cp %s .' % (root_dir+'/RSEM/'+eachsample+'/'+eachsample+'.RSEM.out/'+eachsample+'.FPKM_density_distribution.png'))
assert not os.system('cp /PUBLIC/source/RNA/noRef/readmefiles_v3/8.GeneExprQuantification.README.txt %s' % (root_dir+'/'+project+'_result/'+project+'_results/'+str(result_order)+'.GeneExprQuantification/'+str(result_order)+'.GeneExprQuantification.README.txt'))
os.chdir(root_dir+'/'+project+'_result/'+project+'_results')
GeneExprQuantification_order=result_order
'''

if set([1,2,5]).issubset(includes):
	code+='''
os.chdir(str(GeneExprQuantification_order)+'.GeneExprQuantification')
for eachsample in samples:
	assert not os.system('perl /PUBLIC/source/RNA/noRef/sample_FPKM_annot.pl %s %s %s' % (eachsample+'.Readcount_FPKM.xls',root_dir+'/ANNOTATION_ALL/Annotation_Summary.xls',eachsample+'_annotation.xls'))
if len(samples) > 1:
	assert not os.system('perl /PUBLIC/source/RNA/noRef/sample_only_total_annot_v1.pl -n %s -a %s ' % (sample,root_dir+'/ANNOTATION_ALL/Annotation_Summary.xls'))
os.chdir(root_dir+'/'+project+'_result/'+project+'_results')
''' 

if set([1,2,5,7]).issubset(includes):
	code+='''
result_order+=1
assert not os.system('mkdir '+str(result_order)+'.RNA-seqQC')
os.chdir(''+str(result_order)+'.RNA-seqQC')
assert not os.system('mkdir '+str(result_order)+'.1.SaturationCurve')
os.chdir(''+str(result_order)+'.1.SaturationCurve')
for eachsample in samples:
	assert not os.system('cp %s .' % (root_dir+'/RNA_SEQ_QC/RNA_seq_QC/'+eachsample+'.Saturation_curve.png'))
	assert not os.system('cp %s .' % (root_dir+'/RNA_SEQ_QC/RNA_seq_QC/'+eachsample+'.Saturation_curve.pdf'))
	assert not os.system('cp %s .' % (root_dir+'/RNA_SEQ_QC/RNA_seq_QC/'+eachsample+'.Saturation_curve.txt'))
assert not os.system('cp /PUBLIC/source/RNA/noRef/readmefiles_v3/9.1.SaturationCurve.README.txt %s' % (root_dir+'/'+project+'_result/'+project+'_results/'+str(result_order)+'.RNA-seqQC/'+str(result_order)+'.1.SaturationCurve/'+str(result_order)+'.1.SaturationCurve.README.txt'))
os.chdir('..')

if flag_uniform:
	assert not os.system('mkdir '+str(result_order)+'.2.UniformDistribution')
	os.chdir(''+str(result_order)+'.2.UniformDistribution')
	for eachsample in samples:
		assert not os.system('cp %s .' % (root_dir+'/RNA_SEQ_QC/RNA_seq_QC/'+eachsample+'.Uniform_distribution.png'))
		assert not os.system('cp %s .' % (root_dir+'/RNA_SEQ_QC/RNA_seq_QC/'+eachsample+'.Uniform_distribution.pdf'))
		assert not os.system('cp %s .' % (root_dir+'/RNA_SEQ_QC/RNA_seq_QC/'+eachsample+'.Uniform_distribution.txt'))
	assert not os.system('cp /PUBLIC/source/RNA/noRef/readmefiles_v3/9.2.UniformDistribution.README.txt %s' % (root_dir+'/'+project+'_result/'+project+'_results/'+str(result_order)+'.RNA-seqQC/'+str(result_order)+'.2.UniformDistribution/'+str(result_order)+'.2.UniformDistribution.README.txt'))
os.chdir(root_dir+'/'+project+'_result/'+project+'_results')
rnaseq_qc_order=result_order
'''

if set([1,2,5,6,7]).issubset(includes):
	code+='''
if flag_uniform:
	os.chdir(str(rnaseq_qc_order)+'.RNA-seqQC')
	assert not os.system('mkdir '+str(rnaseq_qc_order)+'.3.Correlation')
	os.chdir(''+str(rnaseq_qc_order)+'.3.Correlation')
	for png in glob.glob(root_dir+'/DIFF_EXP/Diff_analysis.out/corr_plot/*.png'):
		assert not os.system('cp %s .' % png)
	for pdf in glob.iglob(root_dir+'/DIFF_EXP/Diff_analysis.out/corr_plot/*.pdf'):
		assert not os.system('cp %s .' % pdf)
	assert not os.system('cp /PUBLIC/source/RNA/noRef/readmefiles_v3/9.3.Correlation.README.txt %s' % (root_dir+'/'+project+'_result/'+project+'_results/'+str(rnaseq_qc_order)+'.RNA-seqQC/'+str(rnaseq_qc_order)+'.3.Correlation/'+str(result_order)+'.3.Correlation.README.txt'))
	os.chdir(root_dir+'/'+project+'_result/'+project+'_results')
else:
	os.chdir(str(rnaseq_qc_order)+'.RNA-seqQC')
	assert not os.system('mkdir '+str(rnaseq_qc_order)+'.2.Correlation')
	os.chdir(''+str(rnaseq_qc_order)+'.2.Correlation')
	for png in glob.glob(root_dir+'/DIFF_EXP/Diff_analysis.out/corr_plot/*.png'):
		assert not os.system('cp %s .' % png)
	for pdf in glob.iglob(root_dir+'/DIFF_EXP/Diff_analysis.out/corr_plot/*.pdf'):
		assert not os.system('cp %s .' % pdf)
	assert not os.system('cp /PUBLIC/source/RNA/noRef/readmefiles_v3/9.3.Correlation.README.txt %s' % (root_dir+'/'+project+'_result/'+project+'_results/'+str(rnaseq_qc_order)+'.RNA-seqQC/'+str(rnaseq_qc_order)+'.2.Correlation/'+str(result_order)+'.2.Correlation.README.txt'))
	os.chdir(root_dir+'/'+project+'_result/'+project+'_results')
'''

if set([1,5,6]).issubset(includes):
	code+='''
compare_name='%s'
venn_cluster_name='%s'
group_iter='%s'
flag_repeat=%s
groups_iter=group_iter.split(',')
compare_names=[each.split(':') for each in compare_name.split(',')]
venn_cluster_vs_names=[]
flag_venn=False
for each in venn_cluster_name.split(','):
	if each.count(':') >= 2 and each.count(':') <= 5:
		flag_venn=True
		venn_cluster_vs_names.append(each.replace(':','vs'))

os.chdir(str(GeneExprQuantification_order)+'.GeneExprQuantification')
assert not os.system('cp %%s .' %% (root_dir+'/DIFF_EXP/Diff_analysis.out/FPKM_interval.xls'))
os.chdir(root_dir+'/'+project+'_result/'+project+'_results')

result_order+=1
assert not os.system('mkdir '+str(result_order)+'.DiffExprAnalysis')
os.chdir(''+str(result_order)+'.DiffExprAnalysis')
assert not os.system('mkdir '+str(result_order)+'.1.GeneExpContrast')
os.chdir(''+str(result_order)+'.1.GeneExpContrast')
assert not os.system('cp %%s .' %% (root_dir+'/DIFF_EXP/Diff_analysis.out/FPKM_density_distribution.png'))
assert not os.system('cp %%s .' %% (root_dir+'/DIFF_EXP/Diff_analysis.out/FPKM_density_distribution.pdf'))
assert not os.system('cp %%s .' %% (root_dir+'/DIFF_EXP/Diff_analysis.out/FPKM_boxplot.png'))
assert not os.system('cp %%s .' %% (root_dir+'/DIFF_EXP/Diff_analysis.out/FPKM_boxplot.pdf'))
assert not os.system('cp /PUBLIC/source/RNA/noRef/readmefiles_v3/10.1.GeneExpContrast.README.txt %%s' %% (root_dir+'/'+project+'_result/'+project+'_results/'+str(result_order)+'.DiffExprAnalysis/'+str(result_order)+'.1.GeneExpContrast/'+str(result_order)+'.1.GeneExpContrast.README.txt'))
os.chdir('..')

assert not os.system('mkdir '+str(result_order)+'.2.DiffExprAnalysis')
os.chdir(''+str(result_order)+'.2.DiffExprAnalysis')
for each in compare_names:
	temp='vs'.join(each)
	assert not os.system('cp %%s .' %% (root_dir+'/DIFF_EXP/Diff_analysis.out/'+temp+'/'+temp+'.Differential_analysis_results.xls'))
	assert not os.system('cp %%s .' %% (root_dir+'/DIFF_EXP/Diff_analysis.out/'+temp+'/'+temp+'.DEG.xls'))
	assert not os.system('cp %%s .' %% (root_dir+'/DIFF_EXP/Diff_analysis.out/'+temp+'/'+temp+'.DEGlist.txt'))
	assert not os.system('cp %%s .' %% (root_dir+'/DIFF_EXP/Diff_analysis.out/'+temp+'/'+temp+'.DEGlist_up.txt'))
	assert not os.system('cp %%s .' %% (root_dir+'/DIFF_EXP/Diff_analysis.out/'+temp+'/'+temp+'.DEGlist_down.txt'))
	assert not os.system('cp %%s .' %% (root_dir+'/DIFF_EXP/Diff_analysis.out/'+temp+'/'+temp+'.DEG_up.xls'))
	assert not os.system('cp %%s .' %% (root_dir+'/DIFF_EXP/Diff_analysis.out/'+temp+'/'+temp+'.DEG_down.xls'))
assert not os.system('cp /PUBLIC/source/RNA/noRef/readmefiles_v3/10.2.DiffExprAnalysis.README.txt %%s' %% (root_dir+'/'+project+'_result/'+project+'_results/'+str(result_order)+'.DiffExprAnalysis/'+str(result_order)+'.2.DiffExprAnalysis/'+str(result_order)+'.2.DiffExprAnalysis.README.txt'))
os.chdir('..')

assert not os.system('mkdir '+str(result_order)+'.3.DEGsFilter')
os.chdir(''+str(result_order)+'.3.DEGsFilter')
for each in compare_names:
	temp='vs'.join(each)
	assert not os.system('cp %%s .' %% (root_dir+'/DIFF_EXP/Diff_analysis.out/'+temp+'/'+temp+'.Volcanoplot.pdf'))
	assert not os.system('cp %%s .' %% (root_dir+'/DIFF_EXP/Diff_analysis.out/'+temp+'/'+temp+'.Volcanoplot.png'))
assert not os.system('cp /PUBLIC/source/RNA/noRef/readmefiles_v3/10.3.DEGsFilter.README.txt %%s' %% (root_dir+'/'+project+'_result/'+project+'_results/'+str(result_order)+'.DiffExprAnalysis/'+str(result_order)+'.3.DEGsFilter/'+str(result_order)+'.3.DEGsFilter.README.txt'))
os.chdir('..')


if flag_venn:
	assert not os.system('mkdir '+str(result_order)+'.4.VennDiagram')
	os.chdir(''+str(result_order)+'.4.VennDiagram')
	for each in venn_cluster_vs_names:
		assert not os.system('cp -r %%s .' %% (root_dir+'/DIFF_EXP/Diff_analysis.out/'+each+'/'))
		os.chdir(each)
		txt=glob.glob('*.txt')
		for eachtxt in txt:
			assert not os.system('perl /PUBLIC/source/RNA/noRef/Diff_analysis/DE_analysis/venn_annot.pl %%s %%s' %% (eachtxt,root_dir+'/ANNOTATION_ALL/Annotation_Summary.xls'))
		os.chdir('..')

	assert not os.system('cp /PUBLIC/source/RNA/noRef/readmefiles_v3/10.4.VennDiagram.README.txt %%s' %% (root_dir+'/'+project+'_result/'+project+'_results/'+str(result_order)+'.DiffExprAnalysis/'+str(result_order)+'.4.VennDiagram/'+str(result_order)+'.4.VennDiagram.README.txt'))
	os.chdir('..')
	assert not os.system('mkdir '+str(result_order)+'.5.DEGcluster')
	os.chdir(''+str(result_order)+'.5.DEGcluster')
	assert not os.system('cp /PUBLIC/source/RNA/noRef/readmefiles_v3/10.5.DEGcluster.README.txt %%s' %% (root_dir+'/'+project+'_result/'+project+'_results/'+str(result_order)+'.DiffExprAnalysis/'+str(result_order)+'.5.DEGcluster/'+str(result_order)+'.5.DEGcluster.README.txt'))
	assert not os.system('cp %%s DEG_union_for_cluster.txt' %% (root_dir+'/DIFF_EXP/Diff_analysis.out/DEG_union_for_cluster'))
	assert not os.system('cp %%s .' %% (root_dir+'/DIFF_EXP/Diff_analysis.out/Hcluster_heatmap.png'))
	assert not os.system('cp %%s .' %% (root_dir+'/DIFF_EXP/Diff_analysis.out/Hcluster_heatmap.pdf'))
	os.system('cp %%s .' %% (root_dir+'/DIFF_EXP/Diff_analysis.out/Hcluster_heatmap.detail.pdf'))
	assert not os.system('cp -r %%s .' %% (root_dir+'/DIFF_EXP/Diff_analysis.out/H_cluster'))
	assert not os.system('cp -r %%s .' %% (root_dir+'/DIFF_EXP/Diff_analysis.out/SOM_cluster'))
	assert not os.system('cp -r %%s .' %% (root_dir+'/DIFF_EXP/Diff_analysis.out/K_means_cluster'))
	os.chdir(root_dir+'/'+project+'_result/'+project+'_results')
	diff_exp_order=result_order

else:
	assert not os.system('mkdir '+str(result_order)+'.4.DEGcluster')
	os.chdir(''+str(result_order)+'.4.DEGcluster')
	assert not os.system('cp %%s DEG_union_for_cluster.txt' %% (root_dir+'/DIFF_EXP/Diff_analysis.out/DEG_union_for_cluster'))
	assert not os.system('cp %%s .' %% (root_dir+'/DIFF_EXP/Diff_analysis.out/Hcluster_heatmap.png'))
	assert not os.system('cp %%s .' %% (root_dir+'/DIFF_EXP/Diff_analysis.out/Hcluster_heatmap.pdf'))
	os.system('cp %%s .' %% (root_dir+'/DIFF_EXP/Diff_analysis.out/Hcluster_heatmap.detail.pdf'))
	assert not os.system('cp -r %%s .' %% (root_dir+'/DIFF_EXP/Diff_analysis.out/H_cluster'))
	assert not os.system('cp -r %%s .' %% (root_dir+'/DIFF_EXP/Diff_analysis.out/SOM_cluster'))
	assert not os.system('cp -r %%s .' %% (root_dir+'/DIFF_EXP/Diff_analysis.out/K_means_cluster'))
	assert not os.system('cp /PUBLIC/source/RNA/noRef/readmefiles_v3/10.5.DEGcluster.README.txt %%s' %% (root_dir+'/'+project+'_result/'+project+'_results/'+str(result_order)+'.DiffExprAnalysis/'+str(result_order)+'.4.DEGcluster/'+str(result_order)+'.4.DEGcluster.README.txt'))
	os.chdir(root_dir+'/'+project+'_result/'+project+'_results')
	diff_exp_order=result_order

''' % (compare_name,venn_cluster_name,group_iter,flag_repeat)

if set([1,2,5,6]).issubset(includes):
	code+='''
os.chdir(''+str(diff_exp_order)+'.DiffExprAnalysis')
if flag_venn:
	assert not os.system('mkdir '+str(diff_exp_order)+'.6.DEGannotation')
	os.chdir(''+str(diff_exp_order)+'.6.DEGannotation')
	ann_sum=open(root_dir+'/ANNOTATION_ALL/Annotation_Summary.xls').readlines()
	ann_tit=ann_sum.pop(0).strip().split('\t')
	ann={eachann.strip().split('\t')[0]:eachann.strip().split('\t')[1:] for eachann in ann_sum if eachann != ''}
	for each in compare_names:
		temp='vs'.join(each)
		result=[]
		diff_sum=open(root_dir+'/DIFF_EXP/Diff_analysis.out/'+temp+'/'+temp+'.DEG.xls').readlines()
		diff_tit=diff_sum.pop(0).strip()
		for eachLine in diff_sum:
			if (eachLine.strip() != '') and (eachLine.split()[0].strip() in ann):
				result.append(eachLine.strip()+'\\t'+'\\t'.join(ann[eachLine.split()[0].strip()])+'\\n')
		title=[diff_tit+'\\t'+'\\t'.join(ann_tit[1:])+'\\n']	
		open(temp+'.DEGannot.xls','w').writelines(title+result)		
	
		result=[]
		diff_up=open(root_dir+'/DIFF_EXP/Diff_analysis.out/'+temp+'/'+temp+'.DEG_up.xls').readlines()
		diff_tit=diff_up.pop(0).strip()
		for eachLine in diff_up:
			if (eachLine.strip() != '') and (eachLine.split()[0].strip() in ann):
				result.append(eachLine.strip()+'\\t'+'\\t'.join(ann[eachLine.split()[0].strip()])+'\\n')
		title=[diff_tit+'\\t'+'\\t'.join(ann_tit[1:])+'\\n']	
		open(temp+'.DEG_up_annot.xls','w').writelines(title+result)		

		result=[]
		diff_down=open(root_dir+'/DIFF_EXP/Diff_analysis.out/'+temp+'/'+temp+'.DEG_down.xls').readlines()
		diff_tit=diff_down.pop(0).strip()
		for eachLine in diff_down:
			if (eachLine.strip() != '') and (eachLine.split()[0].strip() in ann):
				result.append(eachLine.strip()+'\\t'+'\\t'.join(ann[eachLine.split()[0].strip()])+'\\n')
		title=[diff_tit+'\\t'+'\\t'.join(ann_tit[1:])+'\\n']	
		open(temp+'.DEG_down_annot.xls','w').writelines(title+result)		
	
	assert not os.system('cp /PUBLIC/source/RNA/noRef/readmefiles_v3/10.6.DEGannotation.README.txt %s' % (root_dir+'/'+project+'_result/'+project+'_results/'+str(diff_exp_order)+'.DiffExprAnalysis/'+str(diff_exp_order)+'.6.DEGannotation/'+str(result_order)+'.6.DEGannotation.README.txt'))
	os.chdir(root_dir+'/'+project+'_result/'+project+'_results')

else:
	assert not os.system('mkdir '+str(diff_exp_order)+'.5.DEGannotation')
	os.chdir(''+str(diff_exp_order)+'.5.DEGannotation')
	ann_sum=open(root_dir+'/ANNOTATION_ALL/Annotation_Summary.xls').readlines()
	ann_tit=ann_sum.pop(0).strip().split('\t')
	ann={eachann.strip().split('\t')[0]:eachann.strip().split('\t')[1:] for eachann in ann_sum if eachann != ''}
	for each in compare_names:
		temp='vs'.join(each)
		result=[]
		diff_sum=open(root_dir+'/DIFF_EXP/Diff_analysis.out/'+temp+'/'+temp+'.DEG.xls').readlines()
		diff_tit=diff_sum.pop(0).strip()
		for eachLine in diff_sum:
			if (eachLine.strip() != '') and (eachLine.split()[0].strip() in ann):
				result.append(eachLine.strip()+'\\t'+'\\t'.join(ann[eachLine.split()[0].strip()])+'\\n')
		title=[diff_tit+'\\t'+'\\t'.join(ann_tit[1:])+'\\n']	
		open(temp+'.DEGannot.xls','w').writelines(title+result)		

		result=[]
		diff_up=open(root_dir+'/DIFF_EXP/Diff_analysis.out/'+temp+'/'+temp+'.DEG_up.xls').readlines()
		diff_tit=diff_up.pop(0).strip()
		for eachLine in diff_up:
			if (eachLine.strip() != '') and (eachLine.split()[0].strip() in ann):
				result.append(eachLine.strip()+'\\t'+'\\t'.join(ann[eachLine.split()[0].strip()])+'\\n')
		title=[diff_tit+'\\t'+'\\t'.join(ann_tit[1:])+'\\n']	
		open(temp+'.DEG_up_annot.xls','w').writelines(title+result)		

		result=[]
		diff_down=open(root_dir+'/DIFF_EXP/Diff_analysis.out/'+temp+'/'+temp+'.DEG_down.xls').readlines()
		diff_tit=diff_down.pop(0).strip()
		for eachLine in diff_down:
			if (eachLine.strip() != '') and (eachLine.split()[0].strip() in ann):
				result.append(eachLine.strip()+'\\t'+'\\t'.join(ann[eachLine.split()[0].strip()])+'\\n')
		title=[diff_tit+'\\t'+'\\t'.join(ann_tit[1:])+'\\n']	
		open(temp+'.DEG_down_annot.xls','w').writelines(title+result)		
	
	assert not os.system('cp /PUBLIC/source/RNA/noRef/readmefiles_v3/10.6.DEGannotation.README.txt %s' % (root_dir+'/'+project+'_result/'+project+'_results/'+str(diff_exp_order)+'.DiffExprAnalysis/'+str(diff_exp_order)+'.5.DEGannotation/'+str(result_order)+'.5.DEGannotation.README.txt'))
	os.chdir(root_dir+'/'+project+'_result/'+project+'_results')


'''

if set([1,2,5,6,9]).issubset(includes):
	code+='''
result_order+=1
assert not os.system('mkdir '+str(result_order)+'.DEG_GOenrichment')
os.chdir(''+str(result_order)+'.DEG_GOenrichment')
assert not os.system('mkdir '+str(result_order)+'.1.GOenrichment')
os.chdir(''+str(result_order)+'.1.GOenrichment')
for each in compare_names:
	temp='vs'.join(each)
	os.system('cp %s .' % (root_dir+'/GO_ENRICHMENT/ALL/'+temp+'/'+temp+'.DEG_GO_enrichment_result.xls'))
	os.system('cp %s .' % (root_dir+'/GO_ENRICHMENT/ALL/'+temp+'/'+temp+'.DEG_Enriched_GO_classification.png'))
	os.system('cp %s .' % (root_dir+'/GO_ENRICHMENT/ALL/'+temp+'/'+temp+'.DEG_Enriched_GO_classification.pdf'))
	os.system('cp %s .' % (root_dir+'/GO_ENRICHMENT/ALL/'+temp+'/'+temp+'.DEG_Enriched_GO_classification_gene_count.txt'))
	os.system('cp %s .' % (root_dir+'/GO_ENRICHMENT/UP/'+temp+'/'+temp+'_up.DEG_GO_enrichment_result.xls'))
	os.system('cp %s .' % (root_dir+'/GO_ENRICHMENT/UP/'+temp+'/'+temp+'_up.DEG_Enriched_GO_classification.png'))
	os.system('cp %s .' % (root_dir+'/GO_ENRICHMENT/UP/'+temp+'/'+temp+'_up.DEG_Enriched_GO_classification.pdf'))
	os.system('cp %s .' % (root_dir+'/GO_ENRICHMENT/UP/'+temp+'/'+temp+'_up.DEG_Enriched_GO_classification_gene_count.txt'))
	os.system('cp %s .' % (root_dir+'/GO_ENRICHMENT/DOWN/'+temp+'/'+temp+'_down.DEG_GO_enrichment_result.xls'))
	os.system('cp %s .' % (root_dir+'/GO_ENRICHMENT/DOWN/'+temp+'/'+temp+'_down.DEG_Enriched_GO_classification.png'))
	os.system('cp %s .' % (root_dir+'/GO_ENRICHMENT/DOWN/'+temp+'/'+temp+'_down.DEG_Enriched_GO_classification.pdf'))
	os.system('cp %s .' % (root_dir+'/GO_ENRICHMENT/DOWN/'+temp+'/'+temp+'_down.DEG_Enriched_GO_classification_gene_count.txt'))
assert not os.system('cp /PUBLIC/source/RNA/noRef/readmefiles_v3/11.1.GOenrichment.README.txt %s' % (root_dir+'/'+project+'_result/'+project+'_results/'+str(result_order)+'.DEG_GOenrichment/'+str(result_order)+'.1.GOenrichment/'+str(result_order)+'.1.GOenrichment.README.txt'))
os.chdir('..')

assert not os.system('mkdir '+str(result_order)+'.2.DAG')
os.chdir(''+str(result_order)+'.2.DAG')
for each in compare_names:
	temp='vs'.join(each)
	os.system('cp %s .' % (root_dir+'/GO_ENRICHMENT/ALL/'+temp+'/'+temp+'.DEG_Enriched_GO_bp_DAG.png'))
	os.system('cp %s .' % (root_dir+'/GO_ENRICHMENT/ALL/'+temp+'/'+temp+'.DEG_Enriched_GO_bp_DAG.pdf'))
	os.system('cp %s .' % (root_dir+'/GO_ENRICHMENT/ALL/'+temp+'/'+temp+'.DEG_Enriched_GO_cc_DAG.png'))
	os.system('cp %s .' % (root_dir+'/GO_ENRICHMENT/ALL/'+temp+'/'+temp+'.DEG_Enriched_GO_cc_DAG.pdf'))
	os.system('cp %s .' % (root_dir+'/GO_ENRICHMENT/ALL/'+temp+'/'+temp+'.DEG_Enriched_GO_mf_DAG.png'))
	os.system('cp %s .' % (root_dir+'/GO_ENRICHMENT/ALL/'+temp+'/'+temp+'.DEG_Enriched_GO_mf_DAG.pdf'))
	os.system('cp %s .' % (root_dir+'/GO_ENRICHMENT/UP/'+temp+'/'+temp+'_up.DEG_Enriched_GO_bp_DAG.png'))
	os.system('cp %s .' % (root_dir+'/GO_ENRICHMENT/UP/'+temp+'/'+temp+'_up.DEG_Enriched_GO_bp_DAG.pdf'))
	os.system('cp %s .' % (root_dir+'/GO_ENRICHMENT/UP/'+temp+'/'+temp+'_up.DEG_Enriched_GO_cc_DAG.png'))
	os.system('cp %s .' % (root_dir+'/GO_ENRICHMENT/UP/'+temp+'/'+temp+'_up.DEG_Enriched_GO_cc_DAG.pdf'))
	os.system('cp %s .' % (root_dir+'/GO_ENRICHMENT/UP/'+temp+'/'+temp+'_up.DEG_Enriched_GO_mf_DAG.png'))
	os.system('cp %s .' % (root_dir+'/GO_ENRICHMENT/UP/'+temp+'/'+temp+'_up.DEG_Enriched_GO_mf_DAG.pdf'))
	os.system('cp %s .' % (root_dir+'/GO_ENRICHMENT/DOWN/'+temp+'/'+temp+'_down.DEG_Enriched_GO_bp_DAG.png'))
	os.system('cp %s .' % (root_dir+'/GO_ENRICHMENT/DOWN/'+temp+'/'+temp+'_down.DEG_Enriched_GO_bp_DAG.pdf'))
	os.system('cp %s .' % (root_dir+'/GO_ENRICHMENT/DOWN/'+temp+'/'+temp+'_down.DEG_Enriched_GO_cc_DAG.png'))
	os.system('cp %s .' % (root_dir+'/GO_ENRICHMENT/DOWN/'+temp+'/'+temp+'_down.DEG_Enriched_GO_cc_DAG.pdf'))
	os.system('cp %s .' % (root_dir+'/GO_ENRICHMENT/DOWN/'+temp+'/'+temp+'_down.DEG_Enriched_GO_mf_DAG.png'))
	os.system('cp %s .' % (root_dir+'/GO_ENRICHMENT/DOWN/'+temp+'/'+temp+'_down.DEG_Enriched_GO_mf_DAG.pdf'))
assert not os.system('cp /PUBLIC/source/RNA/noRef/readmefiles_v3/11.2.DAG.README.txt %s' % (root_dir+'/'+project+'_result/'+project+'_results/'+str(result_order)+'.DEG_GOenrichment/'+str(result_order)+'.2.DAG/'+str(result_order)+'.2.DAG.README.txt'))
os.chdir(root_dir+'/'+project+'_result/'+project+'_results')
'''

if set([1,2,5,6,8]).issubset(includes):
	code+='''
result_order+=1
assert not os.system('mkdir '+str(result_order)+'.DEG_KEGGenrichment')
os.chdir(''+str(result_order)+'.DEG_KEGGenrichment')
for each in compare_names:
	temp='vs'.join(each)
	os.system('cp %s .' % (root_dir+'/KEGG/ALL/'+temp+'/'+temp+'.DEG_KEGG_pathway_enrichment_result.xls'))
	os.system('cp %s .' % (root_dir+'/KEGG/ALL/'+temp+'/'+temp+'.DEG_enriched_KEGG_pathway_top20.xls'))
	os.system('cp %s .' % (root_dir+'/KEGG/ALL/'+temp+'/'+temp+'.DEG_enriched_KEGG_pathway_scatterplot.png'))
	os.system('cp %s .' % (root_dir+'/KEGG/ALL/'+temp+'/'+temp+'.DEG_enriched_KEGG_pathway_scatterplot.pdf'))
	os.system('mkdir '+temp+'_kegg_web')
	os.system('cp %s %s' % (root_dir+'/KEGG/ALL/'+temp+'/'+temp+'.DEG_KEGG_pathway_enrichment_add.xls_rendered_html_detail.html','./'+temp+'_kegg_web'+'/'+temp+'.DEG_enriched_KEGG_pathway_API.html'))
	os.system('cp -r %s %s' % (root_dir+'/KEGG/ALL/'+temp+'/src','./'+temp+'_kegg_web'))

	os.system('cp %s .' % (root_dir+'/KEGG/UP/'+temp+'/'+temp+'_up.DEG_KEGG_pathway_enrichment_result.xls'))
	os.system('cp %s .' % (root_dir+'/KEGG/UP/'+temp+'/'+temp+'_up.DEG_enriched_KEGG_pathway_top20.xls'))
	os.system('cp %s .' % (root_dir+'/KEGG/UP/'+temp+'/'+temp+'_up.DEG_enriched_KEGG_pathway_scatterplot.png'))
	os.system('cp %s .' % (root_dir+'/KEGG/UP/'+temp+'/'+temp+'_up.DEG_enriched_KEGG_pathway_scatterplot.pdf'))
	os.system('mkdir '+temp+'_up_kegg_web')
	os.system('cp %s %s' % (root_dir+'/KEGG/UP/'+temp+'/'+temp+'_up.DEG_KEGG_pathway_enrichment_add.xls_rendered_html_detail.html','./'+temp+'_up_kegg_web'+'/'+temp+'_up.DEG_enriched_KEGG_pathway_API.html'))
	os.system('cp -r %s %s' % (root_dir+'/KEGG/UP/'+temp+'/src','./'+temp+'_up_kegg_web'))
	
	os.system('cp %s .' % (root_dir+'/KEGG/DOWN/'+temp+'/'+temp+'_down.DEG_KEGG_pathway_enrichment_result.xls'))
	os.system('cp %s .' % (root_dir+'/KEGG/DOWN/'+temp+'/'+temp+'_down.DEG_enriched_KEGG_pathway_top20.xls'))
	os.system('cp %s .' % (root_dir+'/KEGG/DOWN/'+temp+'/'+temp+'_down.DEG_enriched_KEGG_pathway_scatterplot.png'))
	os.system('cp %s .' % (root_dir+'/KEGG/DOWN/'+temp+'/'+temp+'_down.DEG_enriched_KEGG_pathway_scatterplot.pdf'))
	os.system('mkdir '+temp+'_down_kegg_web')
	os.system('cp %s %s' % (root_dir+'/KEGG/DOWN/'+temp+'/'+temp+'_down.DEG_KEGG_pathway_enrichment_add.xls_rendered_html_detail.html','./'+temp+'_down_kegg_web'+'/'+temp+'_down.DEG_enriched_KEGG_pathway_API.html'))
	os.system('cp -r %s %s' % (root_dir+'/KEGG/DOWN/'+temp+'/src','./'+temp+'_down_kegg_web'))

assert not os.system('cp /PUBLIC/source/RNA/noRef/readmefiles_v3/12.DEG.KEGGenrichment.README.txt %s' % (root_dir+'/'+project+'_result/'+project+'_results/'+str(result_order)+'.DEG_KEGGenrichment/'+str(result_order)+'.DEG.KEGGenrichment.README.txt'))
os.chdir(root_dir+'/'+project+'_result/'+project+'_results')
'''
if set([1,2,5,6,10]).issubset(includes):
        code+='''
result_order+=1
assert not os.system('mkdir '+str(result_order)+'.DEG_PPI')
os.chdir(''+str(result_order)+'.DEG_PPI')
for each in compare_names:
        temp='vs'.join(each)
        os.system('cp %s .' % (root_dir+'/PPI/ALL/'+temp+'/'+temp+'.PPI.txt'))
        os.system('cp %s .' % (root_dir+'/PPI/UP/'+temp+'/'+temp+'_up.PPI.txt'))
        os.system('cp %s .' % (root_dir+'/PPI/DOWN/'+temp+'/'+temp+'_down.PPI.txt'))
assert not os.system('cp /PUBLIC/source/RNA/noRef/readmefiles_v2/13.DEG.PPI.README.txt %s' % (root_dir+'/'+project+'_result/'+project+'_results/'+str(result_order)+'.DEG_PPI/'+str(result_order)+'.DEG.PPI.README.txt'))
assert not os.system('cp /PUBLIC/source/RNA/noRef/readmefiles_v3/CytoscapeQuickStart.pdf %s' % (root_dir+'/'+project+'_result/'+project+'_results/'+str(result_order)+'.DEG_PPI/CytoscapeQuickStart.pdf'))
os.chdir(root_dir+'/'+project+'_result/'+project+'_results')
'''
code+='''
os.chdir(root_dir+'/'+project+'_result')
#assert not os.system('cp /PUBLIC/source/RNA/noRef/Pipeline_noRef/NovoQuery.exe .')
#assert not os.system('cp /PUBLIC/source/RNA/noRef/Pipeline_noRef/NovoQuery_manual.pdf .')
assert not os.system('cp /PUBLIC/source/RNA/noRef/Pipeline_noRef/Novofinder.exe .')
assert not os.system('cp /PUBLIC/source/RNA/noRef/Pipeline_noRef/Novofinder_manual.pdf .')
os.chdir(root_dir)

'''

code+='''

print project+' report generating...'
os.chdir(root_dir)
os.chdir(project+'_result/'+project+'_results')
os.system('tree -d -v -L 2 -o ../../result_tree.html -H ../../'+project+'_results')
os.chdir(root_dir)
# generate the noref html report   (separate the processes from above for clearness)
assert not os.system('mkdir '+project+'_result/'+project+'_report')
os.chdir(project+'_result/'+project+'_report')
assert not os.system('cp -r /PUBLIC/source/RNA/noRef/Pipeline_noRef/report_v2/src .')
render={}
if (ss =='0.5'):
	os.system('cp /PUBLIC/source/RNA/noRef/Pipeline_noRef/method/Methods-noRef.pdf src/Methods-noRef_muti.pdf')
	render['flag_strand']=False
else:
	os.system('cp /PUBLIC/source/RNA/noRef/Pipeline_noRef/method/Methods-noRef_directional.pdf src/Methods-noRef_muti.pdf')
	render['flag_strand']=True

assert not os.system('perl /PUBLIC/source/RNA/noRef/Pipeline_noRef/result_tree.pl ../../result_tree.html src/result_tree.html')
render['title']=project
render['flag_qc']=False
render['flag_trinity']=False
render['flag_annotation']=False
render['flag_cds']=False
render['flag_snp']=False
render['flag_ssr']=False
render['flag_rsem']=False
render['flag_rnaseqqc']=False
render['flag_diffexp']=False
render['flag_go']=False
render['flag_kegg']=False
render['flag_ppi']=False
render['flag_flow']=False

render['flag_pe']=True

render['flag_dge']=False

figure_order=0
table_order=0
html_order=0

render['flag_qc']=True
render['flag_gc']=flag_gc
html_order+=1  #for raw data stat
html_order+=1  #for qc
render['html_order_qc']=str(html_order)
render['html_order_qc']=str(html_order)
figure_order+=1
render['figure_order_qc1']=str(figure_order)
figure_order+=1
if flag_gc:
	render['figure_order_qc2']=str(figure_order)
	figure_order+=1
render['figure_order_qc3']=str(figure_order)
table_order+=1
render['table_order_qc']=str(table_order)

render['figure1']=[]
render['figure2']=[]
render['figure3']=[]
render['table1']=[]
os.chdir(root_dir+'/'+project+'_result/'+project+'_report/src')
for eachsample in samples:
	assert not os.system('convert -resize 600 %s ./images/%s.error_rate_distribution.png' % (root_dir+'/QC/'+eachsample+'/files/'+eachsample+'.Error.png',eachsample))
	assert not os.system('convert -resize 90 %s ./images/%s.error_rate_distribution.JPEG' % (root_dir+'/QC/'+eachsample+'/files/'+eachsample+'.Error.png',eachsample))
	render['figure1'].append(["'"+'images/'+eachsample+'.error_rate_distribution.png'+"'","'"+'images/'+eachsample+'.error_rate_distribution.JPEG'+"'"])
	if flag_gc:
		assert not os.system('convert -resize 600 %s ./images/%s.GC_content_distribution.png' % (root_dir+'/QC/'+eachsample+'/files/'+eachsample+'.GC.png',eachsample))
		assert not os.system('convert -resize 90 %s ./images/%s.GC_content_distribution.JPEG' % (root_dir+'/QC/'+eachsample+'/files/'+eachsample+'.GC.png',eachsample))
		render['figure2'].append(["'"+'images/'+eachsample+'.GC_content_distribution.png'+"'","'"+'images/'+eachsample+'.GC_content_distribution.JPEG'+"'"])
	assert not os.system('convert -resize 600 %s ./images/%s.raw_reads_classification.png' % (root_dir+'/QC/'+eachsample+'/files/'+eachsample+'.pie3d.png',eachsample))
	assert not os.system('convert -resize 90 %s ./images/%s.raw_reads_classification.JPEG' % (root_dir+'/QC/'+eachsample+'/files/'+eachsample+'.pie3d.png',eachsample))
	render['figure3'].append(["'"+'images/'+eachsample+'.raw_reads_classification.png'+"'","'"+'images/'+eachsample+'.raw_reads_classification.JPEG'+"'"])
	
	for eachLine in open(root_dir+'/QC/'+eachsample+'/files/dataTable'):
		if eachLine.strip() != '':
			render['table1'].append(eachLine.strip().split())
'''

if set([1]).issubset(includes):
	code+='''
render['flag_trinity']=True
render['flag_flow']=True
render['minkmercov']=minkmercov
html_order+=1
render['html_order_trinity']=str(html_order)

render['TRINITY']=''.join(fasta_trim(root_dir+'/TRINITY/trinity.out/Trinity.fasta',3,2))


table_order+=1
render['table_order_trinity']=str(table_order)			
render['table2']=open(root_dir+'/TRINITY/assembly_INFO/Transcript_length_stat.txt').readlines()[1].strip().split()
render['table2_1']=open(root_dir+'/TRINITY/assembly_INFO/Unigene_length_stat.txt').readlines()[1].strip().split()

table_order+=1
render['table_order_trinity1']=str(table_order)	
render['table2_2']=open(root_dir+'/TRINITY/assembly_INFO/Transcript_n50.txt').readlines()[1].strip().split()
render['table2_3']=open(root_dir+'/TRINITY/assembly_INFO/Unigene_n50.txt').readlines()[1].strip().split()

figure_order+=1
render['figure_order_trinity']=str(figure_order)
render['figure3_1']=[]
assert not os.system('convert -resize 600 %s ./images/Unigene_Transcript.length_distribution.png' % (root_dir+'/TRINITY/assembly_INFO/Unigene_Transcript.length_distribution.png'))
assert not os.system('convert -resize 90 %s ./images/Unigene_Transcript.length_distribution.JPEG' % (root_dir+'/TRINITY/assembly_INFO/Unigene_Transcript.length_distribution.png'))
assert not os.system('convert -resize 600 %s ./images/Transcript_bar.png' % (root_dir+'/TRINITY/assembly_INFO/Transcript_bar.png'))
assert not os.system('convert -resize 90 %s ./images/Transcript_bar.JPEG' % (root_dir+'/TRINITY/assembly_INFO/Transcript_bar.png'))
assert not os.system('convert -resize 600 %s ./images/Unigene_bar.png' % (root_dir+'/TRINITY/assembly_INFO/Unigene_bar.png'))
assert not os.system('convert -resize 90 %s ./images/Unigene_bar.JPEG' % (root_dir+'/TRINITY/assembly_INFO/Unigene_bar.png'))
render['figure3_1'].append(["'"+'images/Unigene_Transcript.length_distribution.png'+"'","'"+'images/Unigene_Transcript.length_distribution.JPEG'+"'"])
render['figure3_1'].append(["'"+'images/Transcript_bar.png'+"'","'"+'images/Transcript_bar.JPEG'+"'"])
render['figure3_1'].append(["'"+'images/Unigene_bar.png'+"'","'"+'images/Unigene_bar.JPEG'+"'"])

''' 

if set([1,2]).issubset(includes):
	code+='''
render['flag_annotation']=True
html_order+=1
render['html_order_annotation']=str(html_order)

table_order+=1
render['table_order_annotation0']=str(table_order)
render['table3_0']=[]
cont=open(root_dir+'/ANNOTATION_ALL/Annotation.statistics').readlines()
cont.pop(0)
for eachLine in cont:
	if eachLine.strip() != '':
		render['table3_0'].append(eachLine.strip().split('\t'))

		
table_order+=1
render['table_order_annotation']=str(table_order)
render['table3']=[]
k=5
cont=open(root_dir+'/ANNOTATION_ALL/Annotation_Summary.xls').readlines()
cont.pop(0)
for eachLine in cont:
	if eachLine.strip() != '':
		temp=eachLine.split()
		if all([each != '--' for each in temp[:6]]):
			render['table3'].append(temp[:6])
			k-=1
	if k == 0:
		break

figure_order+=1
render['figure_order_annotaion0']=str(figure_order)
render['figure4_1']=[]
assert not os.system('convert -resize 600 %s ./images/Species_classification.3dpie.png' % (root_dir+'/ANNOTATION_ALL/Species_classification.3dpie.png'))
assert not os.system('convert -resize 90 %s ./images/Species_classification.3dpie.JPEG' % (root_dir+'/ANNOTATION_ALL/Species_classification.3dpie.png'))
assert not os.system('convert -resize 600 %s ./images/Evalue_Distribution.3dpie.png' % (root_dir+'/ANNOTATION_ALL/Evalue_Distribution.3dpie.png'))
assert not os.system('convert -resize 90 %s ./images/Evalue_Distribution.3dpie.JPEG' % (root_dir+'/ANNOTATION_ALL/Evalue_Distribution.3dpie.png'))
assert not os.system('convert -resize 600 %s ./images/Similarity_Distribution.3dpie.png' % (root_dir+'/ANNOTATION_ALL/Similarity_Distribution.3dpie.png'))
assert not os.system('convert -resize 90 %s ./images/Similarity_Distribution.3dpie.JPEG' % (root_dir+'/ANNOTATION_ALL/Similarity_Distribution.3dpie.png'))
render['figure4_1'].append(["'"+'images/Species_classification.3dpie.png'+"'","'"+'images/Species_classification.3dpie.JPEG'+"'"])
render['figure4_1'].append(["'"+'images/Evalue_Distribution.3dpie.png'+"'","'"+'images/Evalue_Distribution.3dpie.JPEG'+"'"])
render['figure4_1'].append(["'"+'images/Similarity_Distribution.3dpie.png'+"'","'"+'images/Similarity_Distribution.3dpie.JPEG'+"'"])


figure_order+=1
render['figure_order_annotaion1']=str(figure_order)
assert not os.system('cp %s ./images/GO_classification_bar.png' % (root_dir+'/ANNOTATION_ALL/GO_classification_bar.png'))
# render['figure5']='images/GO_classification_bar.png'

figure_order+=1
render['figure_order_annotaion2']=str(figure_order)
assert not os.system('cp %s ./images/KOG_classification.png' % (root_dir+'/ANNOTATION_ALL/KOG_classification.png'))
# render['figure6']='images/KOG_classification.png'

figure_order+=1
render['figure_order_annotaion3']=str(figure_order)
assert not os.system('cp %s ./images/KEGG_classification.png' % (root_dir+'/ANNOTATION_ALL/KEGG_classification.png'))
# render['figure7']

render['CDS1']=''.join(fasta_trim(root_dir+'/ANNOTATION_ALL/CDSprediction/'+project+'.blast.cds.fasta',3,2))
render['CDS2']=''.join(fasta_trim(root_dir+'/ANNOTATION_ALL/CDSprediction/'+project+'.estscan.cds.fasta',3,2))

render['flag_cds']=True
html_order+=1 #for CDS prediction
render['html_order_cds']=str(html_order)
'''

if set([1,3]).issubset(includes):
	code+='''
render['flag_snp']=True
html_order+=1
render['html_order_snp']=str(html_order)
table_order+=1
render['table_order_snp']=str(table_order)
render['table4']=[]
k=5
render['SNP_head']='</th><th>'.join(samples)
cont=open(root_dir+'/SNP/ResultsQ30/SNPs.xls').readlines()
cont.pop(0)
for eachLine in cont:
	if (not eachLine.startswith('#')) and eachLine.strip() != '':
		temp=eachLine.split()
		snp='</td><td>'.join(temp)
		render['table4'].append(snp)
		k-=1
	if k == 0:
		break

'''

if set([1,4]).issubset(includes):
	code+='''
render['flag_ssr']=True
html_order+=1
render['html_order_ssr']=str(html_order)
table_order+=1
render['table_order_ssr1']=str(table_order)
table_order+=1
render['table_order_ssr2']=str(table_order)
render['table5']=[]
render['table6']=[]
k=5
for eachLine in open(root_dir+'/SSR/SSR_Primer/'+project+'.SSR.pos.misa'):
        if (not eachLine.startswith('ID')) and eachLine.strip() != '':
                temp=[each.strip() for each in eachLine.split()]
                render['table5'].append(temp)
                k-=1    
        if k == 0:
                break
k=5
for eachLine in open(root_dir+'/SSR/SSR_Primer/unigene.fasta.SSR_primer'):
	if (not eachLine.startswith('ID')) and eachLine.strip() != '':
		temp=[each.strip() for each in eachLine.split()]
		if len(temp[0:16]) == 16:
			render['table6'].append([temp[0]]+temp[7:16])
			k-=1	
	if k == 0:
		break

figure_order+=1
render['figure_order_ssr1']=str(figure_order)
assert not os.system('cp %s ./images/unigene.fasta.SSR_motifs_distribution.png' % (root_dir+'/SSR/SSR_Primer/unigene.fasta.SSR_motifs_distribution.png'))
# render['figure9']
'''

if set([1,5]).issubset(includes):
	code+='''
render['flag_rsem']=True
html_order+=1
render['html_order_rsem']=str(html_order)
table_order+=1
render['table_order_mapping']=str(table_order)
os.chdir(root_dir)
assert not os.system('perl /PUBLIC/source/RNA/noRef/RSEM/mapping_rate_trans_v2.pl -n %s' % sample)
assert not os.system('cp %s %s' % (root_dir+'/RSEM/mapping_rate.txt',root_dir+'/'+project+'_result/'+project+'_results/'+str(RSEM_order)+'.GeneExprQuantification/Mapping_rate.xls'))
render['table_mapping']=[]
f_map=open('./RSEM/mapping_rate.txt').readlines()
f_map.pop(0)
for eachLine in f_map:
        render['table_mapping'].append(eachLine.strip().split('\t'))

os.chdir(project+'_result/'+project+'_report')


table_order+=1
render['table_order_rsem']=str(table_order)
render['table7']=[]
sample_name=samples[0]
k=5
for eachLine in open(root_dir+'/RSEM/'+sample_name+'/'+sample_name+'.RSEM.out/'+sample_name+'.Readcount_FPKM.xls'):
	if (not eachLine.startswith('gene_id')) and eachLine.strip() != '':
		render['table7'].append(eachLine.strip().split())
		k-=1
	if k == 0:
		break

figure_order+=1
render['figure_order_rsem']=str(figure_order)
render['figure10']=[]
os.chdir(root_dir+'/'+project+'_result/'+project+'_report/src')
for eachsample in samples:
	assert not os.system('convert -resize 600 %s ./images/%s.FPKM_density_distribution.png' % (root_dir+'/RSEM/'+eachsample+'/'+eachsample+'.RSEM.out/'+eachsample+'.FPKM_density_distribution.png',eachsample))
	assert not os.system('convert -resize 90 %s %s' % (root_dir+'/RSEM/'+eachsample+'/'+eachsample+'.RSEM.out/'+eachsample+'.FPKM_density_distribution.png','./images/'+eachsample+'.FPKM_density_distribution.JPEG'))
	render['figure10'].append(["'"+'images/'+eachsample+'.FPKM_density_distribution.png'+"'","'"+'images/'+eachsample+'.FPKM_density_distribution.JPEG'+"'"])
'''

if set([1,2,5,7]).issubset(includes):
	code+='''
render['flag_rnaseqqc']=True
render['flag_uniform']=flag_uniform
render['flag_rnaseqqc_cor']=False
html_order+=1
render['html_order_rnaseqqc']=str(html_order)
figure_order+=1
render['figure_order_rnaseqqc1']=str(figure_order)
if flag_uniform:
	figure_order+=1
	render['figure_order_rnaseqqc2']=str(figure_order)
render['figure11']=[]
render['figure12']=[]
for eachsample in samples:
	assert not os.system('convert -resize 600 %s ./images/%s.Saturation_curve.png' % (root_dir+'/RNA_SEQ_QC/RNA_seq_QC/'+eachsample+'.Saturation_curve.png',eachsample))
	assert not os.system('convert -resize 90 %s %s' % (root_dir+'/RNA_SEQ_QC/RNA_seq_QC/'+eachsample+'.Saturation_curve.png','./images/'+eachsample+'.Saturation_curve.JPEG'))
	render['figure11'].append(["'"+'images/'+eachsample+'.Saturation_curve.png'+"'","'"+'images/'+eachsample+'.Saturation_curve.JPEG'+"'"])
	if flag_uniform:
		assert not os.system('convert -resize 600 %s ./images/%s.Uniform_distribution.png' % (root_dir+'/RNA_SEQ_QC/RNA_seq_QC/'+eachsample+'.Uniform_distribution.png',eachsample))
		assert not os.system('convert -resize 90 %s %s' % (root_dir+'/RNA_SEQ_QC/RNA_seq_QC/'+eachsample+'.Uniform_distribution.png','./images/'+eachsample+'.Uniform_distribution.JPEG'))
		render['figure12'].append(["'"+'images/'+eachsample+'.Uniform_distribution.png'+"'","'"+'images/'+eachsample+'.Uniform_distribution.JPEG'+"'"])
'''

if set([1,2,5,6,7]).issubset(includes):
	code+='''
render['flag_rnaseqqc_cor']=True
figure_order+=1
render['figure_order_rnaseqqc_cor']=str(figure_order)
render['figure13']=[]
assert not os.system('convert -resize 600 %s %s' % (root_dir+'/DIFF_EXP/Diff_analysis.out/corr_plot/cor_pearson.png','./images/cor_pearson.png'))
assert not os.system('convert -resize 90 %s %s' % (root_dir+'/DIFF_EXP/Diff_analysis.out/corr_plot/cor_pearson.png','./images/cor_pearson.JPEG'))
render['figure13'].append(["'"+'images/cor_pearson.png'+"'","'"+'images/cor_pearson.JPEG'+"'"])
for png in glob.iglob(root_dir+'/DIFF_EXP/Diff_analysis.out/corr_plot/*.png'):
	if png !=root_dir+'/DIFF_EXP/Diff_analysis.out/corr_plot/cor_pearson.png':
		name=re.search(r'(/.*/)(.*)\.png',png).group(2)
		assert not os.system('convert -resize 600 %s %s' % (png,'./images/'+name+'.png'))
		assert not os.system('convert -resize 90 %s %s' % (png,'./images/'+name+'.JPEG'))
		render['figure13'].append(["'"+'images/'+os.path.basename(png)+"'","'"+'images/'+name+'.JPEG'+"'"])
'''

if set([1,5,6]).issubset(includes):
	code+='''

render['flag_diffexp']=True
html_order+=1
render['html_order_diffexp']=str(html_order)
figure_order+=1
render['figure_order_diffexp1']=str(figure_order)
render['figure14_1']=[]
assert not os.system('convert -resize 600 %s ./images/FPKM_density_distribution.png' % (root_dir+'/DIFF_EXP/Diff_analysis.out/FPKM_density_distribution.png'))
assert not os.system('convert -resize 90 %s ./images/FPKM_density_distribution.JPEG' % (root_dir+'/DIFF_EXP/Diff_analysis.out/FPKM_density_distribution.png'))
assert not os.system('convert -resize 600 %s ./images/FPKM_boxplot.png' % (root_dir+'/DIFF_EXP/Diff_analysis.out/FPKM_boxplot.png'))
assert not os.system('convert -resize 90 %s ./images/FPKM_boxplot.JPEG' % (root_dir+'/DIFF_EXP/Diff_analysis.out/FPKM_boxplot.png'))
render['figure14_1'].append(["'"+'images/FPKM_density_distribution.png'+"'","'"+'images/FPKM_density_distribution.JPEG'+"'"])
render['figure14_1'].append(["'"+'images/FPKM_boxplot.png'+"'","'"+'images/FPKM_boxplot.JPEG'+"'"])
table_order+=1
render['table_order_diffexp']=str(table_order)
render['flag_repeat']=flag_repeat
render['flag_venn']=flag_venn

render['table8']=[]
com=compare_names[0]
render['com1']=com[0]
render['com2']=com[1]
com_temp='vs'.join(com)
k=5
for eachLine in open(root_dir+'/DIFF_EXP/Diff_analysis.out/'+com_temp+'/'+com_temp+'.DEG.xls'):
	if (not eachLine.upper().startswith('GENE_ID')) and eachLine.strip() != '':
		temp=eachLine.split()
		if flag_repeat:
			render['table8'].append([temp[0],temp[1],temp[2],temp[3],temp[4],temp[5]])
		else:
			render['table8'].append([temp[0],temp[1],temp[2],temp[3],temp[4],temp[5]])
		k-=1
	if k == 0:
		break

figure_order+=1
render['figure_order_diffexp2']=str(figure_order)
render['figure15']=[]
for each in compare_names:
	temp='vs'.join(each)
	assert not os.system('convert -resize 600 %s ./images/%s.Volcanoplot.png' % (root_dir+'/DIFF_EXP/Diff_analysis.out/'+temp+'/'+temp+'.Volcanoplot.png',temp))
	assert not os.system('convert -resize 90 %s %s' % (root_dir+'/DIFF_EXP/Diff_analysis.out/'+temp+'/'+temp+'.Volcanoplot.png','./images/'+temp+'.Volcanoplot.JPEG'))
	render['figure15'].append(["'"+'images/'+temp+'.Volcanoplot.png'+"'","'"+'images/'+temp+'.Volcanoplot.JPEG'+"'"])

if flag_venn:
	figure_order+=1
	render['figure_order_diffexp3']=str(figure_order)
	render['figure16']=[]
	for each in venn_cluster_vs_names:
		assert not os.system('convert -resize 600 %s ./images/%s.DEG_Venn_diagram.png' % (root_dir+'/DIFF_EXP/Diff_analysis.out/'+each+'/'+each+'.DEG_Venn_diagram.png',each))
		assert not os.system('convert -resize 90 %s %s' % (root_dir+'/DIFF_EXP/Diff_analysis.out/'+each+'/'+each+'.DEG_Venn_diagram.png','./images/'+each+'.DEG_Venn_diagram.JPEG'))
		render['figure16'].append(["'"+'images/'+each+'.DEG_Venn_diagram.png'+"'","'"+'images/'+each+'.DEG_Venn_diagram.JPEG'+"'"])

figure_order+=1
render['figure_order_diffexp4']=str(figure_order)	
render['figure16_1']=[]
assert not os.system('convert -resize 600 %s ./images/Hcluster_heatmap.png' % (root_dir+'/DIFF_EXP/Diff_analysis.out/Hcluster_heatmap.png'))
assert not os.system('convert -resize 90 %s ./images/Hcluster_heatmap.JPEG' % (root_dir+'/DIFF_EXP/Diff_analysis.out/Hcluster_heatmap.png'))
assert not os.system('convert -resize 600 %s ./images/H_cluster.png' % (root_dir+'/DIFF_EXP/Diff_analysis.out/H_cluster/H_cluster.png'))
assert not os.system('convert -resize 90 %s ./images/H_cluster.JPEG' % (root_dir+'/DIFF_EXP/Diff_analysis.out/H_cluster/H_cluster.png'))
render['figure16_1'].append(["'"+'images/Hcluster_heatmap.png'+"'","'"+'images/Hcluster_heatmap.JPEG'+"'"])
render['figure16_1'].append(["'"+'images/H_cluster.png'+"'","'"+'images/H_cluster.JPEG'+"'"])
'''	

if set([1,2,5,6,9]).issubset(includes):
	code+='''
render['flag_go']=True
html_order+=1
render['html_order_go']=str(html_order)
table_order+=1
render['table_order_go']=str(table_order)
render['table9']=[]
k=4
for eachLine in open(root_dir+'/GO_ENRICHMENT/ALL/'+com_temp+'/'+com_temp+'.DEG_GO_enrichment_result.xls'):
	if (not eachLine.startswith('GO_accession')) and eachLine.strip() != '':
		temp=eachLine.split('\t')
		render['table9'].append(temp[:7])
		k-=1
	if k == 0:
		break

figure_order+=1
render['figure_order_go1']=str(figure_order)
render['figure18']=[]
render['figure18_lable']=[]
figure_order+=1
render['figure_order_go2']=str(figure_order)
render['figure19']=[]
render['figure19_lable']=[]

for each in compare_names:
	temp='vs'.join(each)
	os.system('convert -resize 600 %s ./images/%s.DEG_Enriched_GO_classification.png' % (root_dir+'/GO_ENRICHMENT/ALL/'+temp+'/'+temp+'.DEG_Enriched_GO_classification.png',temp))
	os.system('convert -resize 90 %s %s' % (root_dir+'/GO_ENRICHMENT/ALL/'+temp+'/'+temp+'.DEG_Enriched_GO_classification.png','./images/'+temp+'.DEG_Enriched_GO_classification.JPEG'))
	render['figure18'].append(["'"+'images/'+temp+'.DEG_Enriched_GO_classification.png'+"'","'"+'images/'+temp+'.DEG_Enriched_GO_classification.JPEG'+"'"])
	render['figure18_lable'].append(temp)
	os.system('convert -resize 600 %s ./images/%s.DEG_Enriched_GO_bp_DAG.png' % (root_dir+'/GO_ENRICHMENT/ALL/'+temp+'/'+temp+'.DEG_Enriched_GO_bp_DAG.png',temp))
	os.system('convert -resize 90 %s %s' % (root_dir+'/GO_ENRICHMENT/ALL/'+temp+'/'+temp+'.DEG_Enriched_GO_bp_DAG.png','./images/'+temp+'.DEG_Enriched_GO_bp_DAG.JPEG'))
	render['figure19'].append(["'"+'images/'+temp+'.DEG_Enriched_GO_bp_DAG.png'+"'","'"+'images/'+temp+'.DEG_Enriched_GO_bp_DAG.JPEG'+"'"])
	os.system('convert -resize 600 %s ./images/%s.DEG_Enriched_GO_cc_DAG.png' % (root_dir+'/GO_ENRICHMENT/ALL/'+temp+'/'+temp+'.DEG_Enriched_GO_cc_DAG.png',temp))
	os.system('convert -resize 90 %s %s' % (root_dir+'/GO_ENRICHMENT/ALL/'+temp+'/'+temp+'.DEG_Enriched_GO_cc_DAG.png','./images/'+temp+'.DEG_Enriched_GO_cc_DAG.JPEG'))
	render['figure19'].append(["'"+'images/'+temp+'.DEG_Enriched_GO_cc_DAG.png'+"'","'"+'images/'+temp+'.DEG_Enriched_GO_cc_DAG.JPEG'+"'"])
	os.system('convert -resize 600 %s ./images/%s.DEG_Enriched_GO_mf_DAG.png' % (root_dir+'/GO_ENRICHMENT/ALL/'+temp+'/'+temp+'.DEG_Enriched_GO_mf_DAG.png',temp))
	os.system('convert -resize 90 %s %s' % (root_dir+'/GO_ENRICHMENT/ALL/'+temp+'/'+temp+'.DEG_Enriched_GO_mf_DAG.png','./images/'+temp+'.DEG_Enriched_GO_mf_DAG.JPEG'))
	render['figure19'].append(["'"+'images/'+temp+'.DEG_Enriched_GO_mf_DAG.png'+"'","'"+'images/'+temp+'.DEG_Enriched_GO_mf_DAG.JPEG'+"'"])
render['figure18_lable']=','.join(render['figure18_lable'])
render['figure19_lable']=render['figure18_lable']
'''

if set([1,2,5,6,8]).issubset(includes):
	code+='''
render['flag_kegg']=True
html_order+=1
render['html_order_kegg']=str(html_order)

table_order+=1
render['table_order_kegg']=str(table_order)	
flag=0
k=5
sig_map=''
kegg_png=[]
render['table10']=[]
for eachLine in open(root_dir+'/KEGG/ALL/'+com_temp+'/'+com_temp+'.DEG_KEGG_pathway_enrichment_result.xls'):
	if eachLine.startswith('#Term'):
		flag=1
		continue
	if flag == 1:
		if eachLine.strip() != '':
			temp=eachLine.split('\t')
			if sig_map == '':
				sig_map=temp[2].strip()
				sig_map_png=root_dir+'/KEGG/ALL/'+com_temp+'/src/'+sig_map+'.png'
				if os.path.isfile(sig_map_png):
					sig_map=temp[2].strip()
					sig_map_png=root_dir+'/KEGG/ALL/'+com_temp+'/src/'+sig_map+'.png'
					kegg_png.append(sig_map_png)
					sig_map = ''
				else:
					sig_map=''
			render['table10'].append(temp[:7])
			k-=1
		if k == 0:
			break

figure_order+=1
render['figure_order_kegg1']=str(figure_order)
render['figure20']=[]
render['figure20_lable']=[]
render['figure21']=[]
for each in compare_names:
	temp='vs'.join(each)
	os.system('convert -resize 600 %s ./images/%s.DEG_enriched_KEGG_pathway_scatterplot.png' % (root_dir+'/KEGG/ALL/'+temp+'/'+temp+'.DEG_enriched_KEGG_pathway_scatterplot.png',temp))
	os.system('convert -resize 90 %s %s' % (root_dir+'/KEGG/ALL/'+temp+'/'+temp+'.DEG_enriched_KEGG_pathway_scatterplot.png','./images/'+temp+'.DEG_enriched_KEGG_pathway_scatterplot.JPEG'))
	render['figure20'].append(["'"+'images/'+temp+'.DEG_enriched_KEGG_pathway_scatterplot.png'+"'","'"+'images/'+temp+'.DEG_enriched_KEGG_pathway_scatterplot.JPEG'+"'"])
	render['figure20_lable'].append(temp)
render['figure20_lable']=','.join(render['figure20_lable'])

figure_order+=1
render['figure_order_kegg2']=str(figure_order)		
for png in kegg_png:
	name=re.search(r'(/.*/)(.*)\.png',png).group(2)
	assert not os.system('convert -resize 600 %s %s' % (png,'./images/'+name+'.png'))
	assert not os.system('convert -resize 90 %s %s' % (png,'./images/'+name+'.JPEG'))
	render['figure21'].append(["'"+'images/'+os.path.basename(png)+"'","'"+'images/'+name+'.JPEG'+"'"])
'''
if set([1,2,5,6,10]).issubset(includes):
        code+='''
render['flag_ppi']=True
html_order+=1
render['html_order_ppi']=str(html_order)

figure_order+=1
render['figure_order_ppi']=str(figure_order)
'''

#code+='''
'''
id_order=0

id_order+=1    #for sequence
render['id_order_sequence']=str(id_order)
if render['flag_flow']:
        id_order+=1    #for pipeline
        render['id_order_pipeline']=str(id_order)
if render['flag_qc']:
        id_order+=1    #for raw
        render['id_order_raw']=str(id_order)
        id_order+=1    #for qc 
        render['id_order_qc']=str(id_order)
if render['flag_trinity']:
        id_order+=1    #for trinity
        render['id_order_trinity']=str(id_order)
if render['flag_annotation']:
        id_order+=1    #for annotation
        render['id_order_annotation']=str(id_order)
if render['flag_cds']:
        id_order+=1    #for cds
        render['id_order_cds']=str(id_order)
if render['flag_snp']:
        id_order+=1    #for snp
        render['id_order_snp']=str(id_order)
if render['flag_ssr']:
        id_order+=1    #for ssr
        render['id_order_ssr']=str(id_order)
if render['flag_rsem']:
        id_order+=1    #for rsem
        render['id_order_rsem']=str(id_order)
if render['flag_rnaseqqc']:
        id_order+=1    #for rnaseqqc
        render['id_order_rnaseqqc']=str(id_order)
if render['flag_diffexp']:
        id_order+=1    #for diffexp
        render['id_order_diffexp']=str(id_order)
if render['flag_go']:
        id_order+=1    #for go
        render['id_order_go']=str(id_order)
if render['flag_kegg']:
        id_order+=1    #for kegg
        render['id_order_kegg']=str(id_order)
if render['flag_ppi']:
        id_order+=1    #for ppi
        render['id_order_ppi']=str(id_order)
'''
code+='''
os.chdir(root_dir+'/'+project+'_result/'+project+'_report')
t=loader.get_template('src/right.html')
index=loader.get_template('index.html')
top=loader.get_template('src/top.html')
c=Context(render)
righthtml=t.render(c)
indexhtml=index.render(c)
tophtml=top.render(c)
open('src/right.html','w').write(righthtml)
open('src/top.html','w').write(tophtml)
open(project+'_Report.html','w').write(indexhtml)

assert not os.system("/PUBLIC/software/RNA/wkhtmltopdf/wkhtmltopdf --page-width 350mm --page-height 495mm -n --print-media-type --footer-center '[page] / [topage]' " +'src/right.html '+project+'_Report.pdf')
assert not os.system("sed -i 's/#00//g' "+project+'_Report.pdf')

############################################
os.chdir(root_dir+'/QC')
rdmd5=glob.glob('rd_md5.sh.o*')
if not rdmd5:
	assert not os.system('qsub -V -cwd -l vf=0.5G -l p=1 rd_md5.sh')
cdmd5=glob.glob('cd_md5.sh.o*')
if not cdmd5:
	assert not os.system('qsub -V -cwd -l vf=1G -l p=1 cd_md5.sh')

############################################
#noref data sum
os.chdir(root_dir)
assert not os.system('perl /PUBLIC/source/RNA/noRef/data_sum/noRef_sum_report.pl -n %s -p %s' % (sample,project))
'''
open('NOREF_step8.py','w').write(code)
##########################################################################
#for byebye (sh->py only, not run)  root_dir
code='python /PUBLIC/source/RNA/noRef/Pipeline_noRef/byebye.py --project '+project+'\n'

open('byebye.sh','w').write(code)
#for data_given (sh->py only, not run)	root_dir
code='perl /PUBLIC/source/RNA/noRef/Pipeline_noRef/release_data.pl -adir '+root_dir+' -type TR -s '+sample+'\n'
open ('data_give.sh','w').write(code)
