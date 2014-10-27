import sys
import os
import os.path
import re
import argparse
root_dir=os.getcwd()

##########################################################################
##########################################################################

#parse the arguments
parser = argparse.ArgumentParser(description="unigene DGE pipline v1.0")
parser.add_argument('--project',help='project name, maybe same with the name of root dir, which will be displayed in the final report title, [REQUIRED]',required=True)
parser.add_argument('--sample',help="sample names(sample1,sample2,...)  warning: order sensitive!",required=True)
parser.add_argument('--fq',help="the original directory of the raw fastq reads, [REQUIRED]",default=None)
parser.add_argument('--mapfile',help="mapfile files (mapfile1,mapfile2.The first column of mapfile is the NH number ,the second column of mapfile is the sample name)",default=None)
parser.add_argument('--raw_dir',help="the original directory of the raw fastq reads keep order in line with mapfile files (raw_dir1,raw_dir2)",default=None)
parser.add_argument('--ad',help="the original directory of the adapter list, ",default=None)
parser.add_argument('--generate_adapter',help='whether to generate the adpepter list y or n',choices=['y','n'],default='n')
parser.add_argument('--index',help='the index number of the fastq files,if you want to generate the adapter list,the parameter will be necessary,diffrent samples split by ,warning: order sensitive!',default=None)
parser.add_argument('--ss',help="strand specific, [REQUIRED]",choices=['0.0','0.5','1.0'],required=True)
parser.add_argument('--fa',help="the unigene FASTA file, [REQUIRED for mapping]",required=True)
parser.add_argument('--group',help="sample classification, e.g. sample1/sample2,sample3, [REQUIRED]",required=True)
parser.add_argument('--groupname',help="group names, [group1,group2,... ]",required=True)
parser.add_argument('--compare',help="group comparison strategy 1:2,1:3,..., [REQUIRED]",required=True)
parser.add_argument('--venn',help="venn and cluster mode, defult is all groups, 1:2_1:3,1:3_2:3,... ",default=None)
parser.add_argument('--goann',help="gene to GO annotations file, [REQUIRED]",required=True)
parser.add_argument('--koann',help="gene to KEGG annotations file, [REQUIRED]",default=None)
parser.add_argument('--species',help="abbreviation for species, for kegg pathway. (note: animal is 000, plant is 001, all species: /PUBLIC/database/RNA//kobas2.0-data-20120208/seq_pep)[REQUIRED]",default=None)
parser.add_argument('--ppi_number',help="species code, (ref file: /PUBLIC/database/RNA/string_ppi/species.v9.0.txt)",default=None)
parser.add_argument('--ppi_blast',help="whether to run blast to get the protein-protein interactions",choices=['y','n'],default=None)
parser.add_argument('--genenamefile',help="genenamefile, 1st col is geneID, 2nd col is genename",default=None)

# extract, check and display the parameters
argv = vars(parser.parse_args())
project=argv['project'].strip()
display = open(root_dir + '/' + project + '_command.txt','w')
display.write('project: %s\n' % (project))
sample=argv['sample'].strip()
samples=sample.split(',')
samples_tmp=list(set(samples))
assert len(samples)==len(samples_tmp)
display.write('sample: %s\n' % (sample))
#-----------------------------------------------------------------------
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
	fq_tmp=fq+'/'+each+'.fq.gz'
	assert os.path.isfile(fq_tmp)
display.write('fq: %s\n' % (fq))
display.write('adapter: ')
if argv['ad']:
	ad = argv['ad'].strip()
	ad = os.path.abspath(ad)
	for each in samples:
		ad_tmp=ad+'/'+each+'.adapter.list.gz'
		assert os.path.isfile(ad_tmp)
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
	if argv['index']==None:
		print 'Error:  the parameters --index and --generate_adapter are not consensus!\n'
		exit()
else:
	if argv['index'] != None:
		print 'Error:  the parameters --index and --generate_adapter are not consensus!\n'
		exit()
#-----------------------------------------------------------------------

ss = argv['ss'].strip()
display.write('ss: %s\n' % (ss))
fa = argv['fa'].strip()
fa = os.path.abspath(fa)
display.write('fa: %s\n' % (fa))
#-----------------------------------------------------------------------
group = argv['group'].strip()
groupname = argv['groupname'].strip()
compare = argv['compare'].strip()
if argv['venn']:
        venn = argv['venn'].strip()
else:
	venn = compare.replace(',','_')
groups=group.split(',')
groupnames=groupname.split(',')
compares=compare.split(',')
venns=venn.split(',')
for each in groups:
	temp=each.split('/')
	assert len(temp)==len(list(set(temp)))
	for each2 in temp:
		assert each2 in samples
assert len(groups)==len(groupnames)
display.write('\ngroups:\n')
assert len(groupnames)==len(list(set(groupnames)))
for i,each in enumerate(groupnames):
	group_tmp=groups[i]
	group_tmp=group_tmp.replace('/',',')
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

goann = argv['goann'].strip()
goann = os.path.abspath(goann)
assert os.path.isfile(goann)
display.write('goann: %s\n' % (goann))
if argv['koann']:
	if argv['species']:
		print 'Error:  the parameters --koann and --species are not consensus!\n'
		exit()
	else:
		koann=argv['koann'].strip()
		koann=os.path.abspath(koann)
		assert os.path.isfile(koann)
		display.write('koann: %s\n' % (koann))
else:
	if argv['species']==None:
		print 'Error:  the parameters --koann and --species are not consensus!\n'
		exit()
	else:
		species = argv['species'].strip()
		display.write('KEGG species: %s\n' % (species))
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
sample='%s'
fa='%s'
ss='%s'
root_dir='%s'
''' % (fq,sample,fa,ss,root_dir)

if argv['ad']:
	code+='''
ad='%s'
''' % (ad)

if generate_adapter != 'n':
	code+='''
index='%s'
	''' % (index)

if argv['mapfile']:
	code+'''
mapfile='%s'
''' % (root_dir+'/libraryID')	
code+='''
assert not os.system('mkdir QC_UG')
os.chdir('QC_UG')
f=open(root_dir+'/QC_UG/generate_QC.sh','w')
f.write('perl /PUBLIC/source/RNA/QC/QC_v2/allrunQC_v2.pl -fq %s -se-pe se -n %s -o %s -unigenefa %s -ss %s ' %(fq,sample,root_dir+'/QC_UG',fa,ss))
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
assert not os.system('sh %s/generate_QC.sh' %(root_dir+'/QC_UG'))

samples=sample.split(',')
for eachsample in samples:
	os.chdir(root_dir+'/QC_UG/'+eachsample)
	assert not os.system('qsub -V -l vf=5g,p=6 -cwd %s' % (eachsample+'_QC_noref_se.sh'))
os.chdir(root_dir)
'''
open('UG_step1_QC.py','w').write(code)
#################################################################################################
#for QCreport
code='''
import os

sample='%s'
project='%s'
root_dir='%s'
''' % (sample,project,root_dir)

code+='''

assert not os.system('mkdir QC_UG/QCreport')
os.system('sh /PUBLIC/source/RNA/QC/QC_v2/QCreport/UG_QCreport.sh -dir %s -sample %s -title %s -results %s' %(root_dir+'/QC_UG',sample,project,root_dir+'/QC_UG/QCreport'))
'''
open('UG_step1_QCreport.py','w').write(code)
#################################################################################################
#for Diff analysis and Curves
code='''
import re
import os

sample='%s'
group='%s'
groupname='%s'
compare='%s'
venn='%s'
fa='%s'
qc_dir='%s'
root_dir='%s'
''' % (sample,group,groupname,compare,venn,fa,root_dir+'/QC_UG/',root_dir)

if argv['genenamefile']:
	code+='''
genenamefile='%s'
''' % (genenamefile)

code+='''
##Diff
fq=root_dir+'/QC_UG/cleandata'
assert not os.system('mkdir Diff_UG')
os.chdir('Diff_UG')
f=open(root_dir+'/Diff_UG/generate_Diff.sh','w')
f.write('perl /PUBLIC/source/RNA/RefRNA/DGE_unigene/DiffExpression/runDiff_analysis_v3.pl -fq %s -fa %s -o %s -group %s -groupname %s -compare %s -venn %s -qc_dir %s ' % (fq,fa,root_dir+'/Diff_UG',group,groupname,compare,venn,qc_dir))
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
code+='''
##Curve
assert not os.system('mkdir Curve_UG')
os.chdir('Curve_UG')
f=open(root_dir+'/Curve_UG/generate_SatCurve.sh','w')
bam=root_dir+'/QC_UG/bam'
genelength=root_dir+'/QC_UG/gene.length'
f.write('perl /PUBLIC/source/RNA/RefRNA/DGE_unigene/runsaturation.pl -bam %s -n %s -gi %s -o %s\\n' % (bam,sample,genelength,root_dir+'/Curve_UG'))
f.write('sh %s' % (root_dir+'/Curve_UG/runSaturationCurve.sh'))
f.close()

assert not os.system('qsub -V -cwd -l vf=10G -l p=1 generate_SatCurve.sh')
os.chdir(root_dir)

'''
open('UG_step2_Diff_SatCurve.py','w').write(code)
#################################################################################################
#for Enrichment and PPI
code='''
import os

fa='%s'
sample='%s'
goann='%s'
groupname='%s'
compare='%s'
root_dir='%s'
''' % (fa,sample,goann,groupname,compare,root_dir)
if argv['species']:
	code+='''
species='%s'
''' % (species)
else:
	code+='''
koann='%s'
''' % (koann)

if argv['ppi_number']:
	code+='''
ppi_number='%s'
ppi_blast='%s'
''' % (ppi_number,ppi_blast)

code+='''
groupnames=groupname.split(',')
compares=compare.split(',')


##GOSeq
assert not os.system('mkdir GOSeq_UG')
os.chdir('GOSeq_UG')
go=open(root_dir+'/GOSeq_UG/runGOSeq.sh','w')
for each in compares:
	temp=each.split(':')
	dir=groupnames[int(temp[0])-1]+'vs'+groupnames[int(temp[1])-1]
	id=root_dir+'/Diff_UG/Diff/'+dir+'/'+dir+'.diffgeneID'
	out=root_dir+'/GOSeq_UG/'+dir
	length=root_dir+'/Diff_UG/Diff/gene.length'
	result=root_dir+'/GOSeq_UG/'+dir+'/'+dir+'.GO_enrichment_result.xls'
	diff=root_dir+'/Diff_UG/Diff/'+dir+'/'+dir+'.diffgene.xls'
	go.write('###############%s#####################\\n' % (dir))
	go.write('mkdir %s\\n' % (root_dir+'/GOSeq_UG/'+dir))
	go.write('perl /PUBLIC/source/RNA/RefRNA/DGE/Enrichment/goseq_graph_v3.pl -i %s -goann %s -n %s -o %s -length %s\\n' % (id,goann,dir,out,length))
	go.write('perl /PUBLIC/source/RNA/RefRNA/DGE/Enrichment/changeGO_up_down.pl %s %s %s\\n' % (result,diff,root_dir+'/GOSeq_UG/'+dir+'/'+dir+'.GO_enrichment_result_up_down.xls'))
	go.write('/PUBLIC/software/public/System/R-2.15.3/bin/Rscript /PUBLIC/source/RNA/RefRNA/DGE/Enrichment/goBar.R %s %s %s\\n' % (result,out,dir))
	go.write('/PUBLIC/software/public/System/R-2.15.3/bin/Rscript /PUBLIC/source/RNA/RefRNA/DGE/Enrichment/goBar2.R %s %s %s\\n' % (root_dir+'/GOSeq_UG/'+dir+'/'+dir+'.GO_enrichment_result_up_down.xls',out,dir))
go.close()
assert not os.system('qsub -V -cwd -l vf=5G -l p=1 runGOSeq.sh')
os.chdir(root_dir)
'''

if argv['genenamefile'] == None:
	code+='''

#genenamefile
assert not os.system('mkdir Blast_UG')
os.chdir('Blast_UG')
f=open(root_dir+'/Blast_UG/runBlast_swissprot.sh','w')
query=root_dir+'/Diff_UG/Diff/diffgene_union.seq'
outdir1=root_dir+'/Blast_UG/Blast_Swissprot/'
out=root_dir+'/Blast_UG/Blast_Swissprot/diffgene_union.seq.blastout'
f.write('echo start blastx\\ndate\\n')
f.write('mkdir %s\\n' % (outdir1))
f.write('/PUBLIC/software/public/Alignment/ncbi-blast-2.2.28+/bin/blastx -query %s -db /PUBLIC/database/Common/SwissProt/uniprot_sprot.fasta -evalue 1e-5 -outfmt 5 -max_target_seqs 1 -num_threads 10 -out %s\\n' % (query,out))
f.write('echo blastx end\\ndate\\n')
f.write('perl /PUBLIC/source/RNA/RefRNA/DGE/scriptdir/extractIDsEVxml.pl %s %s\\n' % (out,outdir1+'/diffgene_union.genenames'))
compare=root_dir+'/Diff_UG/Diff/compare.txt'
indir=root_dir+'/Diff_UG/Diff/'
outdir2=root_dir+'/Diff_UG/Diff/DiffGeneList'
f.write('mkdir %s\\n' % (outdir2))
f.write('perl /PUBLIC/source/RNA/RefRNA/DGE/scriptdir/getdiffGN.up_down_v2.pl %s %s %s %s\\n' % (indir,compare,outdir1+'/diffgene_union.genenames',outdir2))
f.close()
assert not os.system('qsub -V -cwd -l vf=5G -l p=10 runBlast_swissprot.sh')
os.chdir(root_dir)
'''


if argv['species']:
	code+='''
##KOBAS
assert not os.system('mkdir KOBAS_UG')
os.chdir('KOBAS_UG')
ko=open(root_dir+'/KOBAS_UG/runKOBAS.sh','w')
path=open(root_dir+'/KOBAS_UG/runpathway.sh','w')
query=root_dir+'/Diff_UG/Diff/diffgene_union.seq'
blastout=root_dir+'/Blast_UG/KOBAS_blast.xml'
ko.write('###############Run Blast#####################\\n')
ko.write('mkdir %s\\n' % (root_dir+'/Blast_UG/'))
ko.write('perl /PUBLIC/source/RNA/RefRNA/DGE/KEGG_Enrichment/bin/KEGG_step1_blast.pl %s %s %s %s\\n' % (query,species,blastout,root_dir+'/Blast_UG/KOBAS_blast.sh'))
ko.write('sh %s\\n' % (root_dir+'/Blast_UG/KOBAS_blast.sh'))
for each in compares:
	temp=each.split(':')
	dir=groupnames[int(temp[0])-1]+'vs'+groupnames[int(temp[1])-1]
	id=root_dir+'/Diff_UG/Diff/'+dir+'/'+dir+'.diffgeneID'
	out=root_dir+'/KOBAS_UG/'+dir
	ko.write('###############%s#####################\\n' % (dir))
	path.write('echo "###############%s#####################"\\n' % (dir))
	path.write('cd %s\\n' % (out))
	result=root_dir+'/KOBAS_UG/'+dir+'/'+'add.'+dir+'.identify.xls'
	diff=root_dir+'/Diff_UG/Diff/'+dir+'/'+dir+'.diffgene.xls'
	path.write('python /PUBLIC/source/RNA/RefRNA/DGE/KEGG_Enrichment/bin/pathway_annotation_flow_parallel_annotationfault_tolerant3.pyc --table %s --diff %s\\n' % (result,diff))
	path.write('mv %s %s\\n' % ('add.'+dir+'.identify.xls_rendered_html_detail.html',dir+'.html'))
	ko.write('mkdir %s\\n' % (root_dir+'/KOBAS_UG/'+dir))
	script=root_dir+'/KOBAS_UG/'+dir+'/run.sh'
	ko.write('perl /PUBLIC/source/RNA/RefRNA/DGE/KEGG_Enrichment/bin/KEGG_step2_enrich.pl -id %s -out-dir %s -species %s -blast-result %s -sample-names %s>%s\\n' % (id,out,species,blastout,dir,script))
	ko.write('sh %s\\n' % (script))
path.close()
ko.close()
assert not os.system('qsub -V -cwd -l vf=5G -l p=6 runKOBAS.sh')
os.chdir(root_dir)
'''
else:
	code+='''
#KOBAS
assert not os.system('mkdir KOBAS_UG')
os.chdir('KOBAS_UG')
ko=open(root_dir+'/KOBAS_UG/runKOBAS.sh','w')
path=open(root_dir+'/KOBAS_UG/runpathway.sh','w')
for each in compares:
	temp=each.split(':')
	dir=groupnames[int(temp[0])-1]+'vs'+groupnames[int(temp[1])-1]
	id=root_dir+'/Diff_UG/Diff/'+dir+'/'+dir+'.diffgeneID'
	out=root_dir+'/KOBAS_UG/'+dir
	ko.write('###############%s#####################\\n' % (dir))
	path.write('###############%s#####################\\n' % (dir))
	path.write('cd %s\\n' % (out))
	result=root_dir+'/KOBAS_UG/'+dir+'/'+'add.'+dir+'.identify.xls'
	diff=root_dir+'/Diff_UG/Diff/'+dir+'/'+dir+'.diffgene.xls'
	path.write('python /PUBLIC/source/RNA/RefRNA/DGE/Enrichment/pathway_annotation_flow_parallel_annotationfault_tolerant.pyc --table %s --diff %s\\n' % (result,diff))
	path.write('mv %s %s\\n' % ('add.'+dir+'.identify.xls_rendered_html_detail.html',dir+'.html'))
	ko.write('mkdir %s\\n' % (root_dir+'/KOBAS_UG/'+dir))
	script=root_dir+'/KOBAS_UG/'+dir+'/run.sh'
	ko.write('perl /PUBLIC/source/RNA/noRef/KEGG_enrichment/runKEGG_enrich.pl -diff %s -ko %s -out %s -g %s\\n' % (id,koann,out,dir))
	ko.write('mv %s/%s.DEG_KEGG_pathway_enrichment_add.xls %s\\n' % (out,dir,result))
path.close()
ko.close()
assert not os.system('qsub -V -cwd -l vf=5G -l p=6 runKOBAS.sh')
os.chdir(root_dir)
'''

if argv['ppi_number']:
	code+='''
##PPI
assert not os.system('mkdir PPI_UG')
os.chdir('PPI_UG')
ppi=open(root_dir+'/PPI_UG/runPPI.sh','w')
for each in compares:
	temp=each.split(':')
	dir=groupnames[int(temp[0])-1]+'vs'+groupnames[int(temp[1])-1]
	id=root_dir+'/Diff_UG/Diff/'+dir+'/'+dir+'.diffgeneID'
	seq=root_dir+'/Diff_UG/Diff/Diff_Gene_Seq/'+dir+'.diffgene.seq'
	out=root_dir+'/PPI_UG/'+dir
	ppi.write('mkdir %s\\n' % (out))
'''
	if argv['ppi_blast'] == 'y':
		code+='''
	ppi.write('python /PUBLIC/source/RNA/RefRNA/DGE/Enrichment/BLASTX_TO_PPI_v4.py --species %s --fa %s --outdir %s --name %s\\n' % (ppi_number,seq,out,dir))
'''
	else:
		code+='''
	ppi_dir='/PUBLIC/database/RNA/PPI_ref/PPI_99'
	ppi.write('python /PUBLIC/source/RNA/RefRNA/UG/Enrichment/get_my_PPI.py -p %s -g %s -o %s\\n' % (ppi_dir+'/PPI_'+ppi_number+'.txt',id,out+'/'+dir+'.ppi.txt'))
'''
	code+='''
ppi.close()
assert not os.system('qsub -V -cwd -l vf=5G -l p=4 runPPI.sh')
os.chdir(root_dir)
'''
open('UG_step3_Enrichment.py','w').write(code)

#################################################################################################
#for Report

code='''
import re
import os

sample='%s'
groupname='%s'
compare='%s'
project='%s'
root_dir='%s'
''' % (sample,groupname,compare,project,root_dir)

code+='''
samples=sample.split(',')

##md5
dir=root_dir+'/QC_UG'
md5=open(root_dir+'/QC_UG/md5.sh','w')
os.chdir(dir)
md5.write('cd %s/cleandata\\n' % (dir))
md5.write('echo -e \"md5\\tfile\">md5.txt\\n')
for each in samples:
	md5.write('md5sum %s.clean.fq>>md5.txt\\n' % (each))
	md5.write('unlink %s.clean.fq\\n' %(each))
	temp=root_dir+'/QC_UG/'+each+'/clean_data/'+each+'.clean.fq'
	md5.write('gzip %s\\n' % (temp))
	md5.write('ln -s %s.gz %s\\n' % (temp,dir+'/cleandata/'))
md5.close()
assert not os.system('qsub -V -l vf=3g,p=1 -cwd md5.sh')
os.chdir(root_dir)

##path
path=root_dir+'/KOBAS_UG/runpathway.sh'
assert not os.system('sh %s' % (path))

##Report
groupnames=groupname.split(',')
compares=compare.split(',')
compare_names=[]
for each in compares:
	temp=each.split(':')
	compare_names.append(groupnames[int(temp[0])-1]+'vs'+groupnames[int(temp[1])-1])
compare_names=','.join(compare_names)
assert not os.system('mkdir ' + project + '_results_UG')
os.system('sh /PUBLIC/source/RNA/RefRNA/DGE_unigene/UG_result_v2.sh -dir %s -sample %s -compare %s -title %s -results %s' % (root_dir,sample,compare_names,project,root_dir+'/' + project + '_results_UG'))
'''
open('UG_step4_Report.py','w').write(code)
#################################################################################################
#for data release
code='''
import os

fa='%s'
goann='%s'
sample='%s'
geneinfo='%s'
root_dir='%s'
project='%s'
''' % (fa,goann,sample,root_dir+'/Diff_UG/Diff/gene.length',root_dir,project)

if argv['genenamefile']:
	code+='''
genenamefile='%s'
''' % (genenamefile)
else:
	code+='''
genenamefile='%s'
''' % (root_dir+'/Blast_UG/Blast_Swissprot/diffgene_union.genenames')

code+='''
assert not os.system('mkdir data_release')
data_release=open(root_dir+'/data_release/data_release.sh','w')
data_release.write('mkdir %s\\n' % (root_dir+'/data_release/data_give'))
dir=root_dir+'/data_release/data_give'
data_release.write('echo "###############prepare SuppFiles#####################"\\n')
data_release.write('mkdir %s\\n' % (root_dir+'/'+project+'_results_UG/results/0.SuppFiles'))
SuppFiles_dir=root_dir+'/'+project+'_results_UG/results/0.SuppFiles'
data_release.write('cp /PUBLIC/source/RNA/RefRNA/TransRef/data_release/SuppFiles.README.txt %s\\n' % (SuppFiles_dir))
data_release.write('cp %s %s\\n' % (genenamefile,SuppFiles_dir+'/gene.description.xls'))
data_release.write('perl /PUBLIC/source/RNA/RefRNA/TransRef/data_release/GO_class.pl %s /PUBLIC/software/RNA/GOseq/gene_ontology.1_2.obo.ab %s\\n' % (goann,SuppFiles_dir+'/gene.goannot.xls'))
data_release.write('sed -e \"1i\Gene_ID\\tLength\" %s>%s\\n' % (geneinfo,SuppFiles_dir+'/gene.info.xls'))
data_release.write('cp %s %s\\n' % (fa,SuppFiles_dir+'/gene.fasta'))
data_release.write('gzip %s\\n' % (SuppFiles_dir+'/gene.fasta'))
data_release.write('echo "###############copy NovoQuery#####################"\\n')
data_release.write('cp /PUBLIC/source/RNA/RefRNA/TransRef/data_release/NovoQuery.exe %s\\n' % (root_dir+'/'+project+'_results_UG/'))
data_release.write('cp /PUBLIC/source/RNA/RefRNA/TransRef/data_release/NovoQuery_manual.pdf %s\\n' % (root_dir+'/'+project+'_results_UG/'))
data_release.write('echo "###############tar Report#####################"\\n')
data_release.write('cd %s\\n' % (root_dir))
data_release.write('tar -cvzf %s %s\\n' % ('data_release/data_give/' + project + '_UG.tar.gz',project+'_results_UG/'))
data_release.write('cd %s\\n' % (root_dir))
data_release.write('echo "###############prepare clean data#####################"\\n')
data_release.write('ln -s %s %s\\n' % (root_dir+'/QC_UG/cleandata',dir+'/cleandata'))
data_release.write('echo "###############prepare raw data#####################"\\n')
data_release.write('mkdir %s\\n' % (dir+'/rawdata'))

samples=sample.split(',')
for eachsample in samples:
        data_release.write('ln -s %s %s\\n' % (root_dir+'/QC_UG/'+eachsample+'/'+eachsample+'.fq.gz',dir+'/rawdata/'))

data_release.write('echo "###############prepare IGV_data#####################"\\n')
data_release.write('mkdir %s\\n' % (dir+'/IGV_data'))
data_release.write('cp /PUBLIC/source/RNA/RefRNA/TransRef/data_release/IGV_data.README.txt %s\\n' % (dir+'/IGV_data'))
data_release.write('ln -s %s %s\\n' % (root_dir+'/QC_UG/bam',dir+'/IGV_data/bam'))
data_release.close()

os.chdir('data_release')
assert not os.system('qsub -V -cwd -l vf=1G -l p=1 data_release.sh')
os.chdir(root_dir)
'''

open('UG_step5_data_release.py','w').write(code)
#################################################################################################
#for byebye
byebye=open(root_dir+'/byebye.sh','w')
byebye.write('sh /PUBLIC/source/RNA/RefRNA/TransRef/Pipeline/byebye_v2.1.sh -dir %s -sample %s -pipline UG\n' % (root_dir,sample))
byebye.close()

