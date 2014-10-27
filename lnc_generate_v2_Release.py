import sys
import os
import os.path
import re
import argparse
root_dir=os.getcwd()
parser = argparse.ArgumentParser(description="lncRNA pipline v1.0")
parser.add_argument('--project',help='project name, maybe same with the name of root dir, which will be displayed in the final report title, [REQUIRED]',required=True)
parser.add_argument('--sample',help="sample names(sample1,sample2,...)  warning: order sensitive!",required=True)
parser.add_argument('--fq',help="the original directory of the raw fastq reads, [REQUIRED]",default=None)
parser.add_argument('--mapfile',help="mapfile files (mapfile1,mapfile2.The first column of mapfile is the NH number ,the second column of mapfile is the sample name)",default=None)
parser.add_argument('--number',help="chromosome number, [REQUIRED for density plot]",required=True)
parser.add_argument('--raw_dir',help="the original directory of the raw fastq reads keep order in line with mapfile files (raw_dir1,raw_dir2)",default=None)
parser.add_argument('--ad',help="the original directory of the adapter list, ",default=None)
parser.add_argument('--generate_adapter',help='whether to generate the adpepter list y or n',choices=['y','n'],default='n')
parser.add_argument('--index',help='the index number of the fastq files,if you want to generate the adapter list,the parameter will be necessary,diffrent samples split by ,warning: order sensitive!',default=None)
parser.add_argument('--length',help="the length of sequenced reads,[defalut=100]",default='100')
parser.add_argument('--fa',help="the reference FASTA file, [REQUIRED for mapping]",required=True)
parser.add_argument('--gtf',help="the annotation GTF file, [REQUIRED for mapping]",required=True)
parser.add_argument('--group',help="sample classification, e.g. sample1/sample2,sample3, [REQUIRED]",default=None)
parser.add_argument('--groupname',help="group names, [group1,group2,... ]",default=None)
parser.add_argument('--compare',help="group comparison strategy 1:2,1:3,..., [REQUIRED]",default=None)
parser.add_argument('--venn',help="venn and cluster mode, defult is all groups, 1:2_1:3,1:3_2:3,... ",default=None)
parser.add_argument('--goann',help="gene to GO annotations file, [REQUIRED]",default=None)
parser.add_argument('--species',help="abbreviation for species, for kegg pathway. (note: all species: /PROJ/RNA/share/software/kobas2.0-data-20120208/seq_pepi_v2), [defalut=ko]",default='kaas')
parser.add_argument('--ppi_number',help="species code, (ref file: /PROJ/RNA/share/database/string_ppi/species.v9.0.txt)",default=None)
parser.add_argument('--ppi_blast',help="whether to run blast to get the protein-protein interactions",choices=['y','n'],default='y')
parser.add_argument('--genenamefile',help="genenamefile, 1st col is geneID, 2nd col is genename",default=None)
parser.add_argument('--maf',help="maffile",default=None)
parser.add_argument('--ss',help="strand specific, [REQUIRED]",choices=['no','yes','reverse'],default='reverse')
parser.add_argument('--abbr',help="abbreviation",default=True)
parser.add_argument('--CPC',help="whether use CPC",choices=['yes','no'],default='yes')
parser.add_argument('--CNCI',help="whether use CNCI",choices=['yes','no'],default='yes')
parser.add_argument('--PFAM',help="whether use PFAM",choices=['yes','no'],default='yes')
parser.add_argument('--phyloCSF',help="whether use phyloCSF",choices=['yes','no'],default='yes')
parser.add_argument('--gff3',help="gff or gtf with CDS,old gff",required=True)
parser.add_argument('--fa3',help="old fasta file",default=None)
parser.add_argument('--DEXseq',help="whether to do DEXseq",choices=['yes','no'],required=True)
parser.add_argument('--type',help="lncRNA type selected to analysis",choices=['u','x','i'],default='u,i,x')
parser.add_argument('--results',help="generate results",choices=['y','n'],default='y')
argv = vars(parser.parse_args())
project=argv['project'].strip()
lnctype=argv['type'].strip()
phyloCSF=argv['phyloCSF'].strip()
abbr=argv['abbr'].strip()
CPC=argv['CPC'].strip()
CNCI=argv['CNCI'].strip()
PFAM=argv['PFAM'].strip()
gff3=argv['gff3'].strip()
#fa3=argv['fa3'].strip()
DEXseq=argv['DEXseq'].strip()
gen_results=argv['results'].strip()
display = open(root_dir + '/' + 'lncRNA_command.txt','w')
display.write('project: %s\n' % (project))
if phyloCSF=='yes' or phyloCSF=='y':
    maf=argv['maf'].strip()
    display.write('maf file: %s\n' % (maf))
display.write('the abbreviation for species : %s\n' % (abbr))
display.write('whether use CPC : %s\n' % (CPC))
display.write('whether use CNCI : %s\n' % (CNCI))
display.write('whether use PFAM : %s\n' % (PFAM))
display.write('whether use phyloCSF : %s\n' % (phyloCSF))
sample=argv['sample'].strip()
samples=sample.split(',')
samples_tmp=list(set(samples))
assert len(samples)==len(samples_tmp)
display.write('sample: %s\n' % (sample))
mapfiles=[each.strip() for each in argv['mapfile'].strip().split(',') if each.strip() != '']
mapfile=' '.join(mapfiles)
assert not os.system('cat %s |sort -u |sed \'/^$/d\' >%s' % (mapfile,root_dir+'/libraryID'))
if argv['fq']:
        fq = argv['fq'].strip()
        fq = os.path.abspath(fq)
else:
        os.system('mkdir raw_data')
        assert not os.system('perl /PUBLIC/source/RNA/RefRNA/ln_raw_data.pl %s %s pe raw_data' %(argv['mapfile'],argv['raw_dir']))
        fq = root_dir + '/raw_data'
for each in samples:
        fq_tmp1=fq+'/'+each+'_1.fq.gz'
        fq_tmp2=fq+'/'+each+'_2.fq.gz'
        assert os.path.isfile(fq_tmp1)
        assert os.path.isfile(fq_tmp2)
display.write('fq: %s\n' % (fq))
display.write('adapter: ')
if argv['ad']:
        ad = argv['ad'].strip()
        ad = os.path.abspath(ad)
        for each in samples:
                ad_tmp1=ad+'/'+each+'_1.adapter.list.gz'
                ad_tmp2=ad+'/'+each+'_2.adapter.list.gz'
                assert os.path.isfile(ad_tmp1)
                assert os.path.isfile(ad_tmp2)
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
if argv['group'] == None:
	groups=samples
	groups_iter=samples
	group=sample
	flag_repeat=False
else:
	group=argv['group'].strip()
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
	groupname=argv['groupname'].strip()
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
venn_plot='no'
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
	venn_plot='yes'
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
#--------------------------------------------------------------------------------------------------
goann = argv['goann'].strip()
goann = os.path.abspath(goann)
assert os.path.isfile(goann)
display.write('goann: %s\n' % (goann))
if argv['species']:
	species = argv['species'].strip()
else:
	species = 'kaas'
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

if ss == 'no':
    lib='fr-unstranded'
elif ss == 'reverse':
    lib = 'fr-firststrand'
else:
    lib = 'fr-secondstrand'

flag_uniform = False
if gen_results == 'y':
	flag_results = True
else:
	flag_results = False
################################################################################################################
###for QC
#######################################

code='''
import os
import os.path
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
assert not os.system('mkdir 1.QC')
os.chdir('1.QC')
f=open(root_dir+'/1.QC/generate_QC.sh','w')
f.write('awk \\'{if($3==\\"exon\\"){print $0}}\\' %s > %s\\n' % (gtf,root_dir+'/1.QC/exon.gtf'))
f.write('msort -k mf1 -k nf4 %s > %s\\n' % (root_dir+'/1.QC/exon.gtf',root_dir+'/1.QC/sorted.gtf'))
f.write('perl /PUBLIC/source/RNA/QC/QC_v2/gtf2bed.pl %s %s\\n' % (root_dir+'/1.QC/sorted.gtf',root_dir+'/1.QC/sorted.bed'))
f.write('perl /PUBLIC/source/RNA/lncRNA/lncRNA_v1/QC/lnc_runQC_v2.pl  -fq %s -se-pe pe -n %s -o %s -spe %s -R %s -G %s -bed %s ' %(fq,sample,root_dir+'/1.QC',ss,fa,gtf,root_dir+'/1.QC/sorted.bed'))
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
assert not os.system('sh %s/generate_QC.sh' %(root_dir+'/1.QC'))

samples=sample.split(',')
for eachsample in samples:
    os.chdir(root_dir+'/1.QC/'+eachsample)
    assert not os.system('qsub -V -l vf=15g,p=4 -cwd %s' % (eachsample+'_QC.sh'))
os.chdir(root_dir)
'''
open('lnc_step1_QC.py','w').write(code)

###########################################################################################################
######QCreport generate
###########################################################################################################
code='''
import os
sample='%s'
project='%s'
root_dir='%s'
''' % (sample,project,root_dir)

code+='''

assert not os.system('mkdir 1.QC/QCreport')
os.system('sh /PUBLIC/source/RNA/lncRNA/lncRNA_v1/QC/lnc_QCreport.sh -dir %s -sample %s -title %s -results %s' %(root_dir+'/1.QC',sample,project,root_dir+'/1.QC/QCreport'))
'''
open('lnc_step1_QCreport.py','w').write(code)




############################################################################################################
###
###for Assembly and RNAseq Advanced QC
bam = root_dir+"/1.QC/bam"
sam = root_dir+"/1.QC/sam"


code='''
import os
import os.path
import re
import argparse
import sys
import misopy.cluster_utils as cluster_utils
from os.path    import isfile, join
from os         import listdir

parser = argparse.ArgumentParser(description='run assemble & advanced QC and resource list required')
parser.add_argument('--n',default='0',help="the split number")
parser.add_argument('--h',default='3',help="the fat node for qsub")
argv=vars(parser.parse_args())
h=argv['h']
n=argv['n']
n=int(n)
sample='%s'
fa='%s'
gtf='%s'
bam='%s'
sam='%s'
group='%s'
groupname='%s'
root_dir='%s'
lib='%s'
length='%s'
number='%s'
abbr='%s'
''' % (sample,fa,gtf,bam,sam,group,groupname,root_dir,lib,length,number,abbr)

code+='''
#for assemble
#cufflinks

def loadfile ( path, suffix ):
    allfiles = [ join(path,f) for f in listdir(path) if join(path,f).endswith(suffix) ]
    return allfiles

def generate_qsub(sh, script_dir, vf, p):
    cluster_script = 'qsub -cwd -l vf=%dG,p=%d -o \"%s\" -e \"%s\" -V %s' % (vf, p, script_dir, script_dir, sh)
    return cluster_script

samples=sample.split(',')

assert not os.system('mkdir 2.lnc_assemble') 
os.chdir('2.lnc_assemble')

##cufflinks
assert not os.system('mkdir cufflinks')
cufflinks_dir=os.path.abspath("cufflinks")

f=open(root_dir+'/2.lnc_assemble/generate_cufflinks.sh','w')
f.write('perl /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Assemble/cufflinks/cufflinks_assembly_quatification.pl -i %s -p %s -o %s -n 4 -lib %s \\n' % (bam,sample,cufflinks_dir,lib))
f.close()
os.system('sh %s' % (root_dir+'/2.lnc_assemble/generate_cufflinks.sh'))

for eachsample in samples:
        os.chdir(root_dir+'/2.lnc_assemble/cufflinks/'+eachsample)
        assert not os.system('qsub -V -l vf=5g,p=4 -cwd %s' % (eachsample+'.sh'))
os.chdir(root_dir)

##scripture
os.chdir('2.lnc_assemble')
assert not os.system('mkdir scripture')
scripture_dir=os.path.abspath("scripture")
os.chdir(scripture_dir)
os.system('mkdir data_prepare')
data_prepare_dir=os.path.abspath("data_prepare")
f=open(scripture_dir+"/run_scripture_split.sh",'w')
f.write('python  /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Assemble/scripture/bin/fa_split_length_stat.py  --inputfasta %s --outputlength %s --subfa %s \\n' % (fa,data_prepare_dir+"/genome_length.txt",data_prepare_dir+"/subfa"))
f.close()
scripture_split_Qsub = generate_qsub(scripture_dir+"/run_scripture_split.sh", scripture_dir+'/data_prepare', 1, 1)
scripture_split_Job = cluster_utils.launch_job(scripture_split_Qsub,cmd_name = 'qsub')
cluster_utils.wait_on_jobs([scripture_split_Job], cluster_cmd = 'qsub')

gener_script_jobs=[]
for eachsample in samples:
	os.chdir(scripture_dir)
	assert not os.system('mkdir %s' % (eachsample))
	sample_dir=os.path.abspath(eachsample)
	os.system('mkdir %s' % (sample_dir+"/temp"))
	f=open(sample_dir+"/"+"generatescript_"+eachsample+".sh",'w')
	f.write('python /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Assemble/scripture/bin/bam_stats.py  %s >%s \\n' % (root_dir+"/1.QC/bam/"+eachsample+".bam",sample_dir+"/temp/"+eachsample+"_bam_stat.xls"))
	#f.write('wait\\n')
	f.write('python /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Assemble/scripture/bin/scripture_script_generate.py --bam_stat %s --output %s --genome_length %s --subfa %s --bam %s --sample %s\\n' % (sample_dir+"/temp/"+eachsample+"_bam_stat.xls",sample_dir,data_prepare_dir+"/genome_length.txt",data_prepare_dir+"/subfa",bam+"/"+eachsample+".bam",eachsample))
	f.write('python /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Assemble/scripture/bin/assemble_split.py %s %s\\n' % (sample_dir+"/script/"+eachsample+".sh",sample_dir+"/script/split"))
	f.close()
	gener_script_Qsub = generate_qsub(sample_dir+"/"+"generatescript_"+eachsample+".sh", sample_dir, 1, 1)
	gener_script_jobs.append(cluster_utils.launch_job(gener_script_Qsub,cmd_name = 'qsub'))
cluster_utils.wait_on_jobs(gener_script_jobs, cluster_cmd = 'qsub')

os.chdir(scripture_dir)
for eachsample in samples:
	scripture_sh_list = loadfile (scripture_dir+'/'+eachsample+'/script/split','.sh')
	for sub_sh_file in scripture_sh_list:
		os.system('qsub -cwd -l vf=15G -V -l p=1 %s' % (sub_sh_file))

os.chdir(root_dir)

assert not os.system('mkdir 3.RNA_advanced_QC')
Curve_dir=os.path.abspath('3.RNA_advanced_QC')

'''
if flag_uniform:
	code+='''
os.chdir('3.RNA_advanced_QC')
assert not os.system('mkdir UniformDistribution')
UniformDistribution_dir=os.path.abspath('UniformDistribution')
f=open(Curve_dir+"/run_MC.sh",'w')
samples=sample.split(',')
for eachsample  in samples:
	os.system('mkdir %s\\n' % ( UniformDistribution_dir+"/"+eachsample))
	f.write('cd %s\\n' % (UniformDistribution_dir+"/"+eachsample))
	os.chdir(UniformDistribution_dir+"/"+eachsample)
	os.system('ln -s %s .' % (bam+"/"+eachsample+".bam"))
	os.system('ln -s %s %s' % (fa,abbr+".fa"))
	os.system('ln -s %s %s' % (gtf,abbr+".gtf"))
	f.write('perl /PUBLIC/source/RNA/QC/QCscripts/runRNAmatrix.pl -I %s -R %s -G %s\\n' % (UniformDistribution_dir+"/"+eachsample+"/"+eachsample+".bam",UniformDistribution_dir+"/"+eachsample+"/"+abbr+".fa",UniformDistribution_dir+"/"+eachsample+"/"+abbr+".gtf"))
	f.write('sh %s/run_RNAmatrix.sh\\n' % (UniformDistribution_dir+"/"+eachsample))
	os.chdir(UniformDistribution_dir)
f.close()
os.chdir(Curve_dir)
cluster_num=re.search(r'([1-9])',h).group(1)
P_mem='mem'+str(cluster_num)
q_mem=P_mem+'.q'
if cluster_num==str(3) or cluster_num==str(4):
    P_mem='mem3'
    q_mem='mem3.q'
assert not os.system('/opt/gridengine/bin/lx26-amd64/qsub -cwd -l vf=30G -V -l p=5  -P %s -q %s run_MC.sh' % (P_mem,q_mem))
'''
code+='''
os.chdir(Curve_dir)
assert not os.system('mkdir Density')
Density_dir=os.path.abspath('Density')
f=open(Curve_dir+"/run_density.sh",'w')
for eachsample in samples:
	os.system('mkdir %s\\n' % (Density_dir+"/"+eachsample))
	f.write('perl /PUBLIC/source/RNA/RefRNA/DGE/Curve/reads_density.pl -bam %s -fa %s -n %s -r %s -sample %s -o %s\\n' % (bam+"/"+eachsample+".bam",fa,number,length,eachsample,Density_dir+"/"+eachsample))
f.close()
os.chdir(Curve_dir)
assert not os.system('qsub -cwd -l vf=5G -l p=2 -V %s' % (Curve_dir+"/run_density.sh")) 


'''
open('lnc_step2_Assembly_AdvancedQC.py','w').write(code)


##################################################################################################
#####mRNA cufflinks & coding_noncoding_filt  
##########################
cufflinks_dir=os.path.abspath(root_dir+"/2.lnc_assemble/cufflinks")
scripture_dir=os.path.abspath(root_dir+"/2.lnc_assemble/scripture")

code='''
import os
from optparse   import OptionParser
from sys        import stderr
from os.path    import splitext
from os         import chdir
from os         import listdir
from os         import system
from os         import getcwd
from os.path    import isfile, join
from subprocess import Popen
from subprocess import PIPE
import misopy.cluster_utils as cluster_utils


RUBY_PATH       = "/PUBLIC/software/public/ruby-2.1.1/bin/ruby"
RUBY_SCRIPT     = "/PUBLIC/software/public/Annotation/kaas_sa2/script/multi_cpu.rb"

sample='%s'
fa='%s'
gtf='%s'
bam='%s'
sam='%s'
group='%s'
groupname='%s'
root_dir='%s'
abbr='%s'
CPC='%s'
CNCI='%s'
PFAM='%s'
phyloCSF='%s'
scripture_dir='%s'
cufflinks_dir='%s'
gff3='%s'
lnctype='%s'
lib='%s'
''' % (sample,fa,gtf,bam,sam,group,groupname,root_dir,abbr,CPC,CNCI,PFAM,phyloCSF,scripture_dir,cufflinks_dir,gff3,lnctype,lib)

if argv['maf']:
    code+='''
maf='%s'
''' % (maf)
code+='''

###########################mRNA Cufflinks########################################
def generate_qsub(sh, script_dir, vf, p):
    cluster_script = 'qsub -cwd -l vf=%dG,p=%d -o \"%s\" -e \"%s\" -V %s' % (vf, p, script_dir, script_dir, sh)
    return cluster_script

os.chdir(root_dir)
assert not os.system('mkdir 4.mRNA_Cufflinks')
os.chdir('4.mRNA_Cufflinks')
f=open(root_dir+'/4.mRNA_Cufflinks/generate_mRNA_Cufflinks.sh','w')
f.write('awk \\'{if (/gene_type \\"mRNA\\"|gene_type \\"protein_coding\\"/) {print $0}}\\' %s > %s\\n' % (gtf,root_dir+'/4.mRNA_Cufflinks/mRNA.gtf'))
f.write('perl /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Assemble/cufflinks/cufflinks_assembly_quatification.pl -gtf %s -i %s -p %s -o %s -n 4 -lib %s \\n' % (root_dir+'/4.mRNA_Cufflinks/mRNA.gtf',bam,sample,root_dir+'/4.mRNA_Cufflinks/Cufflinks',lib))
f.close()
os.system('sh %s' % (root_dir+'/4.mRNA_Cufflinks/generate_mRNA_Cufflinks.sh'))
samples=sample.split(',')
mRNA_cufflinks_job=[]
for eachsample in samples:
        os.chdir(root_dir+'/4.mRNA_Cufflinks/Cufflinks/'+eachsample)
        mRNA_cufflinks_Qsub = generate_qsub(eachsample+'.sh',root_dir+'/4.mRNA_Cufflinks/Cufflinks/'+eachsample, 5, 4)
        mRNA_cufflinks_job.append( cluster_utils.launch_job(mRNA_cufflinks_Qsub,cmd_name = 'qsub'))
os.chdir(root_dir)

#########################coding-noncoding filter####################################
os.system('mkdir 5.coding_noncoding_filt')
coding_noncoding_filt_dir=os.path.abspath('5.coding_noncoding_filt')
os.chdir('5.coding_noncoding_filt')
os.system('mkdir step1')
os.system('mkdir step2')
os.system('mkdir step3')
os.system('mkdir step4')
os.system('mkdir step5')
os.system('mkdir step6')
os.system('mkdir step7_compare')
step1_dir=os.path.abspath("step1")
step2_dir=os.path.abspath("step2")
step3_dir=os.path.abspath("step3")
step4_dir=os.path.abspath("step4")
step5_dir=os.path.abspath("step5")
step6_dir=os.path.abspath("step6")
step7_dir=os.path.abspath("step7_compare")


cufflinks=open(root_dir+'/5.coding_noncoding_filt/run_cufflinks_prepare.sh','w')
cufflinks.write('cd %s\\n ' % (coding_noncoding_filt_dir))
cufflinks.write('ln -s %s \\n' % (fa))
cufflinks.write('ln -s %s \\n' % (gtf))
cufflinks.write('ln -s %s \\n' % (maf))
cufflinks.write('perl /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Filter/bin/GTF_mRNA_others.pl %s %s/mRNA.gtf %s/others.gtf\\n' % (gtf,coding_noncoding_filt_dir,coding_noncoding_filt_dir))
cufflinks.close()
CufflinksPre_Qsub = generate_qsub(root_dir+'/5.coding_noncoding_filt/run_cufflinks_prepare.sh', root_dir+'/5.coding_noncoding_filt', 1, 1)
CufflinksPre_Job = cluster_utils.launch_job(CufflinksPre_Qsub,cmd_name = 'qsub')
cluster_utils.wait_on_jobs([CufflinksPre_Job], cluster_cmd = 'qsub')

samples=sample.split(',')
os.chdir(step1_dir)
for eachsample in samples:
	os.system('mkdir %s' % (eachsample))
	eachsample_dir=os.path.abspath(eachsample)
	os.system('mkdir %s' % ( eachsample+"/cufflinks"))
	os.system('mkdir %s' % ( eachsample+"/scripture"))
	os.system('mkdir %s' % (eachsample+"/scripture/quatification"))
cufflinks=open(step1_dir+'/run_cufflinks.sh','w')
cufflinks.write('perl /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Filter/bin/cufflinks_assembly_quatification.pl -i %s -p %s -o %s -n 4 -lib %s -cufflinks %s  -scripture %s\\n' % (bam,sample,step1_dir,lib,cufflinks_dir,scripture_dir))
cufflinks.close()
os.system('sh %s' % (step1_dir+'/run_cufflinks.sh'))

Cufflinks_all_Job=[]
for eachsample in samples:
	os.chdir(step1_dir+'/'+eachsample)
	Cufflinks_Qsub = generate_qsub(eachsample+'.sh',step1_dir+'/'+eachsample, 5, 4)
	Cufflinks_all_Job.append( cluster_utils.launch_job(Cufflinks_Qsub,cmd_name = 'qsub'))
	os.chdir(step1_dir)
cluster_utils.wait_on_jobs(Cufflinks_all_Job, cluster_cmd = 'qsub')

os.chdir(step1_dir)
cuffcompare=open(root_dir+'/5.coding_noncoding_filt/run_cuffcompare.sh','w')
for eachsample in samples:
	eachsample_dir=os.path.abspath(eachsample)
	quatification_dir=os.path.abspath(eachsample+"/scripture/quatification")
	cuffcompare.write('perl /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Filter/bin/GTF_multi_single_exon_length.pl %s %s %s \\n' % (quatification_dir+"/transcripts.gtf",eachsample+".scripture",eachsample_dir+"/scripture"))

os.chdir(step2_dir)
for eachsample in samples:
	os.system('mkdir %s ' % (eachsample))
	eachsample_dir=os.path.abspath(eachsample)
	os.system('mkdir %s ' % (eachsample+"/cufflinks"))
	os.system('mkdir %s ' % (eachsample+"/scripture"))
	cuffcompare.write('perl /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Filter/bin/filter_cov.pl %s %s \\n' % (step1_dir+"/"+eachsample+"/cufflinks/"+eachsample+".cufflinks.step1_length_multi_exon.gtf",eachsample_dir+"/cufflinks/"+eachsample+".cufflinks.step2_cov.gtf"))
	cuffcompare.write('perl /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Filter/bin/filter_cov.pl %s %s \\n' % (step1_dir+"/"+eachsample+"/scripture/"+eachsample+".scripture.step1_length_multi_exon.gtf",eachsample_dir+"/scripture/"+eachsample+".scripture.step2_cov.gtf"))

os.chdir(step3_dir)
cuffcompare.write('ls %s/*step2_cov.gtf >%s \\n' % (step2_dir+"/*/*",step3_dir+"/cufflinks_list.txt"))
cuffcompare.write('/PUBLIC/software/RNA/cufflinks-2.1.1/cuffcompare -R  -o %s/all -i %s\\n' % (step3_dir,step3_dir+"/cufflinks_list.txt"))
cuffcompare.write('python /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Filter/bin/get_reconstruction_gtf.py %s/all.tracking %s/all.combined.gtf >%s/all.step3_reconstruction.gtf \\n' % (step3_dir,step3_dir,step3_dir))

os.chdir(step4_dir)
cuffcompare.write('/PUBLIC/software/RNA/cufflinks-2.1.1/cuffcompare -r %s -R %s -o %s/cuffcompare_others\\n' % (coding_noncoding_filt_dir+"/others.gtf",step3_dir+"/all.step3_reconstruction.gtf",step4_dir))
cuffcompare.write('perl /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Filter/bin/GTF_delete_c_j_eq.pl -in_gtf %s -compared_gtf %s -out_gtf %s -transcript_class_code %s \\n' % (step3_dir+"/all.step3_reconstruction.gtf",step4_dir+"/cuffcompare_others.combined.gtf",step4_dir+"/cuffcompare_others.step4_delete_c_j_eq.gtf",step4_dir+"/cuffcompare_others.step4_class_code.list.txt"))

os.chdir(step5_dir)
cuffcompare.write('/PUBLIC/software/RNA/cufflinks-2.1.1/cuffcompare -r %s -R %s -o %s\\n' % (coding_noncoding_filt_dir+"/mRNA.gtf",step4_dir+"/cuffcompare_others.step4_delete_c_j_eq.gtf",step5_dir+"/cuffcompare_mRNA"))
cuffcompare.write('perl /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Filter/bin/GTF_extract_classcode.pl -in_gtf %s -compared_gtf %s -extract_class_code %s -out_gtf %s -transcript_class_code %s\\n' % (step4_dir+"/cuffcompare_others.step4_delete_c_j_eq.gtf",step5_dir+"/cuffcompare_mRNA.combined.gtf",lnctype,step5_dir+"/cuffcompare_mRNA.step5_extract_class_code.gtf",step5_dir+"/cuffcompare_mRNA.step5_class_code.list.txt"))
cuffcompare.write('mkdir %s/result\\n' % (step5_dir))
cuffcompare.write('cp %s %s\\n' % (step5_dir+"/cuffcompare_mRNA.step5_extract_class_code.gtf",step5_dir+"/result/step5.gtf"))
cuffcompare.write('perl /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Enrichment/bin/extractTranscriptsfromFA.pl %s %s %s\\n' % (step5_dir+"/result/step5.gtf",fa,step5_dir+"/result/step5.fa"))
cuffcompare.write('python /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Filter/bin/picture/step1-step5_stat.py %s %s >%s\\n' % (sample,root_dir,step5_dir+"/result/step_num_stat.xls"))
cuffcompare.write('Rscript /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Filter/bin/picture/step1-step5_stat_v2.R %s %s\\n' % (step5_dir+"/result/step_num_stat.xls",step5_dir+"/result"))
cuffcompare.write('Rscript /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Filter/bin/picture/class_code_v2.R %s %s\\n' % (step5_dir+"/cuffcompare_mRNA.step5_class_code.list.txt",step5_dir+"/result"))
cuffcompare.close()
Cuffcompare_Qsub = generate_qsub(root_dir+'/5.coding_noncoding_filt/run_cuffcompare.sh', root_dir+'/5.coding_noncoding_filt/step5',4,2)
Cuffcompare_Job = cluster_utils.launch_job(Cuffcompare_Qsub,cmd_name = 'qsub')
cluster_utils.wait_on_jobs([Cuffcompare_Job], cluster_cmd = 'qsub')

os.chdir(step6_dir)
os.system('mkdir CNCI')
os.system('mkdir CPC')
os.system('mkdir PFAM')
os.system('mkdir phyloCSF')
os.system('mkdir venn')
CNCI_dir=os.path.abspath('CNCI')
CPC_dir=os.path.abspath('CPC')
PFAM_dir=os.path.abspath('PFAM')
phyloCSF_dir=os.path.abspath('phyloCSF')
venn_dir=os.path.abspath('venn')

venn_plot=open(venn_dir+"/run_venn.sh",'w')
venn_plot.write('perl /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Filter/venn/CPC_CNCI_PFAM_CSF.venn.pl -gtf %s -outdir %s  ' % (step5_dir+"/result/step5.gtf",venn_dir))

Coding_filter=[]
if CNCI=='yes' or CNCI=='y':
	cnci=open(CNCI_dir+"/run_CNCI.sh",'w')
	cnci.write('perl /PUBLIC/source/RNA/lncRNA/Filter/CNCI/CNCI.pl -f %s -p 5 -l /PUBLIC/software/RNA/CNCI_package/libsvm-3.0 -o %s -b\\n' % (step5_dir+"/result/step5.fa",CNCI_dir))
	cnci.write('sed -i -e \\'s/>//g\\' -e \\'s/score: //g\\' -e \\'s/start: //g\\' -e \\'s/stop: //g\\'  %s\\n' % (CNCI_dir+"/cnci.result.file"))
	cnci.write('cat /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Filter/CNCI/bin/name.xls %s >%s\\n' % (CNCI_dir+"/cnci.result.file",CNCI_dir+"/CNCI.result.xls"))
	cnci.write('awk \\'{if($2 == \\"noncoding\\"){print $1}}\\' %s > %s\\n' % (CNCI_dir+"/CNCI.result.xls",CNCI_dir+"/CNCI.id.txt")) 
	cnci.close()
	CNCI_Qsub = generate_qsub(CNCI_dir+"/run_CNCI.sh", CNCI_dir,2,4)
	CNCI_Job = cluster_utils.launch_job(CNCI_Qsub,cmd_name = 'qsub')
	Coding_filter.append(CNCI_Job)
	venn_plot.write('-CNCI %s  ' % (CNCI_dir+"/CNCI.id.txt"))
if CPC=='yes' or CPC=='y':
	cpc=open(CPC_dir+"/run_CPC.sh",'w')
	cpc.write('perl /PUBLIC/source/RNA/noRef/ORF_predictor/runCPC.pl -seq %s -db eu -strand plus -outdir %s\\n' % (step5_dir+"/result/step5.fa",CPC_dir))
	cpc.write('cat /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Filter/CPC/bin/name.xls %s > %s \\n' % (CPC_dir+"/step5.cpc.txt",CPC_dir+"/CPC.result.xls"))
	cpc.write('awk \\'{if($3 == \\"noncoding\\"){print $1}}\\' %s > %s\\n' %  (CPC_dir+"/CPC.result.xls",CPC_dir+"/CPC.id.txt"))
	cpc.close()
	CPC_Qsub = generate_qsub(CPC_dir+"/run_CPC.sh", CPC_dir,6,4)
	CPC_Job = cluster_utils.launch_job(CPC_Qsub,cmd_name = 'qsub')
	Coding_filter.append(CPC_Job)
	venn_plot.write('-CPC %s  ' % (CPC_dir+"/CPC.id.txt"))
if PFAM=='yes' or PFAM=='y':
	pfam=open(PFAM_dir+"/run_PFAM.sh",'w')
	pfam.write('perl /PUBLIC/source/RNA/lncRNA/Filter/PFAM/seq2protein.pl %s >%s\\n' % (step5_dir+"/result/step5.fa",PFAM_dir+"/lnc.protein.fa"))
	pfam.write('perl /PUBLIC/source/RNA/lncRNA/Filter/PFAM/pfam_scan.pl -fasta %s -dir /PUBLIC/database/Common/PFAM  -outfile %s -pfamB -cpu 3\\n' % (PFAM_dir+"/lnc.protein.fa",PFAM_dir+"/lnc.pfam_scan.out"))
	pfam.write("less %s |sed -e \\'/^#/d\\' -e \\'/^$/d\\' |awk \\'OFS=\\"\\\\t\\"{print ($1,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15)}\\' > %s \\n" % (PFAM_dir+"/lnc.pfam_scan.out",PFAM_dir+"/lnc.pfam_scan.out.temp"))	
	pfam.write('cat /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Filter/PFAM/bin/name.xls %s > %s\\n' % (PFAM_dir+"/lnc.pfam_scan.out.temp",PFAM_dir+"/PFAM.result.xls"))
	pfam.write('python /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Filter/PFAM/bin/select_pfam_id.py %s %s %s\\n' % (step5_dir+"/result/step5.gtf",PFAM_dir+"/lnc.pfam_scan.out.temp",PFAM_dir+"/PFAM.id.txt"))
	pfam.write('rm %s\\n' % (PFAM_dir+"/lnc.pfam_scan.out.temp"))
	pfam.close()
	Pfam_Qsub = generate_qsub(PFAM_dir+"/run_PFAM.sh", PFAM_dir,1,1)
	Pfam_Job = cluster_utils.launch_job(Pfam_Qsub,cmd_name = 'qsub')
	venn_plot.write('-PFAM %s  ' % (PFAM_dir+"/PFAM.id.txt"))
	Coding_filter.append(Pfam_Job)
if phyloCSF=='yes' or phyloCSF=='y':
	phylocsf_step1=open(phyloCSF_dir+"/run_phyloCSF_step1.sh",'w')
	phylocsf_step2=open(phyloCSF_dir+"/run_phyloCSF_step2.sh",'w')
	phylocsf_step1.write('perl /PUBLIC/source/RNA/lncRNA/Filter/phyloCSF/bin/gtf2bed.pl %s %s\\n' % (step5_dir+"/result/step5.gtf",phyloCSF_dir+"/bedfile"))
	phylocsf_step1.write('python /PUBLIC/software/RNA/galaxy/galaxy-dist2/galaxy-dist/tools/maf/interval_maf_to_merged_fasta.py -d %s -c 1 -s 2 -e 3 -S 6 -m %s -i %s -t user -o %s\\n' % (abbr,maf,phyloCSF_dir+"/bedfile",phyloCSF_dir+"/maf.fa"))
	phylocsf_step1.write('sed -e \\'s/%s./%s|/g\\' %s >%s\\n' % (abbr,abbr,phyloCSF_dir+"/maf.fa",phyloCSF_dir+"/maf.rename.fa"))
	phylocsf_step1.write('awk \\'{print $4\\"\\t\\"$1\\"(\\"$6\\"):\\"$2\\"-\\"$3}\\' %s > %s\\n' % (phyloCSF_dir+"/bedfile",phyloCSF_dir+"/id-co.table"))
	phylocsf_step1.write('mkdir %s\\n' % (phyloCSF_dir+"/MAF_split_fa"))
	phylocsf_step1.write('perl /PUBLIC/source/RNA/lncRNA/Filter/phyloCSF/bin/pre_phylocsf.modi.pl %s %s %s %s \\n' % (phyloCSF_dir+"/id-co.table",abbr,phyloCSF_dir+"/maf.rename.fa",phyloCSF_dir+"/MAF_split_fa"))
	phylocsf_step1.write('perl /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Filter/phyloCSF/bin/generate_phyloCSF.pl -split_fa_dir %s -table %s -outdir %s  -maf_num  %s\\n' % (phyloCSF_dir+"/MAF_split_fa",phyloCSF_dir+"/id-co.table",phyloCSF_dir,abbr))
	phylocsf_step1.write('mkdir %s\\n' % (phyloCSF_dir+"/script/split"))
	phylocsf_step1.write('python /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Assemble/scripture/bin/phylocsf_split.py %s %s\\n' % (phyloCSF_dir+"/script/CSF.sh",phyloCSF_dir+"/script/split"))
	phylocsf_step1.close()
	phylo_step1_Qsub = generate_qsub(phyloCSF_dir+"/run_phyloCSF_step1.sh", phyloCSF_dir,4,2)
	phylo_step1_Job = cluster_utils.launch_job(phylo_step1_Qsub,cmd_name = 'qsub')
	cluster_utils.wait_on_jobs([phylo_step1_Job], cluster_cmd = 'qsub')

	phylocsf_line_num = len(open(phyloCSF_dir+"/script/CSF.sh",'rU').readlines())
	if phylocsf_line_num%1000 == 0:
		n = phylocsf_line_num/1000
	else:
		n = phylocsf_line_num/1000 + 1

	rubyjob_ids=[]
	for each in range(0,n):
		phylocsf_ruby_file=phyloCSF_dir+"/ruby_"+str(each)+".sh"
		phylocsf_ruby=open(phylocsf_ruby_file,'w')
		phylocsf_ruby.write(RUBY_PATH+" "+RUBY_SCRIPT+" "+str(5) +" "+'%s\\n' % (phyloCSF_dir+"/script/split/CSF_"+str(each)+".sh"))
		phylocsf_ruby.close()
		phylo_ruby_Qsub = generate_qsub(phylocsf_ruby_file,phyloCSF_dir+"/script/split", 10, 5)
		rubyjob_ids.append( cluster_utils.launch_job(phylo_ruby_Qsub,cmd_name = 'qsub'))
	cluster_utils.wait_on_jobs(rubyjob_ids, cluster_cmd = 'qsub')
	phylocsf_step2.write('sh %s\\n' % (phyloCSF_dir+"/cat.sh"))
	phylocsf_step2.write('mv %s %s\\n' % (phyloCSF_dir+"/CSF_result.id",phyloCSF_dir+"/phyloCSF.id.txt"))
	phylocsf_step2.write('python /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Filter/phyloCSF/bin/format.change.py %s > %s\\n' % (phyloCSF_dir+"/CSF_result",phyloCSF_dir+"/phyloCSF.result.xls"))
	phylocsf_step2.close()
	phylo_step2_Qsub = generate_qsub(phyloCSF_dir+"/run_phyloCSF_step2.sh",phyloCSF_dir, 4, 2)
	phylo_step2_Job = cluster_utils.launch_job(phylo_step2_Qsub,cmd_name = 'qsub')
	Coding_filter.append(phylo_step2_Job)
	cluster_utils.wait_on_jobs(Coding_filter, cluster_cmd = 'qsub')
	venn_plot.write('-CSF %s ' % (phyloCSF_dir+"/phyloCSF.id.txt"))


venn_plot.write('\\n python /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Filter/venn/get_final_gtf.py %s %s %s >%s\\n' % (venn_dir+"/coding_potential.result.id",step5_dir+"/result/step5.gtf",step5_dir+"/cuffcompare_mRNA.step5_class_code.list.txt",venn_dir+"/lncRNA.gtf"))
venn_plot.close()
venn_plot_Qsub = generate_qsub(venn_dir+"/run_venn.sh",venn_dir, 1, 1)
venn_plot_Job = cluster_utils.launch_job(venn_plot_Qsub,cmd_name = 'qsub')
cluster_utils.wait_on_jobs([venn_plot_Job], cluster_cmd = 'qsub')
os.chdir(step7_dir)
os.system('mkdir exon')
os.system('mkdir length')
os.system('mkdir orf')
os.system('mkdir orf/mRNA')
os.system('mkdir orf/lncRNA')
os.system('mkdir seq')
os.system('mkdir conserve')
exon_dir=os.path.abspath('exon')
length_dir=os.path.abspath('length')
orf_dir=os.path.abspath('orf')
seq_dir=os.path.abspath('seq')
conserve_dir=os.path.abspath('conserve')
lnc_conserve=open(step7_dir+'/conserve.sh','w')
lnc_conserve.write('python /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Filter/bin/picture/get_length_exon.py %s %s %s\\n' % (venn_dir+"/lncRNA.gtf",length_dir+"/lncRNA",exon_dir+"/lncRNA"))
lnc_conserve.write('python /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Filter/bin/picture/get_length_exon.py %s %s %s\\n' % (root_dir+"/5.coding_noncoding_filt/mRNA.gtf",length_dir+"/mRNA",exon_dir+"/mRNA"))
lnc_conserve.write('Rscript /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Filter/bin/picture/length_compare_v2.R %s %s %s\\n' % (length_dir+"/mRNA_length.stat.xls",length_dir+"/lncRNA_length.stat.xls",length_dir))
lnc_conserve.write('Rscript /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Filter/bin/picture/exon_compare_v2.R %s %s %s\\n' % (exon_dir+"/mRNA_exon.stat.xls",exon_dir+"/lncRNA_exon.stat.xls",exon_dir))

lnc_conserve.write('perl /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Enrichment/bin/extractTranscriptsfromFA.pl %s %s %s\\n' % (venn_dir+"/lncRNA.gtf",fa,seq_dir+"/lncRNA.seq"))
lnc_conserve.write('perl /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Enrichment/bin/extractTranscriptsfromFA.pl %s %s %s\\n' % (root_dir+"/5.coding_noncoding_filt/mRNA.gtf",fa,seq_dir+"/mRNA.seq"))
lnc_conserve.write('/opt/bio/EMBOSS/bin/getorf -sequence %s -outseq %s\\n' % (seq_dir+"/lncRNA.seq",orf_dir+"/lncRNA/lncRNA.orf.seq"))
lnc_conserve.write('python /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Filter/bin/picture/lnc_orf_stat.py %s >%s\\n' % (orf_dir+"/lncRNA/lncRNA.orf.seq",orf_dir+"/lncRNA/lncRNA.orf.length.xls"))
##f.write('/PUBLIC/software/RNA/estscan-3.0.3/estscan  -M  /PUBLIC/software/RNA/estscan-3.0.3/%s.smat -o %s -t %s %s\\n' % (matrix,orf_dir+"/mRNA/mRNA.estscan_raw.cds.fasta",orf_dir+"/mRNA/mRNA.estscan_raw.pep.fasta",seq_dir+"/mRNA.seq"))
lnc_conserve.write('/PUBLIC/software/RNA/cufflinks-2.1.1/gffread %s -F -g %s -y %s\\n' % (gff3,fa,orf_dir+"/mRNA/mRNA.orf.seq"))
###f.write('perl /PUBLIC/source/RNA/noRef/ORF_predictor/format_estscan.pl %s %s %s %s\\n' % (orf_dir+"/mRNA/mRNA.estscan_raw.cds.fasta",orf_dir+"/mRNA/mRNA.estscan_raw.pep.fasta",orf_dir+"/mRNA/mRNA.estscan.cds.fasta",orf_dir+"/mRNA/mRNA.estscan.pep.fasta"))
#f.write('python /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Filter/bin/picture/mRNA_orf_stat.py %s > %s\\n' % (orf_dir+"/mRNA/mRNA.orf.seq",orf_dir+"/mRNA/mRNA.orf.length.xls"))
lnc_conserve.write('perl /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Filter/bin/picture/orf_seq.pl %s  %s\\n' % (orf_dir+"/mRNA/mRNA.orf.seq",orf_dir+"/mRNA/mRNA.orf.length.xls"))
lnc_conserve.write('Rscript /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Filter/bin/picture/orf_compare_v2.R %s %s %s\\n' % (orf_dir+"/mRNA/mRNA.orf.length.xls",orf_dir+"/lncRNA/lncRNA.orf.length.xls",orf_dir))
lnc_conserve.write('perl /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Conserve/bin/gtf2bed %s > %s\\n' % (venn_dir+"/lncRNA.gtf",conserve_dir+"/lncRNA.bed"))
lnc_conserve.write('perl /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Conserve/bin/gtf2bed %s > %s\\n' % (root_dir+"/5.coding_noncoding_filt/mRNA.gtf",conserve_dir+"/mRNA.bed"))
lnc_conserve.write('sh /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Conserve/bin/sort-maf.sh %s >%s\\n' % (maf,conserve_dir+"/sorted.maf"))
#lnc_conserve.write("awk '{if(/chr/){print $1}}' %s | sort -u > %s" % (conserve_dir+"/mRNA.bed",conserve_dir+"/chr_list"))
#lnc_conserve.write("perl /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Conserve/bin/maf_getchr.pl oviAri3.chr sorted.maf >sorted.chr.maf" % (conserve_dir+"/chr_list",conserve_dir+"/sorted.maf",conserve_dir+"/sorted.chr.maf"))
#lnc_conserve.write("awk '{if(/^chr/){print}}' %s > %s" % (conserve_dir+"/mRNA.bed",conserve_dir+"/mRNA.chr.bed"))
#lnc_conserve.write("awk '{if(/^chr/){print}}' %s > %s" % (conserve_dir+"/lncRNA.bed",conserve_dir+"/lncRNA.chr.bed"))
lnc_conserve.write('/PUBLIC/software/RNA/phast/phast-1.3/bin/phyloP --method SCORE -f %s /PUBLIC/source/RNA/lncRNA/Filter/phyloCSF/mlin-PhyloCSF-983a652/PhyloCSF_Parameters/%s.mod %s >%s\\n' % (conserve_dir+"/mRNA.bed",abbr,conserve_dir+"/sorted.maf",conserve_dir+"/mRNA.score"))
lnc_conserve.write('/PUBLIC/software/RNA/phast/phast-1.3/bin/phyloP --method SCORE -f %s /PUBLIC/source/RNA/lncRNA/Filter/phyloCSF/mlin-PhyloCSF-983a652/PhyloCSF_Parameters/%s.mod %s >%s\\n' % (conserve_dir+"/lncRNA.bed",abbr,conserve_dir+"/sorted.maf",conserve_dir+"/lncRNA.score"))
lnc_conserve.write('perl /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Conserve/bin/conservationPlot.pl %s %s %s\\n' % (conserve_dir+"/mRNA.score",conserve_dir+"/lncRNA.score",conserve_dir+"/mRNAvslncRNA"))
lnc_conserve.close()
lnc_conserve_Qsub = generate_qsub(step7_dir+'/conserve.sh',step7_dir, 5, 2)
lnc_conserve_Job = cluster_utils.launch_job(lnc_conserve_Qsub,cmd_name = 'qsub')
os.chdir(root_dir+'/5.coding_noncoding_filt')

cluster_utils.wait_on_jobs(mRNA_cufflinks_job, cluster_cmd = 'qsub')

os.chdir(root_dir)
os.chdir('4.mRNA_Cufflinks')
f=open(root_dir+'/4.mRNA_Cufflinks/AS_TSS.sh','w')
f.write('perl /PUBLIC/source/RNA/lncRNA/lncRNA_v1/mRNA/CAN/runCAN_TSS_AS_v2.pl -R %s -G %s -lib %s -i %s -sam %s -o %s  -group %s -groupname %s \\n' % (fa,root_dir+'/4.mRNA_Cufflinks/mRNA.gtf',lib,bam,sam,root_dir+'/4.mRNA_Cufflinks/Cufflinks',group,groupname))
f.write('sh %s/runCAN.sh\\n' % (root_dir+'/4.mRNA_Cufflinks'))
f.close()
assert not os.system('qsub -V -l vf=5g,p=2 -cwd %s/AS_TSS.sh' % (root_dir+'/4.mRNA_Cufflinks'))
os.chdir(root_dir)

'''	
open('lnc_step3_mRNACuff_ncfilter.py','w').write(code)
 
##########################################################################################################
############## DEXseq & cuffdiff
code='''
import os
import os.path
import re
import sys
root_dir='%s'
sample='%s'
bam='%s'
fa='%s'
compare='%s'
group='%s'
groupname='%s'
venn='%s'
lib='%s'
gtf='%s'
sam='%s'
ss='%s'
allsamples_gtf='%s'
allsamples_tmap='%s'
''' % (root_dir,sample,root_dir+"/1.QC/bam",fa,compare,group,groupname,venn,lib,gtf,sam,ss,root_dir+'/4.mRNA_Cufflinks/Cufflinks/merged_gtf_tmap/allsamples.gtf',root_dir+'/4.mRNA_Cufflinks/Cufflinks/merged_gtf_tmap/allsamples.allsamples.gtf.tmap')

flag_DEXseq=True

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

#for DEXSeq
assert not os.system('mkdir 6.mRNA_DEXseq')
'''
if flag_DEXseq:
	if DEXseq_compare:
		code+='''
DEXseq_compare='%s'
''' % (DEXseq_compare)
		code+='''
os.chdir('6.mRNA_DEXseq')
f=open(root_dir+'/6.mRNA_DEXseq/generate_DEXseq.sh','w')

f.write('/PUBLIC/software/public/System/Python-2.7.6/bin/python /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Diff_analysis/bin/runDEXseq_rhd.py -G %s -i %s -g %s -n %s -c %s -o %s -p yes -s %s -a 10 -m %s \\n' % (allsamples_gtf,sam,group,groupname,DEXseq_compare,root_dir+'/6.mRNA_DEXseq/DEXseq',ss,allsamples_tmap))
f.close()
assert not os.system('sh %s' % (root_dir+'/6.mRNA_DEXseq/generate_DEXseq.sh >log 2>error &'))
os.chdir(root_dir)
'''
code+='''

#for SNP
assert not os.system('mkdir 7.mRNA_SNP')
os.chdir('7.mRNA_SNP')
f=open(root_dir+'/7.mRNA_SNP/generate_SNP.sh','w')
f.write('perl /PUBLIC/source/RNA/RefRNA/TransRef/GATKsnp/runGATK2.pl -R %s -t bam -i %s -o %s -b %s -n %s\\n' % (fa,bam,root_dir+'/7.mRNA_SNP/SNP',sample,sample))
f.write('mv %s/workflow.sh %s/workflow.sh\\n' % (root_dir+'/7.mRNA_SNP/SNP',root_dir+'/7.mRNA_SNP'))
f.close()
assert not os.system('sh %s' % (root_dir+'/7.mRNA_SNP/generate_SNP.sh'))
assert not os.system('qsub -V -l vf=10g,p=6 -cwd %s/workflow.sh' % (root_dir+'/7.mRNA_SNP'))
os.chdir(root_dir)



###Cuffdiff
##

os.chdir(root_dir)
os.system('mkdir 8.Diff_analysis')
Diff_analysis_dir=os.path.abspath('8.Diff_analysis')
os.chdir('8.Diff_analysis')
os.system('mkdir Diff')
os.system('mkdir Diff/Diff_Gene_Seq')
os.system('mkdir Diff/expression_compare')
Diff_dir=os.path.abspath('Diff')
groupnames=groupname.strip().split(",")
compare_list=[]
for eachline in compare.strip().split(","):
	each=eachline.strip().split(":")
	compare_list.append(groupnames[int(each[0])-1]+"_vs_"+groupnames[int(each[1])-1])

f=open(Diff_analysis_dir+"/run_Diff_analysis.sh",'w')
f.write('cat %s %s >%s\\n' % (root_dir+"/5.coding_noncoding_filt/mRNA.gtf",root_dir+"/5.coding_noncoding_filt/step6/venn/lncRNA.gtf",Diff_dir+"/mRNA_lncRNA.gtf"))
f.write('perl /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Diff_analysis/bin/getid_transcript.pl %s %s\\n' % (Diff_dir+"/mRNA_lncRNA.gtf",Diff_dir+"/geneid_transcriptid.list"))
'''

if flag_repeat == True:
    code+='''

f.write('python /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Diff_analysis/bin/runCuffdiff_v2.py --gtf %s --is-cufflinks --out-dir %s --bam-dir %s --labels %s --groups %s --groupnames %s --compare %s --library-type %s --thread 2\\n' % (Diff_dir+"/mRNA_lncRNA.gtf",Diff_dir,bam,sample,group,groupname,compare,lib))
for each in compare_list:
	f.write('python /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Diff_analysis/bin/get_mRNA_lncRNA_gene.py %s %s %s %s %s\\n' % (Diff_dir+"/"+each+"/"+each+".Diff_isoform.genelist.xls",root_dir+"/5.coding_noncoding_filt/mRNA.gtf",root_dir+"/5.coding_noncoding_filt/step6/venn/lncRNA.gtf",Diff_dir+"/"+each+"/"+each+".mRNA.genelist.xls",Diff_dir+"/"+each+"/"+each+".lncRNA.genelist.xls"))
	f.write('perl /PUBLIC/source/RNA/RefRNA/DGE/Enrichment/getDiffGene.pl %s %s %s\\n' % (Diff_analysis_dir+"/Diff/Diff_Gene_Seq/gene.seq",Diff_dir+"/"+each+"/"+each+".mRNA.genelist.xls",Diff_analysis_dir+"/Diff/Diff_Gene_Seq/"+each+".mRNA.diffgene.seq"))
	f.write('Rscript /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Diff_analysis/bin/Volcanoplot_v2.R %s %s %s %s\\n' %  (Diff_dir+"/"+each+"/"+each+".isoforms.diff.xls",Diff_dir+"/"+each,each,0.05))

'''

if flag_repeat == False:
    code+='''
f.write('python /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Diff_analysis/bin/runCuffdiff_v2_edgeR.py --gtf %s --is-cufflinks --out-dir %s --bam-dir %s --labels %s --groups %s --groupnames %s --compare %s --library-type %s --thread 2\\n' % (Diff_dir+"/mRNA_lncRNA.gtf",Diff_dir,bam,sample,group,groupname,compare,lib))
f.write('cd %s \\n' %(Diff_dir))
f.write("awk 'NR>1 {print $1}' genes_rowmeans.FPKM.xls | sort > geneIDList\\n")
f.write('perl /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Diff_analysis/bin/extractInfo_v2.pl -i %s -g %s -o %s -l %s\\n' % (Diff_dir+"/mRNA_lncRNA.gtf",Diff_dir+"/geneIDList",Diff_dir+"/geneInfo",Diff_dir+"/genelength"))
f.write('perl /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Diff_analysis/bin/runDE_analysis_edgeR_v2.pl -i %s -r %s -c %s -group %s -groupname %s -g %s -venn_cluster %s -o %s\\n' % (Diff_dir+"/genes.readcount.xls",Diff_dir+"/genes_rowmeans.FPKM.xls",Diff_dir+"/geneInfo",group,groupname,compare,venn,Diff_dir))
f.write('sh Diff_Analysis_edgeR.sh\\n')
f.write('cat *_union_cluster | sort -u > union_for_cluster_nohead\\n')
f.write('head -1 genes_rowmeans.FPKM.xls> head_union_for_cluster\\n')
f.write('cat head_union_for_cluster union_for_cluster_nohead >union_for_cluster\\n')
f.write('rm head_union_for_cluster union_for_cluster_nohead\\n')
f.write('Rscript /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Diff_analysis/bin/cluster.R\\n')
f.write('perl /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Diff_analysis/bin/plot.Venn_Cluster_edgeR.pl -in %s -r %s -fpkm %s -compare %s -venn_cluster %s  -groupname %s -output %s\\n' % (Diff_dir,Diff_dir+"/isoforms_rowmeans.FPKM.xls",Diff_dir+"/isoforms.FPKM.xls",Diff_dir+"/compare.txt",venn,groupname,Diff_dir))

for each in compare_list:
    f.write('python /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Diff_analysis/bin/get_mRNA_lncRNA_gene.py %s %s %s %s %s\\n' % (Diff_dir+"/"+each+"/"+each+".DEGlist.txt",root_dir+"/5.coding_noncoding_filt/mRNA.gtf",root_dir+"/5.coding_noncoding_filt/step6/venn/lncRNA.gtf",Diff_dir+"/"+each+"/"+each+".mRNA.genelist.xls",Diff_dir+"/"+each+"/"+each+".lncRNA.genelist.xls"))
    f.write('perl /PUBLIC/source/RNA/RefRNA/DGE/Enrichment/getDiffGene.pl %s %s %s\\n' % (Diff_analysis_dir+"/Diff/Diff_Gene_Seq/gene.seq",Diff_dir+"/"+each+"/"+each+".mRNA.genelist.xls",Diff_analysis_dir+"/Diff/Diff_Gene_Seq/"+each+".mRNA.diffgene.seq"))

'''
code+='''
f.write('perl /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Diff_analysis/bin/plot.Venn_Cluster.pl -in %s -r %s -fpkm %s -compare %s -venn_cluster %s  -groupname %s -output %s\\n' % (Diff_dir,Diff_dir+"/isoforms_rowmeans.FPKM.xls",Diff_dir+"/isoforms.FPKM.xls",Diff_dir+"/compare.txt",venn,groupname,Diff_dir))
f.write('perl /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Diff_analysis/bin/den_boxplot_v3.R.pl -r %s -rpkm %s -group %s -output %s\\n' % (Diff_analysis_dir+"/Diff/isoforms_rowmeans.FPKM.xls",Diff_analysis_dir+"/Diff/isoforms.FPKM.xls",Diff_analysis_dir+"/Diff/group.txt",Diff_analysis_dir+"/Diff"))
f.write('perl /PUBLIC/source/RNA/RefRNA/DGE/DiffExpression/DendorPCAv2.R.pl %s %s %s\\n' % (Diff_analysis_dir+"/Diff/isoforms.FPKM.xls",Diff_analysis_dir+"/Diff/group.txt",Diff_analysis_dir+"/Diff"))
f.write('perl /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Enrichment/bin/extractcDNAfromFA.pl %s %s %s\\n' % (Diff_dir+"/mRNA_lncRNA.gtf",fa,Diff_analysis_dir+"/Diff/Diff_Gene_Seq/gene.seq"))
f.write('python /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Diff_analysis/bin/get_mRNA_lncRNA_transcript.py %s %s %s %s %s\\n' % (Diff_dir+"/isoforms.FPKM.xls",root_dir+"/5.coding_noncoding_filt/mRNA.gtf",root_dir+"/5.coding_noncoding_filt/step6/venn/lncRNA.gtf",Diff_dir+"/expression_compare/mRNA_FPKM.xls",Diff_dir+"/expression_compare/lncRNA_FPKM.xls"))
f.write('Rscript /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Diff_analysis/bin/expression_compare_v2.R %s %s %s\\n' % (Diff_dir+"/expression_compare/mRNA_FPKM.xls",Diff_dir+"/expression_compare/lncRNA_FPKM.xls",Diff_dir+"/expression_compare"))
f.close()
'''


code+='''
####target
os.chdir(root_dir)
os.system('mkdir 9.Target_predict')
os.system('mkdir 9.Target_predict/cis')
target_dir=os.path.abspath('9.Target_predict')
os.chdir('9.Target_predict')
cis_dir=os.path.abspath('cis')
f=open(cis_dir+"/run_cis.sh",'w')
f.write('python /PUBLIC/source/RNA/lncRNA/lncRNA_v1/cis_trans/cis/bin/gtf_format_change.py %s %s\\n' % (root_dir+"/8.Diff_analysis/Diff/mRNA_lncRNA.gtf",cis_dir+"/mRNA_lncRNA_change.gtf"))
f.write('perl /PUBLIC/source/RNA/lncRNA/lncRNA_v1/cis_trans/cis/bin/exonGTF2geneGTF.pl %s %s\\n' % (cis_dir+"/mRNA_lncRNA_change.gtf",cis_dir+"/mRNA_lncRNA_gene.gtf"))
f.write('python /PUBLIC/source/RNA/lncRNA/lncRNA_v1/cis_trans/cis/bin/find_geneid.py %s > %s\\n' % (root_dir+"/5.coding_noncoding_filt/step6/venn/lncRNA.gtf",cis_dir+"/lncRNA_geneid.txt"))
f.write('perl /PUBLIC/source/RNA/lncRNA/lncRNA_v1/cis_trans/cis/bin/find_cis.pl %s %s 10000 %s/10k.cis %s/10k.cis.o2\\n' % (cis_dir+"/lncRNA_geneid.txt",cis_dir+"/mRNA_lncRNA_gene.gtf",cis_dir,cis_dir))
f.write('perl /PUBLIC/source/RNA/lncRNA/lncRNA_v1/cis_trans/cis/bin/find_cis.pl %s %s 100000 %s/100k.cis %s/100k.cis.o2\\n' % (cis_dir+"/lncRNA_geneid.txt",cis_dir+"/mRNA_lncRNA_gene.gtf",cis_dir,cis_dir))
f.close()

os.chdir(target_dir)
samples=sample.split(',')



if len(samples)>4:
	os.system('mkdir trans')
	trans_dir=os.path.abspath('trans')
	os.system('mkdir trans/trans_out')
	f=open(trans_dir+"/run_trans.sh",'w')
	f.write('python  /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Diff_analysis/bin/get_mRNA_lncRNA_gene_trans.py %s %s %s %s %s\\n' % (root_dir+"/8.Diff_analysis/Diff/genes.FPKM.xls",root_dir+"/5.coding_noncoding_filt/mRNA.gtf",root_dir+"/5.coding_noncoding_filt/step6/venn/lncRNA.gtf",trans_dir+"/mRNA_gene_FPKM.xls",trans_dir+"/lncRNA_gene_FPKM.xls"))
	f.write('python /PUBLIC/source/RNA/lncRNA/lncRNA_v1/cis_trans/trans/bin/lncRNA_trans_predict.py --lnc  %s  --mrna %s --outdir %s --corr 0.95\\n' % (trans_dir+"/lncRNA_gene_FPKM.xls",trans_dir+"/mRNA_gene_FPKM.xls",trans_dir+"/trans_out"))
	f.close()
os.chdir(Diff_analysis_dir)
f=open(Diff_analysis_dir+"/run_Cuffdiff_target.sh",'w')
f.write('sh %s\\n' % (Diff_analysis_dir+"/run_Diff_analysis.sh"))
f.write('sh %s\\n' % (cis_dir+"/run_cis.sh"))
if len(samples)>4:
	f.write('sh %s\\n' % (trans_dir+"/run_trans.sh"))
f.close()
os.chdir(Diff_analysis_dir)
os.system('qsub -cwd -l vf=10G -V -l p=2 %s' % (Diff_analysis_dir+"/run_Cuffdiff_target.sh"))

'''
open('lnc_step4_DEU_Diff_Target.py','w').write(code)


##############################################################################################
######################Enrichment
code='''
import os
import os.path
root_dir='%s'
sample='%s'
fa='%s'
compare='%s'
group='%s'
groupname='%s'
species='%s'
ppi_number='%s'
goann='%s'
gtf='%s'
''' % (root_dir,sample,fa,compare,group,groupname,species,ppi_number,goann,gtf)
code+='''
#groupnames=groupname.strip().split(",")
groupnames=groupname.split(',')
compares=compare.split(',')


compare_list=[]
for eachline in compare.strip().split(","):
        each=eachline.strip().split(":")
        compare_list.append(groupnames[int(each[0])-1]+"_vs_"+groupnames[int(each[1])-1])
compare_list_join=",".join(compare_list)


os.chdir(root_dir)
os.system('mkdir 10.Enrichment')
os.system('mkdir 10.Enrichment/lncRNA')
os.system('mkdir 10.Enrichment/mRNA')
os.system('mkdir 10.Enrichment/tmp')
Enrichment_dir=os.path.abspath("10.Enrichment")


os.system('sed \\'1d\\' %s >%s\\n' % (root_dir+"/9.Target_predict/cis/100k.cis",Enrichment_dir+"/tmp/cis_lncRNA_targets.pairs.xls"))

f=open(Enrichment_dir+"/mRNA/generate_mRNA_enrichment.sh",'w')
f.write('perl /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Enrichment/bin/enrich_ref_mRNA.pl -i %s -c %s -species %s -genome %s -goann %s -gtf %s \\n' % (root_dir+"/8.Diff_analysis/Diff",compare_list_join,species,fa,goann,gtf))
f.close()
os.chdir(Enrichment_dir+"/mRNA")
assert not os.system('sh %s' % (Enrichment_dir+"/mRNA/generate_mRNA_enrichment.sh"))
assert not os.system('nohup sh %s 2>%s.e &' % (Enrichment_dir+"/mRNA/mRNA_enrich.sh",Enrichment_dir+"/mRNA/mRNA_enrich.sh"))

os.chdir(Enrichment_dir+"/lncRNA")
f=open(Enrichment_dir+"/lncRNA/generate_lncRNA_cis_enrichment.sh","w")
f.write('perl /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Enrichment/bin/enrich_ref_lncRNA.pl -i %s -c %s -t %s -species %s -genome %s -goann %s -gtf %s -o cis\\n' % (root_dir+"/8.Diff_analysis/Diff",compare_list_join,Enrichment_dir+"/tmp/cis_lncRNA_targets.pairs.xls",species,fa,goann,gtf))
f.close()
assert not os.system('sh %s' % (Enrichment_dir+"/lncRNA/generate_lncRNA_cis_enrichment.sh"))
assert not os.system('nohup sh %s 2>%s.e &' % (Enrichment_dir+"/lncRNA/cis_enrich.sh",Enrichment_dir+"/lncRNA/cis_enrich.sh"))

samples=sample.split(',')
if len(samples)>4:
	os.system('cut -f1,2 %s |sed \\'1d\\'  >%s\\n' % (root_dir+"/9.Target_predict/trans/trans_out/lncRNA_mRNA.trans.xls",Enrichment_dir+"/tmp/trans_lncRNA_targets.pairs.xls"))
	f=open(Enrichment_dir+"/lncRNA/generate_lncRNA_trans_enrichment.sh","w")
	f.write('perl /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Enrichment/bin/enrich_ref_lncRNA.pl -i %s -c %s -t %s -species %s -genome %s -goann %s -gtf %s -o trans\\n' % (root_dir+"/8.Diff_analysis/Diff",compare_list_join,Enrichment_dir+"/tmp/trans_lncRNA_targets.pairs.xls",species,fa,goann,gtf))
	f.close()
	assert not os.system('sh %s' % (Enrichment_dir+"/lncRNA/generate_lncRNA_trans_enrichment.sh"))
	assert not os.system('nohup sh %s 2>%s.e &' % (Enrichment_dir+"/lncRNA/trans_enrich.sh",Enrichment_dir+"/lncRNA/trans_enrich.sh"))

os.chdir(root_dir)
os.system('mkdir 11.PPI')
PPI_dir=os.path.abspath('11.PPI')
os.chdir('11.PPI')
os.system('mkdir mRNA')
os.system('mkdir lncRNA')
os.system('mkdir lncRNA/cis')

mRNA_dir=os.path.abspath('mRNA')
cis_dir=os.path.abspath('lncRNA/cis')
os.chdir('mRNA')
for each in compare_list:
	os.system('mkdir %s' % (each))
	os.system('perl /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Diff_analysis/bin/getDiffGene.pl %s %s %s' % (root_dir+"/8.Diff_analysis/Diff/Diff_Gene_Seq/gene.seq",root_dir+"/8.Diff_analysis/Diff/"+each+"/"+each+".mRNA.genelist.xls",mRNA_dir+"/"+each+"/"+each+".Diffgene.seq"))
	f=open(mRNA_dir+"/"+each+"/run_PPI_"+each+".sh",'w')
	f.write('python /PUBLIC/source/RNA/RefRNA/DGE/Enrichment/BLASTX_TO_PPI_v4.py --species %s  --fa %s  --outdir %s  --name %s' % (ppi_number,mRNA_dir+"/"+each+"/"+each+".Diffgene.seq",mRNA_dir+"/"+each,each))
	f.close()
	os.chdir(mRNA_dir+"/"+each)
	assert not os.system('qsub -cwd -l vf=2G -V %s' % (mRNA_dir+"/"+each+"/run_PPI_"+each+".sh"))
	os.chdir(mRNA_dir)

os.chdir(cis_dir)
for each in compare_list:
	os.system('mkdir %s' % (each))
	os.system('python /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Enrichment/bin/get_lncRNA_target_seq.py %s %s %s >%s' % (root_dir+"/8.Diff_analysis/Diff/Diff_Gene_Seq/gene.seq",root_dir+"/8.Diff_analysis/Diff/"+each+"/"+each+".lncRNA.genelist.xls",root_dir+"/9.Target_predict/cis/100k.cis",cis_dir+"/"+each+"/"+each+".Diffgene.seq"))
	f=open(cis_dir+"/"+each+"/run_PPI_"+each+".cis.sh","w")
	f.write('python /PUBLIC/source/RNA/RefRNA/DGE/Enrichment/BLASTX_TO_PPI_v4.py --species %s  --fa %s  --outdir %s  --name %s' % (ppi_number,cis_dir+"/"+each+"/"+each+".Diffgene.seq",cis_dir+"/"+each,each))
	f.close()
	os.chdir(cis_dir+"/"+each)
	assert not os.system('qsub -cwd -l vf=2G -V %s' % (cis_dir+"/"+each+"/run_PPI_"+each+".cis.sh"))
        os.chdir(cis_dir)
os.chdir(PPI_dir)
if len(samples)>4:
	os.system('mkdir lncRNA/trans')
	trans_dir=os.path.abspath('lncRNA/trans')
	os.chdir(trans_dir)
	for each in compare_list:
		os.system('mkdir %s' % (each))
		os.system('python /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Enrichment/bin/get_lncRNA_target_seq.py %s %s %s >%s' % (root_dir+"/8.Diff_analysis/Diff/Diff_Gene_Seq/gene.seq",root_dir+"/8.Diff_analysis/Diff/"+each+"/"+each+".lncRNA.genelist.xls",root_dir+"/9.Target_predict/trans/trans_out/lncRNA_mRNA.trans.xls",trans_dir+"/"+each+"/"+each+".Diffgene.seq"))
		f=open(trans_dir+"/"+each+"/run_PPI_"+each+".trans.sh","w")
		f.write('python /PUBLIC/source/RNA/RefRNA/DGE/Enrichment/BLASTX_TO_PPI_v4.py --species %s  --fa %s  --outdir %s  --name %s' % (ppi_number,trans_dir+"/"+each+"/"+each+".Diffgene.seq",trans_dir+"/"+each,each))
		f.close()
		os.chdir(trans_dir+"/"+each)
		assert not os.system('qsub -cwd -l vf=2G -V %s' % (trans_dir+"/"+each+"/run_PPI_"+each+".trans.sh"))
		os.chdir(trans_dir)


'''
open('lnc_step5_Enrichment_PPI.py','w').write(code)

###################################################################################
##################report & result



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
settings.configure(DEBUG=True, TEMPLATE_DEBUG=True,TEMPLATE_DIRS=('/PUBLIC/source/RNA/lncRNA/lncRNA_v1/report/new',))
root_dir='%s'
sample='%s'
samples=sample.split(',')
project='%s'
DEXseq='%s'
ss='%s'
phyloCSF='%s'
venn='%s'
CNCI='%s'
CPC='%s'
PFAM='%s'
venn_plot='%s'
group='%s'
groupname='%s'
compare='%s'
venn_cluster_name='%s'
flag_repeat = %s
DEXseq_compare ='%s'
''' % (root_dir,sample,project,DEXseq,ss,phyloCSF,venn,CNCI,CPC,PFAM,venn_plot,group,groupname,compare,venn_cluster_name,flag_repeat,DEXseq_compare)

code+='''
os.chdir(root_dir)
result_dir=root_dir+"/"+project+'_lncRNA_result/'+project+'_results'

groupnames=groupname.split(',')
DEXseq_compare_names=[]
if DEXseq_compare !='':
        DEXseq_compares=DEXseq_compare.split(',')
        for each in DEXseq_compares:
                tmp=each.split(':')
                DEXseq_compare_names.append(groupnames[int(tmp[0])-1]+'vs'+groupnames[int(tmp[1])-1])


groupnames=groupname.split(',')
compares=compare.split(',')
compare_names=[]
for each in compares:
        temp=each.split(':')
        compare_names.append(groupnames[int(temp[0])-1]+'_vs_'+groupnames[int(temp[1])-1])
venn_cluster_vs_names=[]
flag_venn=False
for each in venn_cluster_name.split(','):
        if each.count(':') >= 2 and each.count(':') <= 5:
                flag_venn=True
                venn_cluster_vs_names.append(each.replace(':','_vs_'))

'''

if flag_results == True:
	code+='''
dir=root_dir+'/1.QC'
cd_md5=open(root_dir+'/1.QC/cd_md5.sh','w')
rd_md5=open(root_dir+'/1.QC/rd_md5.sh','w')
os.chdir(dir)
cd_md5.write('cd %s/cleandata\\n' % (dir))
cd_md5.write('echo -e \"md5\\tfile\">cd_md5.txt\\n')
rd_md5.write('echo -e \"md5\\tfile\">rd_md5.txt\\n')
for each in samples:
        cd_md5.write('unlink %s_1.clean.fq\\n' %(each))
        cd_md5.write('unlink %s_2.clean.fq\\n' %(each))
        cd_temp1=root_dir+'/1.QC/'+each+'/clean_data/'+each+'_1.clean.fq'
        cd_temp2=root_dir+'/1.QC/'+each+'/clean_data/'+each+'_2.clean.fq'
        cd_md5.write('gzip %s\\n' % (cd_temp1))
        cd_md5.write('gzip %s\\n' % (cd_temp2))
        cd_md5.write('ln -s %s.gz %s\\n' % (cd_temp1,dir+'/cleandata/'))
        cd_md5.write('ln -s %s.gz %s\\n' % (cd_temp2,dir+'/cleandata/'))
        cd_md5.write('md5sum %s_1.clean.fq.gz>>cd_md5.txt\\n' % (each))
        cd_md5.write('md5sum %s_2.clean.fq.gz>>cd_md5.txt\\n' % (each))
        rd_md5.write('cd %s\\n' % (root_dir+'/1.QC/'+each))
        rd_md5.write('md5sum %s>>%srd_md5.txt\\n' % (each+'_1.fq.gz',root_dir+'/1.QC/'))
        rd_md5.write('md5sum %s>>%srd_md5.txt\\n' % (each+'_2.fq.gz',root_dir+'/1.QC/'))
cd_md5.close()
rd_md5.close()
os.chdir(root_dir)

os.chdir(root_dir+'/1.QC')
rdmd5=glob.glob('rd_md5.txt')
if not rdmd5:
        assert not os.system('qsub -V -l vf=1g,p=1 -cwd rd_md5.sh')
os.chdir(root_dir+'/1.QC/cleandata')
cdmd5=glob.glob('cd_md5.txt')
if not cdmd5:
        os.chdir(root_dir+'/1.QC')
        assert not os.system('qsub -V -l vf=1g,p=1 -cwd cd_md5.sh')
os.chdir(root_dir)
	

for each in compare_names:
    os.chdir(root_dir+"/10.Enrichment/mRNA/"+each+"/Pathway")
    os.system('python /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Enrichment/bin/pathway_annotation_flow_parallel_simple_tolerant.pyc --table %s' % (root_dir+"/10.Enrichment/mRNA/"+each+"/Pathway/add."+each+".indentify.xls"))
    os.chdir(root_dir+"/10.Enrichment/lncRNA/cis/"+each+"/Pathway")
    os.system('python /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Enrichment/bin/pathway_annotation_flow_parallel_simple_tolerant.pyc --table %s' % (root_dir+"/10.Enrichment/lncRNA/cis/"+each+"/Pathway/add."+each+".indentify.xls"))
    samples=sample.split(',')
    if len(samples)>4:
        os.chdir(root_dir+"/10.Enrichment/lncRNA/trans/"+each+"/Pathway")
        os.system('python /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Enrichment/bin/pathway_annotation_flow_parallel_simple_tolerant.pyc --table %s' % (root_dir+"/10.Enrichment/lncRNA/trans/"+each+"/Pathway/add."+each+".indentify.xls"))

os.chdir(root_dir)
assert not os.system('mkdir '+project+'_lncRNA_result')
assert not os.system('mkdir '+project+'_lncRNA_result/'+project+'_results')
os.chdir(result_dir)
result_order=0
result_order+=1
assert not os.system('mkdir '+str(result_order)+'.QC')
qc_order=str(result_order)
os.chdir(''+str(result_order)+'.QC')
assert not os.system('mkdir '+str(result_order)+'.1.ErrorRate')
os.chdir(''+str(result_order)+'.1.ErrorRate')
for eachsample in samples:
        assert not os.system('cp %s %s' % (root_dir+'/1.QC/'+eachsample+'/clean_data/'+eachsample+'.Error.png',eachsample+'.error_rate_distribution.png'))
        assert not os.system('cp %s %s' % (root_dir+'/1.QC/'+eachsample+'/clean_data/'+eachsample+'.Error.pdf',eachsample+'.error_rate_distribution.pdf'))
assert not os.system('cp /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Readme/Error.README.txt %s' % (result_dir+"/"+str(result_order)+'.QC/'+str(result_order)+'.1.ErrorRate/'+str(result_order)+'.1ErrorRate.README.txt'))
os.chdir('..')

assert not os.system('mkdir '+str(result_order)+'.2.GC')
os.chdir(''+str(result_order)+'.2.GC')
for eachsample in samples:
        assert not os.system('cp %s %s' % (root_dir+'/1.QC/'+eachsample+'/clean_data/'+eachsample+'.GC.png',eachsample+'.GC_content_distribution.png'))
        assert not os.system('cp %s %s' % (root_dir+'/1.QC/'+eachsample+'/clean_data/'+eachsample+'.GC.pdf',eachsample+'.GC_content_distribution.pdf'))
assert not os.system('cp /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Readme/GC.README.txt %s' % (result_dir+"/"+str(result_order)+'.QC/'+str(result_order)+'.2.GC/'+str(result_order)+'.2GC.README.txt'))
os.chdir('..')


assert not os.system('mkdir '+str(result_order)+'.3.ReadsClassification')
os.chdir(''+str(result_order)+'.3.ReadsClassification')
for eachsample in samples:
        assert not os.system('cp %s %s' % (root_dir+'/1.QC/'+eachsample+'/files/'+eachsample+'.pie3d.png',eachsample+'.raw_reads_classification.png'))
        assert not os.system('cp %s %s' % (root_dir+'/1.QC/'+eachsample+'/files/'+eachsample+'.pie3d.pdf',eachsample+'.raw_reads_classification.pdf'))
assert not os.system('cp /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Readme/Filter.README.txt %s' % (result_dir+"/"+str(result_order)+'.QC/'+str(result_order)+'.3.ReadsClassification/'+str(result_order)+'.2ReadsClassification.README.txt'))
os.chdir('..')

assert not os.system('mkdir '+str(result_order)+'.4.DataTable')
os.chdir(''+str(result_order)+'.4.DataTable')
assert not os.system("cat %s|cut -f 1-8 > %s " % (root_dir+'/1.QC/*/files/dataTable','datatable.xls.1'))       
assert not os.system('cat %s %s > %s' % ('/PUBLIC/source/RNA/lncRNA/lncRNA_v1/Readme/datatable','datatable.xls.1','datatable.xls'))
assert not os.system('rm %s' % ('datatable.xls.1'))
assert not os.system('cp /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Readme/DataTable.README.txt %s' % (result_dir+"/"+str(result_order)+'.QC/'+str(result_order)+'.4.DataTable/'+str(result_order)+'.3DataTable.README.txt'))
os.chdir(result_dir)

result_order+=1
assert not os.system('mkdir '+str(result_order)+'.Mapping')
mp_order=str(result_order)
os.chdir(''+str(result_order)+'.Mapping')
assert not os.system('mkdir '+str(result_order)+'.1.MapStat')
os.chdir(''+str(result_order)+'.1.MapStat')
for eachsample in samples:
        assert not os.system('cut -f 2 %s > %s' % (root_dir+'/1.QC/'+eachsample+'/files/'+eachsample+'.stat',eachsample+'.stat.2'))
assert not os.system('cut -f 1 %s > %s' % (root_dir+'/1.QC/'+samples[0]+'/files/'+samples[0]+'.stat','stat.1'))
assert not os.system('paste %s %s > %s' % ('stat.1','*.stat.2','MapStat.xls'))
assert not os.system("sed -i -e '/^Reads mapped in proper pairs/d' %s" % ('MapStat.xls'))
assert not os.system("sed -i -e '/^Proper-paired reads map to different chrom/d' %s" % ('MapStat.xls'))
assert not os.system('rm stat.1')
assert not os.system('rm *.stat.2')
assert not os.system('cp /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Readme/MapStat.README.txt %s' % (result_dir+"/"+str(result_order)+'.Mapping/'+str(result_order)+'.1.MapStat/'+str(result_order)+'.1MapStat.README.txt'))
os.chdir('..')

assert not os.system('mkdir '+str(result_order)+'.2.MapReg')
mr_order=str(result_order)
os.chdir(''+str(result_order)+'.2.MapReg')
for eachsample in samples:
        assert not os.system('cp %s %s' % (root_dir+'/1.QC/'+eachsample+'/HTseq/'+eachsample+'.MR_feature.png',eachsample+'.Mapped_Region.png'))
        assert not os.system('cp %s %s' % (root_dir+'/1.QC/'+eachsample+'/files/'+eachsample+'.MR_feature.pdf',eachsample+'.Mapped_Region.pdf'))
	assert not os.system('cp %s %s' % (root_dir+'/1.QC/'+eachsample+'/HTseq/mappedfeature.txt',eachsample+'.Mapped_Region.xls'))
	assert not os.system('python /PUBLIC/source/RNA/lncRNA/lncRNA_v1/QC/bin/feature_stat_combine.py %s %s >all.mappedfeature.xls' % (root_dir,sample))
assert not os.system('cp /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Readme/MapReg.README.txt %s' % (result_dir+"/"+str(result_order)+'.Mapping/'+str(result_order)+'.2.MapReg/'+str(result_order)+'.2MappedRegion.README.txt'))
os.chdir('..')

assert not os.system('mkdir '+str(result_order)+'.3.ChrDen')
cd_order=str(result_order)
os.chdir(''+str(result_order)+'.3.ChrDen')
for eachsample in samples:
        assert not os.system('cp %s .' % (root_dir+'/3.RNA_advanced_QC/Density/'+eachsample+'/'+eachsample+'.density.png'))
        assert not os.system('cp %s .' % (root_dir+'/3.RNA_advanced_QC/Density/'+eachsample+'/'+eachsample+'.density.pdf'))
assert not os.system('cp /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Readme/ChrDen.README.txt %s' % (result_dir+"/"+str(result_order)+'.Mapping/'+str(result_order)+'.3.ChrDen/'+str(result_order)+'.3ChrDen.README.txt'))
os.chdir('..')

assert not os.system('mkdir '+str(result_order)+'.4.IGV')
igv_order=str(result_order)
os.chdir(''+str(result_order)+'.4.IGV')
assert not os.system('cp /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Readme/IGVQuickStart.pdf %s' % (result_dir+"/"+str(result_order)+'.Mapping/'+str(result_order)+'.4.IGV'))
assert not os.system('cp /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Readme/IGV.README.txt %s' % (result_dir+"/"+str(result_order)+'.Mapping/'+str(result_order)+'.4.IGV/'+str(result_order)+'.4IGV.README.txt'))
os.chdir(result_dir)

result_order+=1
assert not os.system('mkdir '+str(result_order)+'.AS')
as_order=str(result_order)
os.chdir(''+str(result_order)+'.AS')
assert not os.system('cp %s .' % (root_dir+'/4.mRNA_Cufflinks/Cufflinks/ASprofile/AS.png'))
assert not os.system('cp %s .' % (root_dir+'/4.mRNA_Cufflinks/Cufflinks/ASprofile/AS.pdf'))
for eachsample in samples:
        assert not os.system('cp %s .' % (root_dir+'/4.mRNA_Cufflinks/Cufflinks/ASprofile/'+eachsample+'/'+eachsample+'.fpkm.xls'))
assert not os.system('cp %s .' % (root_dir+'/4.mRNA_Cufflinks/Cufflinks/ASprofile/ASevent.stat.xls'))
assert not os.system('cp /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Readme/AS.README.txt %s' % (result_dir+"/"+str(result_order)+'.AS/'+str(result_order)+'AS.README.txt'))
os.chdir(result_dir)

result_order+=1
assert not os.system('mkdir '+str(result_order)+'.SNP_InDel')
snp_order=str(result_order)
os.chdir(''+str(result_order)+'.SNP_InDel')
assert not os.system('cp %s .' % (root_dir+'/7.mRNA_SNP/SNP/ResultsQ30/InDels.xls'))
assert not os.system('cp %s .' % (root_dir+'/7.mRNA_SNP/SNP/ResultsQ30/SNPs.xls'))
assert not os.system('cp /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Readme/SNP.README.txt %s' % (result_dir+"/"+str(result_order)+'.SNP_InDel/'+str(result_order)+'SNP_InDel.README.txt'))
os.chdir(result_dir)
'''
if flag_results == True:
	if flag_DEXseq:
		if DEXseq=='yes' or DEXseq=='y':
			code+='''
DEXseq_compare = '%s'
''' % (DEXseq_compare)
		code+='''
groupnames=groupname.split(',')
DEXseq_compare_names=[]
if DEXseq_compare !='':
        DEXseq_compares=DEXseq_compare.split(',')
        for each in DEXseq_compares:
                tmp=each.split(':')
                DEXseq_compare_names.append(groupnames[int(tmp[0])-1]+'vs'+groupnames[int(tmp[1])-1])

result_order+=1
assert not os.system('mkdir '+str(result_order)+'.DEU')
os.chdir(''+str(result_order)+'.DEU')
for each in DEXseq_compare_names:
        assert not os.system('cp -r %s .' % (root_dir+'/6.mRNA_DEXseq/DEXseq/'+each+'/'+each+'.DEU'))
        assert not os.system('cp %s %s' % (root_dir+'/6.mRNA_DEXseq/DEXseq/'+each+'/'+each+'.DEUdiff.xls',each+'.DEU'))
assert not os.system('cp /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Readme/DEU.README.txt %s' % (result_dir+"/"+str(result_order)+'.DEU/'+str(result_order)+'.DEU.README.txt'))
os.chdir(result_dir)
'''
if flag_results == True:
	if ss=='yes' or ss=='reverse':
		code+='''
os.chdir(result_dir)
result_order+=1
assert not os.system('mkdir '+str(result_order)+'.AntiTrans')
os.chdir(''+str(result_order)+'.AntiTrans')
for eachsample in samples:
        assert not os.system('cp %s .' % (root_dir+'/4.mRNA_Cufflinks/Cufflinks/stranded_specific/cisNAT/'+eachsample+'.nat.xls'))
        assert not os.system('cp %s .' % (root_dir+'/4.mRNA_Cufflinks/Cufflinks/stranded_specific/cisNAT/'+eachsample+'.nat.stat.xls'))
assert not os.system('cp /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Readme/AntiTrans.README.txt %s' % (result_dir+"/"+str(result_order)+'.AntiTrans/'+str(result_order)+'.AntiTrans.README.txt'))
os.chdir(result_dir)

result_order+=1
assert not os.system('mkdir '+str(result_order)+'.TSS_TTS')
os.chdir(''+str(result_order)+'.TSS_TTS')
assert not os.system('cp %s %s' % (root_dir+'/4.mRNA_Cufflinks/Cufflinks/stranded_specific/TSS.xls','TSS_TTS.xls'))
assert not os.system('cp /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Readme/TSS_TTS.README.txt %s' % (result_dir+"/"+str(result_order)+'.TSS_TTS/'+str(result_order)+'.TSS_TTS.README.txt'))
os.chdir(result_dir)
'''
if flag_results == True:
	code+='''
result_order+=1
assert not os.system('mkdir '+str(result_order)+'.lncRNA_Assembly')
os.chdir(''+str(result_order)+'.lncRNA_Assembly')
assert not os.system('mkdir cufflinks')
assert not os.system('mkdir scripture')
for eachsample in samples:
	assert not os.system('cp %s %s' % (root_dir+'/2.lnc_assemble/cufflinks/'+eachsample+"/transcripts.gtf",result_dir+"/"+str(result_order)+'.lncRNA_Assembly/cufflinks/'+eachsample+"_cufflinks.gtf"))
	assert not os.system('cp %s %s' % (root_dir+'/2.lnc_assemble/scripture/'+eachsample+"/"+eachsample+".all.bed",result_dir+"/"+str(result_order)+'.lncRNA_Assembly/scripture/'+eachsample+"_scripture.bed"))
assert not os.system('cp /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Readme/lncRNA_Assembly.README.txt %s' % (result_dir+"/"+str(result_order)+'.lncRNA_Assembly/'+str(result_order)+".lncRNA_Assembly.README.txt"))
os.chdir(result_dir)

result_order+=1
assert not os.system('mkdir '+str(result_order)+'.lncRNA_Filter')
os.chdir(''+str(result_order)+'.lncRNA_Filter')
assert not os.system('mkdir '+str(result_order)+'.1Basic_Filter')
assert not os.system('mkdir '+str(result_order)+'.2Coding_Potentiality_Filter')
assert not os.system('mkdir '+str(result_order)+'.3Final_result')
assert not os.system('cp %s %s' % (root_dir+"/5.coding_noncoding_filt/step5/result/step_num_stat.xls",str(result_order)+'.1Basic_Filter'))
assert not os.system('cp %s %s' % (root_dir+"/5.coding_noncoding_filt/step5/result/step_number.pdf",str(result_order)+'.1Basic_Filter'))
assert not os.system('cp %s %s' % (root_dir+"/5.coding_noncoding_filt/step5/result/step_number.png",str(result_order)+'.1Basic_Filter'))
assert not os.system('cp %s %s' % (root_dir+"/5.coding_noncoding_filt/step5/result/class_code.png",str(result_order)+'.1Basic_Filter'))
assert not os.system('cp %s %s' % (root_dir+"/5.coding_noncoding_filt/step5/result/class_code.pdf",str(result_order)+'.1Basic_Filter'))
if CNCI=='yes' or CNCI=='y':
	assert not os.system('mkdir '+str(result_order)+'.2Coding_Potentiality_Filter/CNCI')
	assert not os.system('cp %s %s' % (root_dir+"/5.coding_noncoding_filt/step6/CNCI/CNCI.result.xls",str(result_order)+'.2Coding_Potentiality_Filter/CNCI'))
	assert not os.system('cp %s %s' % (root_dir+"/5.coding_noncoding_filt/step6/CNCI/CNCI.id.txt",str(result_order)+'.2Coding_Potentiality_Filter/CNCI'))
if CPC=='yes' or CPC=='y':
	assert not os.system('mkdir '+str(result_order)+'.2Coding_Potentiality_Filter/CPC')
        assert not os.system('cp %s %s' % (root_dir+"/5.coding_noncoding_filt/step6/CPC/CPC.result.xls",str(result_order)+'.2Coding_Potentiality_Filter/CPC'))
        assert not os.system('cp %s %s' % (root_dir+"/5.coding_noncoding_filt/step6/CPC/CPC.id.txt",str(result_order)+'.2Coding_Potentiality_Filter/CPC'))
if PFAM=='yes' or PFAM=='y':
	assert not os.system('mkdir '+str(result_order)+'.2Coding_Potentiality_Filter/PFAM')
        assert not os.system('cp %s %s' % (root_dir+"/5.coding_noncoding_filt/step6/PFAM/PFAM.result.xls",str(result_order)+'.2Coding_Potentiality_Filter/PFAM'))
        assert not os.system('cp %s %s' % (root_dir+"/5.coding_noncoding_filt/step6/PFAM/PFAM.id.txt",str(result_order)+'.2Coding_Potentiality_Filter/PFAM'))
if phyloCSF=='yes' or phyloCSF=='y':
	assert not os.system('mkdir '+str(result_order)+'.2Coding_Potentiality_Filter/phyloCSF')
        assert not os.system('cp %s %s' % (root_dir+"/5.coding_noncoding_filt/step6/phyloCSF/phyloCSF.result.xls",str(result_order)+'.2Coding_Potentiality_Filter/phyloCSF'))
        assert not os.system('cp %s %s' % (root_dir+"/5.coding_noncoding_filt/step6/phyloCSF/phyloCSF.id.txt",str(result_order)+'.2Coding_Potentiality_Filter/phyloCSF'))
assert not os.system('mkdir '+str(result_order)+'.2Coding_Potentiality_Filter/Venn')
assert not os.system('cp %s %s' % (root_dir+"/5.coding_noncoding_filt/step6/venn/coding_potential.result.id",str(result_order)+'.2Coding_Potentiality_Filter/Venn'))
assert not os.system('cp %s %s' % (root_dir+"/5.coding_noncoding_filt/step6/venn/coding_potential.result.xls",str(result_order)+'.2Coding_Potentiality_Filter/Venn'))
assert not os.system('cp %s %s' % (root_dir+"/5.coding_noncoding_filt/step6/venn/union.venn.png",str(result_order)+'.2Coding_Potentiality_Filter/Venn'))
assert not os.system('cp %s %s' % (root_dir+"/5.coding_noncoding_filt/step6/venn/union.venn.pdf",str(result_order)+'.2Coding_Potentiality_Filter/Venn'))
assert not os.system('cp %s %s' % (root_dir+"/5.coding_noncoding_filt/step6/venn/lncRNA.gtf",str(result_order)+'.3Final_result'))
assert not os.system('cp %s %s' % (root_dir+"/5.coding_noncoding_filt/step7_compare/seq/lncRNA.seq",str(result_order)+'.3Final_result'))
#assert not os.system('cp %s %s' % (root_dir+"/5.coding_noncoding_filt/step7_compare/seq/mRNA.seq",str(result_order)+'.3Final_result'))
assert not os.system('cp /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Readme/lncRNA_Filter.README.txt %s' % (result_dir+"/"+str(result_order)+'.lncRNA_Filter/'+str(result_order)+".lncRNA_Filter.README.txt"))
os.chdir(result_dir)

result_order+=1
assert not os.system('mkdir '+str(result_order)+'.Quatification')
os.chdir(''+str(result_order)+'.Quatification')
assert not os.system('cp %s .' % (root_dir+"/8.Diff_analysis/Diff/expression_compare/lncRNA_FPKM.xls"))
assert not os.system('cp %s .' % (root_dir+"/8.Diff_analysis/Diff/expression_compare/mRNA_FPKM.xls"))
assert not os.system('cp %s %s' % (root_dir+"/8.Diff_analysis/Diff/isoforms.FPKM.xls",result_dir+"/"+str(result_order)+'.Quatification/FPKM.xls'))
assert not os.system('cp /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Readme/Quatification.README.txt %s' % (result_dir+"/"+str(result_order)+'.Quatification/'+str(result_order)+'.Quatification.README.txt'))
os.chdir(result_dir)

result_order+=1
assert not os.system('mkdir '+str(result_order)+'.AdvancedQC')
os.chdir(''+str(result_order)+'.AdvancedQC')
assert not os.system('mkdir '+str(result_order)+'.1.DupCorr')
assert not os.system('mkdir '+str(result_order)+'.3.MeanCov')
assert not os.system('mkdir '+str(result_order)+'.2.ExpLev')
for png in glob.glob(root_dir+'/8.Diff_analysis/Diff/corr_plot/*.png'):
                assert not os.system('cp %s %s' % (png,str(result_order)+'.1.DupCorr'))
for pdf in glob.glob(root_dir+'/8.Diff_analysis/Diff/corr_plot/*.pdf'):
                assert not os.system('cp %s %s' % (pdf,str(result_order)+'.1.DupCorr'))
assert not os.system('cp /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Readme/DupCorr.README.txt %s' % (str(result_order)+'.1.DupCorr/'+str(result_order)+'.1DupCorr.README.txt'))
'''
if flag_results == True:
    if flag_uniform:
        code+='''
for eachsample in samples:
	assert not os.system('perl /PUBLIC/source/RNA/lncRNA/lncRNA_v1/QC/bin/plot_MC_lncRNA.pl %s %s' % (eachsample,root_dir+"/3.RNA_advanced_QC/UniformDistribution/"+eachsample+"/"+eachsample+".report/"))
	assert not os.system('cp %s %s' % (root_dir+"/3.RNA_advanced_QC/UniformDistribution/"+eachsample+"/"+eachsample+".report/"+eachsample+".MC.png",str(result_order)+'.3.MeanCov/'+eachsample+'.Mean_coverage_distribution.png'))
	assert not os.system('cp %s %s' % (root_dir+"/3.RNA_advanced_QC/UniformDistribution/"+eachsample+"/"+eachsample+".report/"+eachsample+".MC.pdf",str(result_order)+'.3.MeanCov/'+eachsample+'.Mean_coverage_distribution.pdf'))
assert not os.system('cp /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Readme/MeanCov.README.txt %s' % (str(result_order)+'.3.MeanCov/'+str(result_order)+'.3.MeanCov.README.txt'))
'''
if flag_results == True:
	code+='''
assert not os.system('cp %s %s' % (root_dir+"/8.Diff_analysis/Diff/boxplot.pdf",str(result_order)+'.2.ExpLev'))
assert not os.system('cp %s %s' % (root_dir+"/8.Diff_analysis/Diff/density.pdf",str(result_order)+'.2.ExpLev'))
assert not os.system('cp %s %s' % (root_dir+"/8.Diff_analysis/Diff/boxplot.png",str(result_order)+'.2.ExpLev'))
assert not os.system('cp %s %s' % (root_dir+"/8.Diff_analysis/Diff/density.png",str(result_order)+'.2.ExpLev'))
assert not os.system('cp /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Readme/ExpLev.README.txt %s' % (str(result_order)+'.2.ExpLev/'+str(result_order)+'.2.ExpLev.README.txt'))
os.chdir(result_dir)

groupnames=groupname.split(',')
compares=compare.split(',')
compare_names=[]
for each in compares:
        temp=each.split(':')
        compare_names.append(groupnames[int(temp[0])-1]+'_vs_'+groupnames[int(temp[1])-1])
venn_cluster_vs_names=[]
flag_venn=False
for each in venn_cluster_name.split(','):
        if each.count(':') >= 2 and each.count(':') <= 5:
                flag_venn=True
                venn_cluster_vs_names.append(each.replace(':','_vs_'))
result_order+=1
assert not os.system('mkdir '+str(result_order)+'.DiffExprAnalysis')
os.chdir(''+str(result_order)+'.DiffExprAnalysis')
assert not os.system('mkdir '+str(result_order)+'.1.DEGsList')
os.chdir(''+str(result_order)+'.1.DEGsList')
'''
if flag_results == True:
	if flag_repeat == True:
		code+='''
for each in compare_names:
        assert not os.system('cp %s %s' % (root_dir+'/8.Diff_analysis/Diff/'+each+"/"+each+'.isoforms.diff.xls',each+'.Differential_analysis_results.xls'))
        assert not os.system('cp %s %s' % (root_dir+'/8.Diff_analysis/Diff/'+each+"/"+each+'.isoforms.filter.xls',each+'.DE.xls'))
        assert not os.system('cp %s %s' % (root_dir+'/8.Diff_analysis/Diff/'+each+"/"+each+'.Diff_isoform.list.xls',each+'.DEG_list.xls'))
        assert not os.system('cp %s %s' % (root_dir+'/8.Diff_analysis/Diff/'+each+"/"+each+'.Diff_isoform.genelist.xls',each+'.DEG_genelist.xls'))
        assert not os.system('cp %s %s' % (root_dir+'/8.Diff_analysis/Diff/'+each+"/"+each+'.lncRNA.genelist.xls',each+'.DEG_lncRNA_genelist.xls'))
        assert not os.system('cp %s %s' % (root_dir+'/8.Diff_analysis/Diff/'+each+"/"+each+'.mRNA.genelist.xls',each+'.DEG_mRNA_genelist.xls'))
assert not os.system('cp /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Readme/DiffList.README.txt %s' % (result_dir+"/"+str(result_order)+'.DiffExprAnalysis/'+str(result_order)+'.1.DEGsList/'+str(result_order)+'.1.DEGsList.README.txt'))
'''
if flag_results == True:
	if flag_repeat == False:
		code+='''
for each in compare_names:
        assert not os.system('cp %s %s' % (root_dir+'/8.Diff_analysis/Diff/'+each+"/"+each+'.Differential_analysis_results.xls',each+'.Differential_analysis_results.xls'))
        assert not os.system('cp %s %s' % (root_dir+'/8.Diff_analysis/Diff/'+each+"/"+each+'.DEG_up.xls',each+'.DEG_up.xls'))
        assert not os.system('cp %s %s' % (root_dir+'/8.Diff_analysis/Diff/'+each+"/"+each+'.DEG_down.xls',each+'.DEG_down.xls'))
        assert not os.system('cp %s %s' % (root_dir+'/8.Diff_analysis/Diff/'+each+"/"+each+'.DEG.xls',each+'.DEG.xls'))
        assert not os.system('cp %s %s' % (root_dir+'/8.Diff_analysis/Diff/'+each+"/"+each+'.DEGlist_up.txt',each+'.DEG_list_up.txt'))
        assert not os.system('cp %s %s' % (root_dir+'/8.Diff_analysis/Diff/'+each+"/"+each+'.DEGlist_down.txt',each+'.DEG_list_down.txt'))
        assert not os.system('cp %s %s' % (root_dir+'/8.Diff_analysis/Diff/'+each+"/"+each+'.DEGlist.txt',each+'.DEG_list.txt'))
        assert not os.system('cp %s %s' % (root_dir+'/8.Diff_analysis/Diff/'+each+"/"+each+'.gene.xls',each+'.DEG_gene.xls'))
        assert not os.system('cp %s %s' % (root_dir+'/8.Diff_analysis/Diff/'+each+"/"+each+'.lncRNA.genelist.xls',each+'.DEG_lncRNA_genelist.xls'))
        assert not os.system('cp %s %s' % (root_dir+'/8.Diff_analysis/Diff/'+each+"/"+each+'.mRNA.genelist.xls',each+'.DEG_mRNA_genelist.xls'))
assert not os.system('cp /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Readme/DiffList.README.txt %s' % (result_dir+"/"+str(result_order)+'.DiffExprAnalysis/'+str(result_order)+'.1.DEGsList/'+str(result_order)+'.1.DEGsList.README.txt'))
'''
if flag_results == True:
	code+='''
os.chdir('..')

assert not os.system('mkdir '+str(result_order)+'.2.DEGsFilter')
os.chdir(''+str(result_order)+'.2.DEGsFilter')
for each in compare_names:
        assert not os.system('cp %s %s' % (root_dir+'/8.Diff_analysis/Diff/'+each+'/'+each+'.Volcanoplot.pdf',each+'.Volcanoplot.pdf'))
        assert not os.system('cp %s %s' % (root_dir+'/8.Diff_analysis/Diff/'+each+'/'+each+'.Volcanoplot.png',each+'.Volcanoplot.png'))
assert not os.system('cp /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Readme/DEGsFilter.README.txt %s' % (result_dir+"/"+str(result_order)+'.DiffExprAnalysis/'+str(result_order)+'.2.DEGsFilter/'+str(result_order)+'.2.DEGsFilter.README.txt'))
os.chdir('..')

assert not os.system('mkdir '+str(result_order)+'.3.DEGcluster')
os.chdir(''+str(result_order)+'.3.DEGcluster')
assert not os.system('cp %s %s' % (root_dir+'/8.Diff_analysis/Diff/union_for_cluster','DEG_union_for_cluster.xls'))
assert not os.system('cp %s .' % (root_dir+'/8.Diff_analysis/Diff/Hcluster_heatmap.png'))
assert not os.system('cp %s .' % (root_dir+'/8.Diff_analysis/Diff/Hcluster_heatmap.pdf'))
os.system('cp %s .' % (root_dir+'/8.Diff_analysis/Diff/Hcluster_heatmap.detail.pdf'))
assert not os.system('cp -r %s .' % (root_dir+'/8.Diff_analysis/Diff/H_cluster'))
assert not os.system('cp -r %s .' % (root_dir+'/8.Diff_analysis/Diff/SOM_cluster'))
assert not os.system('cp -r %s .' % (root_dir+'/8.Diff_analysis/Diff/K_means_cluster'))
assert not os.system('cp /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Readme/DEGcluster.README.txt %s' % (str(result_order)+'.3.DEGcluster.README.txt'))
os.chdir('..')

assert not os.system('cp  %s .' % (root_dir+'/8.Diff_analysis/Diff/geneid_transcriptid.list'))
assert not os.system('cp  /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Readme/DiffExprAnalysis.README.txt %s' % (str(result_order)+'.DiffExprAnalysis.README.txt'))
if  flag_venn:
	assert not os.system('mkdir '+str(result_order)+'.4.VennDiagram')
	os.chdir(''+str(result_order)+'.4.VennDiagram')
	for each in venn_cluster_vs_names:
		assert not os.system('cp -r %s .' % (root_dir+"/8.Diff_analysis/Diff/"+each+"/"))
	assert not os.system('cp /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Readme/VennDiagram.README.txt %s' %(str(result_order)+'.4.VennDiagram.README.txt'))
os.chdir(result_dir)

result_order+=1
assert not os.system('mkdir '+str(result_order)+'.Target')
os.chdir(''+str(result_order)+'.Target')
assert not os.system('mkdir '+str(result_order)+'.1.cis')
assert not os.system('cp %s %s' % (root_dir+"/9.Target_predict/cis/10k.cis",result_dir+"/"+str(result_order)+'.Target/'+str(result_order)+'.1.cis/'))
assert not os.system('cp %s %s' % (root_dir+"/9.Target_predict/cis/10k.cis.o2",result_dir+"/"+str(result_order)+'.Target/'+str(result_order)+'.1.cis/'))
assert not os.system('cp %s %s' % (root_dir+"/9.Target_predict/cis/100k.cis",result_dir+"/"+str(result_order)+'.Target/'+str(result_order)+'.1.cis/'))
assert not os.system('cp %s %s' % (root_dir+"/9.Target_predict/cis/100k.cis.o2",result_dir+"/"+str(result_order)+'.Target/'+str(result_order)+'.1.cis/'))

if len(samples)>4:
	assert not os.system('mkdir '+str(result_order)+'.2.trans')
	assert not os.system('cut -f1,2,3 %s >%s' % (root_dir+"/9.Target_predict/trans/trans_out/lncRNA_mRNA.cortest.xls",result_dir+"/"+str(result_order)+'.Target/'+str(result_order)+'.2.trans/lncRNA_mRNA.cortest.xls'))
	assert not os.system('cut -f1,2,3 %s >%s' % (root_dir+"/9.Target_predict/trans/trans_out/lncRNA_mRNA.trans.xls",result_dir+"/"+str(result_order)+'.Target/'+str(result_order)+'.2.trans/lncRNA_mRNA.trans.xls'))
	assert not os.system('cp %s %s' % (root_dir+"/9.Target_predict/trans/trans_out/lncRNA_mRNA.trans.list",result_dir+"/"+str(result_order)+'.Target/'+str(result_order)+'.2.trans/lncRNA_mRNA.trans.list'))

assert not os.system('cp /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Readme/Target.README.txt %s' % (result_dir+"/"+str(result_order)+'.Target/'+str(result_order)+".Target.README.txt"))

os.chdir(result_dir)
result_order+=1
assert not os.system('mkdir '+str(result_order)+".mRNA_Enrichment")
os.chdir(""+str(result_order)+".mRNA_Enrichment")
assert not os.system('cp /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Readme/mRNA_Enrichment.README.txt %s' % (result_dir+"/"+str(result_order)+".mRNA_Enrichment/"+str(result_order)+".mRNA_Enrichment.README.txt"))
assert not os.system('mkdir '+str(result_order)+".1.GOEnrichment")
assert not os.system('mkdir '+str(result_order)+".2.KEGGEnrichment")
os.chdir(""+str(result_order)+".1.GOEnrichment")
for each in compare_names:
	assert not os.system('mkdir %s' % (each))
	os.system('cp %s %s' % (root_dir+"/10.Enrichment/mRNA/"+each+"/GO2/"+each+".CAD_GO_enrichment_result.xls",each+"/"+each+".GO_Enrichment_result.xls"))
	os.system('cp %s %s' % (root_dir+"/10.Enrichment/mRNA/"+each+"/GO2/"+each+".CAD_Enriched_GO_classification.pdf",each+"/"+each+".Enriched_GO_classification.pdf"))
	os.system('cp %s %s' % (root_dir+"/10.Enrichment/mRNA/"+each+"/GO2/"+each+".CAD_Enriched_GO_classification.png",each+"/"+each+".Enriched_GO_classification.png"))
	os.system('cp %s %s' % (root_dir+"/10.Enrichment/mRNA/"+each+"/GO2/"+each+".CAD_Enriched_GO_bp_DAG.pdf",each+"/"+each+".Enriched_GO_bp_DAG.pdf"))
	os.system('cp %s %s' % (root_dir+"/10.Enrichment/mRNA/"+each+"/GO2/"+each+".CAD_Enriched_GO_bp_DAG.png",each+"/"+each+".Enriched_GO_bp_DAG.png"))
	os.system('cp %s %s' % (root_dir+"/10.Enrichment/mRNA/"+each+"/GO2/"+each+".CAD_Enriched_GO_cc_DAG.pdf",each+"/"+each+".Enriched_GO_cc_DAG.pdf"))
	os.system('cp %s %s' % (root_dir+"/10.Enrichment/mRNA/"+each+"/GO2/"+each+".CAD_Enriched_GO_cc_DAG.png",each+"/"+each+".Enriched_GO_cc_DAG.png"))
	os.system('cp %s %s' % (root_dir+"/10.Enrichment/mRNA/"+each+"/GO2/"+each+".CAD_Enriched_GO_mf_DAG.pdf",each+"/"+each+".Enriched_GO_mf_DAG.pdf"))
	os.system('cp %s %s' % (root_dir+"/10.Enrichment/mRNA/"+each+"/GO2/"+each+".CAD_Enriched_GO_mf_DAG.png",each+"/"+each+".Enriched_GO_mf_DAG.png"))
	os.system('cp %s %s' % (root_dir+"/10.Enrichment/mRNA/"+each+"/GO2/"+each+".CAD_Enriched_GO_classification_gene_count.txt",each+"/"+each+".Enriched_GO_classification_gene_count.txt"))
os.chdir("..")
os.chdir(""+str(result_order)+".2.KEGGEnrichment")
for each in compare_names:
	assert not os.system('mkdir %s' % (each))
	os.system('cp %s %s' % (root_dir+"/10.Enrichment/mRNA/"+each+"/Pathway/"+each+".DEG_enriched_KEGG_pathway_scatterplot.pdf",each+"/"+each+".Enriched_KEGG_Pathway_Scatterplot.pdf"))
	os.system('cp %s %s' % (root_dir+"/10.Enrichment/mRNA/"+each+"/Pathway/"+each+".DEG_enriched_KEGG_pathway_scatterplot.png",each+"/"+each+".Enriched_KEGG_Pathway_Scatterplot.png"))
	os.system('cp %s %s' % (root_dir+"/10.Enrichment/mRNA/"+each+"/Pathway/add."+each+".indentify.xls_rendered_html_detail.html",each+"/"+each+".KEGG.html"))
	os.system('cp %s %s' % (root_dir+"/10.Enrichment/mRNA/"+each+"/Pathway/add."+each+".indentify.xls",each+"/"+each+".DEG_KEGG_pathway_enrichment_result.xls"))
	os.system('cp -r %s %s' % (root_dir+"/10.Enrichment/mRNA/"+each+"/Pathway/src",each+"/"))
	os.system('cp %s %s' % (root_dir+"/10.Enrichment/mRNA/"+each+"/Pathway/top_20.add."+each+".indentify.xls",each+"/"+each+".enriched_KEGG_pathway_top20.xls"))
os.chdir(result_dir)


result_order+=1
assert not os.system('mkdir '+str(result_order)+".mRNA_ProteinNetwork")
os.chdir(''+str(result_order)+".mRNA_ProteinNetwork")
assert not os.system('cp /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Readme/mRNA_ProteinNetwork.README.txt %s' % (str(result_order)+".mRNA_ProteinNetwork.README.txt"))
assert not os.system('cp /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Readme/CytoscapeQuickStart.pdf .')
for each in compare_names:
	assert not os.system('cp %s .' % (root_dir+"/11.PPI/mRNA/"+each+"/"+each+".ppi.txt"))
os.chdir(result_dir)

result_order+=1
assert not os.system('mkdir '+str(result_order)+".lncRNA_Enrichment")
os.chdir(''+str(result_order)+".lncRNA_Enrichment")
assert not os.system('cp /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Readme/lncRNA_Enrichment.README.txt %s' % (result_dir+"/"+str(result_order)+".lncRNA_Enrichment/"+str(result_order)+".lncRNA_Enrichment.README.txt"))
assert not os.system('mkdir '+str(result_order)+".1.cis")
os.chdir(str(result_order)+".1.cis")
assert not os.system('mkdir '+str(result_order)+".1.1.GOEnrichment")
assert not os.system('mkdir '+str(result_order)+".1.2.KEGGEnrichment")
os.chdir(""+str(result_order)+".1.1.GOEnrichment")
for each in compare_names:
	assert not os.system('mkdir %s' % (each))
	os.system('cp %s %s' % (root_dir+"/10.Enrichment/lncRNA/cis/"+each+"/GO2/"+each+".CAD_GO_enrichment_result.xls",each+"/"+each+".GO_Enrichment_result.xls"))
	os.system('cp %s %s' % (root_dir+"/10.Enrichment/lncRNA/cis/"+each+"/GO2/"+each+".CAD_Enriched_GO_classification.pdf",each+"/"+each+".Enriched_GO_classification.pdf"))
	os.system('cp %s %s' % (root_dir+"/10.Enrichment/lncRNA/cis/"+each+"/GO2/"+each+".CAD_Enriched_GO_classification.png",each+"/"+each+".Enriched_GO_classification.png"))
	os.system('cp %s %s' % (root_dir+"/10.Enrichment/lncRNA/cis/"+each+"/GO2/"+each+".CAD_Enriched_GO_bp_DAG.pdf",each+"/"+each+".Enriched_GO_bp_DAG.pdf"))
	os.system('cp %s %s' % (root_dir+"/10.Enrichment/lncRNA/cis/"+each+"/GO2/"+each+".CAD_Enriched_GO_bp_DAG.png",each+"/"+each+".Enriched_GO_bp_DAG.png"))
	os.system('cp %s %s' % (root_dir+"/10.Enrichment/lncRNA/cis/"+each+"/GO2/"+each+".CAD_Enriched_GO_cc_DAG.pdf",each+"/"+each+".Enriched_GO_cc_DAG.pdf"))
	os.system('cp %s %s' % (root_dir+"/10.Enrichment/lncRNA/cis/"+each+"/GO2/"+each+".CAD_Enriched_GO_cc_DAG.pdf",each+"/"+each+".Enriched_GO_cc_DAG.pdf"))
	os.system('cp %s %s' % (root_dir+"/10.Enrichment/lncRNA/cis/"+each+"/GO2/"+each+".CAD_Enriched_GO_mf_DAG.png",each+"/"+each+".Enriched_GO_mf_DAG.png"))
	os.system('cp %s %s' % (root_dir+"/10.Enrichment/lncRNA/cis/"+each+"/GO2/"+each+".CAD_Enriched_GO_mf_DAG.png",each+"/"+each+".Enriched_GO_mf_DAG.png"))
	os.system('cp %s %s' % (root_dir+"/10.Enrichment/lncRNA/cis/"+each+"/GO2/"+each+".CAD_Enriched_GO_classification_gene_count.txt",each+"/"+each+".Enriched_GO_classification_gene_count.txt"))
	assert not os.system('cp %s %s' % (root_dir+"/10.Enrichment/lncRNA/cis/"+each+"/"+each+".difflncRNA-gene.pairs",each+'/'+each+"difflncRNA_mRNA_cis.pairs"))
os.chdir("..")
os.chdir(""+str(result_order)+".1.2.KEGGEnrichment")
for each in compare_names:
	assert not os.system('mkdir %s' % (each))
	os.system('cp %s %s' % (root_dir+"/10.Enrichment/lncRNA/cis/"+each+"/Pathway/"+each+".DEG_enriched_KEGG_pathway_scatterplot.pdf",each+"/"+each+".Enriched_KEGG_Pathway_Scatterplot.pdf"))
	os.system('cp %s %s' % (root_dir+"/10.Enrichment/lncRNA/cis/"+each+"/Pathway/"+each+".DEG_enriched_KEGG_pathway_scatterplot.png",each+"/"+each+".Enriched_KEGG_Pathway_Scatterplot.png"))
	os.system('cp %s %s' % (root_dir+"/10.Enrichment/lncRNA/cis/"+each+"/Pathway/add."+each+".indentify.xls_rendered_html_detail.html",each+"/"+each+".KEGG.html"))
	os.system('cp %s %s' % (root_dir+"/10.Enrichment/lncRNA/cis/"+each+"/Pathway/add."+each+".indentify.xls",each+"/"+each+".DEG_KEGG_pathway_enrichment_result.xls"))
	os.system('cp -r %s %s' % (root_dir+"/10.Enrichment/lncRNA/cis/"+each+"/Pathway/src",each+"/"))
	os.system('cp %s %s' % (root_dir+"/10.Enrichment/lncRNA/cis/"+each+"/Pathway/top_20.add."+each+".indentify.xls",each+"/"+each+".enriched_KEGG_pathway_top20.xls"))
os.chdir(result_dir+"/"+str(result_order)+".lncRNA_Enrichment")

if len(samples)>4:
	assert not os.system('mkdir '+str(result_order)+".2.trans")
	os.chdir(str(result_order)+".2.trans")
	assert not os.system('mkdir '+str(result_order)+".2.1.GOEnrichment")
	assert not os.system('mkdir '+str(result_order)+".2.2.KEGGEnrichment")
	os.chdir(""+str(result_order)+".2.1.GOEnrichment")
	for each in compare_names:
		assert not os.system('mkdir %s' % (each))
		os.system('cp %s %s' % (root_dir+"/10.Enrichment/lncRNA/trans/"+each+"/GO2/"+each+".CAD_GO_enrichment_result.xls",each+"/"+each+".GO_Enrichment_result.xls"))
		os.system('cp %s %s' % (root_dir+"/10.Enrichment/lncRNA/trans/"+each+"/GO2/"+each+".CAD_Enriched_GO_classification.pdf",each+"/"+each+".Enriched_GO_classification.pdf"))
		os.system('cp %s %s' % (root_dir+"/10.Enrichment/lncRNA/trans/"+each+"/GO2/"+each+".CAD_Enriched_GO_classification.png",each+"/"+each+".Enriched_GO_classification.png"))
		os.system('cp %s %s' % (root_dir+"/10.Enrichment/lncRNA/trans/"+each+"/GO2/"+each+".CAD_Enriched_GO_bp_DAG.pdf",each+"/"+each+".Enriched_GO_bp_DAG.pdf"))
		os.system('cp %s %s' % (root_dir+"/10.Enrichment/lncRNA/trans/"+each+"/GO2/"+each+".CAD_Enriched_GO_bp_DAG.png",each+"/"+each+".Enriched_GO_bp_DAG.png"))
		os.system('cp %s %s' % (root_dir+"/10.Enrichment/lncRNA/trans/"+each+"/GO2/"+each+".CAD_Enriched_GO_cc_DAG.pdf",each+"/"+each+".Enriched_GO_cc_DAG.pdf"))
		os.system('cp %s %s' % (root_dir+"/10.Enrichment/lncRNA/trans/"+each+"/GO2/"+each+".CAD_Enriched_GO_cc_DAG.pdf",each+"/"+each+".Enriched_GO_cc_DAG.pdf"))
		os.system('cp %s %s' % (root_dir+"/10.Enrichment/lncRNA/trans/"+each+"/GO2/"+each+".CAD_Enriched_GO_mf_DAG.png",each+"/"+each+".Enriched_GO_mf_DAG.png"))
		os.system('cp %s %s' % (root_dir+"/10.Enrichment/lncRNA/trans/"+each+"/GO2/"+each+".CAD_Enriched_GO_mf_DAG.png",each+"/"+each+".Enriched_GO_mf_DAG.png"))
		os.system('cp %s %s' % (root_dir+"/10.Enrichment/lncRNA/trans/"+each+"/GO2/"+each+".CAD_Enriched_GO_classification_gene_count.txt",each+"/"+each+".Enriched_GO_classification_gene_count.txt"))
		assert not os.system('cp %s %s' % (root_dir+"/10.Enrichment/lncRNA/trans/"+each+"/"+each+".difflncRNA-gene.pairs",each+"/"+each+"difflncRNA_mRNA_trans.pairs"))
	os.chdir("..")
	os.chdir(""+str(result_order)+".2.2.KEGGEnrichment")
	for each in compare_names:
		assert not os.system('mkdir %s' % (each))
		os.system('cp %s %s' % (root_dir+"/10.Enrichment/lncRNA/trans/"+each+"/Pathway/"+each+".DEG_enriched_KEGG_pathway_scatterplot.pdf",each+"/"+each+".Enriched_KEGG_Pathway_Scatterplot.pdf"))
		os.system('cp %s %s' % (root_dir+"/10.Enrichment/lncRNA/trans/"+each+"/Pathway/"+each+".DEG_enriched_KEGG_pathway_scatterplot.png",each+"/"+each+".Enriched_KEGG_Pathway_Scatterplot.png"))
		os.system('cp %s %s' % (root_dir+"/10.Enrichment/lncRNA/trans/"+each+"/Pathway/add."+each+".indentify.xls_rendered_html_detail.html",each+"/"+each+".KEGG.html"))
		os.system('cp %s %s' % (root_dir+"/10.Enrichment/lncRNA/trans/"+each+"/Pathway/add."+each+".indentify.xls",each+"/"+each+".DEG_KEGG_pathway_enrichment_result.xls"))
		os.system('cp -r %s %s' % (root_dir+"/10.Enrichment/lncRNA/trans/"+each+"/Pathway/src",each+"/"))
		os.system('cp %s %s' % (root_dir+"/10.Enrichment/lncRNA/trans/"+each+"/Pathway/top_20.add."+each+".indentify.xls",each+"/"+each+".enriched_KEGG_pathway_top20.xls"))

os.chdir(result_dir)
result_order+=1

assert not os.system('mkdir '+str(result_order)+".lncRNA_ProteinNetwork")
os.chdir(''+str(result_order)+".lncRNA_ProteinNetwork")
assert not os.system('cp /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Readme/lncRNA_ProteinNetwork.README.txt %s' % (str(result_order)+".lncRNA_ProteinNetwork.README.txt"))
assert not os.system('cp /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Readme/CytoscapeQuickStart.pdf .')
assert not os.system('mkdir '+str(result_order)+".1.cis")
for each in compare_names:
	assert not os.system('cp %s %s' % (root_dir+"/11.PPI/lncRNA/cis/"+each+"/"+each+".ppi.txt",str(result_order)+".1.cis"))
if len(samples)>4:
	assert not os.system('mkdir '+str(result_order)+".2.trans")
	for each in compare_names:
		assert not os.system('cp %s %s' % (root_dir+"/11.PPI/lncRNA/trans/"+each+"/"+each+".ppi.txt",str(result_order)+".2.trans"))

os.chdir(result_dir)
result_order+=1
assert not os.system('mkdir '+str(result_order)+".Constructure_Compare")
os.chdir(''+str(result_order)+".Constructure_Compare")
assert not os.system('cp /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Readme/Constructure_Compare.README.txt %s' % (str(result_order)+".Constructure_Compare.README.txt"))
assert not os.system('cp %s .' % (root_dir+"/5.coding_noncoding_filt/step7_compare/length/length_compare.pdf"))
assert not os.system('cp %s .' % (root_dir+"/5.coding_noncoding_filt/step7_compare/length/length_compare.png"))
assert not os.system('cp %s .' % (root_dir+"/5.coding_noncoding_filt/step7_compare/length/lncRNA_length.stat.xls"))
assert not os.system('cp %s .' % (root_dir+"/5.coding_noncoding_filt/step7_compare/length/mRNA_length.stat.xls"))
assert not os.system('cp %s .' % (root_dir+"/5.coding_noncoding_filt/step7_compare/exon/exon_compare.pdf"))
assert not os.system('cp %s .' % (root_dir+"/5.coding_noncoding_filt/step7_compare/exon/exon_compare.png"))
assert not os.system('cp %s .' % (root_dir+"/5.coding_noncoding_filt/step7_compare/exon/lncRNA_exon.stat.xls"))
assert not os.system('cp %s .' % (root_dir+"/5.coding_noncoding_filt/step7_compare/exon/mRNA_exon.stat.xls"))
assert not os.system('cp %s .' % (root_dir+"/5.coding_noncoding_filt/step7_compare/orf/orf_compare.pdf"))
assert not os.system('cp %s .' % (root_dir+"/5.coding_noncoding_filt/step7_compare/orf/orf_compare.png"))
assert not os.system('cp %s .' % (root_dir+"/5.coding_noncoding_filt/step7_compare/orf/lncRNA/lncRNA.orf.length.xls"))
assert not os.system('cp %s .' % (root_dir+"/5.coding_noncoding_filt/step7_compare/orf/mRNA/mRNA.orf.length.xls"))


os.chdir(result_dir)
result_order+=1
assert not os.system('mkdir '+str(result_order)+".Expression_Compare")
os.chdir(''+str(result_order)+".Expression_Compare")
assert not os.system('cp %s .' % (root_dir+"/8.Diff_analysis/Diff/expression_compare/lncRNA_FPKM.xls"))
assert not os.system('cp %s .' % (root_dir+"/8.Diff_analysis/Diff/expression_compare/mRNA_FPKM.xls"))
assert not os.system('cp %s .' % (root_dir+"/8.Diff_analysis/Diff/expression_compare/mRNAvslncRNA.boxplot.pdf"))
assert not os.system('cp %s .' % (root_dir+"/8.Diff_analysis/Diff/expression_compare/mRNAvslncRNA.boxplot.png"))
assert not os.system('cp %s .' % (root_dir+"/8.Diff_analysis/Diff/expression_compare/mRNAvslncRNA.density.pdf"))
assert not os.system('cp %s .' % (root_dir+"/8.Diff_analysis/Diff/expression_compare/mRNAvslncRNA.density.png"))
assert not os.system('cp %s .' % (root_dir+"/8.Diff_analysis/Diff/expression_compare/mRNAvslncRNA.violin.pdf"))
assert not os.system('cp %s .' % (root_dir+"/8.Diff_analysis/Diff/expression_compare/mRNAvslncRNA.violin.png"))
assert not os.system('cp /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Readme/Expression_Compare.README.txt %s' % (str(result_order)+".Expression_Compare.README.txt"))
os.chdir(result_dir)
result_order+=1
assert not os.system('mkdir '+str(result_order)+".Conversation")
os.chdir(''+str(result_order)+".Conversation")
assert not os.system('mkdir '+str(result_order)+".1.Site")
assert not os.system('cp /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Readme/UCSC.pdf %s' % (str(result_order)+".1.Site"))
assert not os.system('cp /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Readme/Site.Readme.txt %s' % (str(result_order)+".1.Site/"+str(result_order)+".1.Site.README.txt"))

if phyloCSF=='yes' or phyloCSF=='y':
	assert not os.system('mkdir '+str(result_order)+".2.Sequence")
	assert not os.system('cp /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Readme/Sequence.Readme.txt %s' % (str(result_order)+".2.Sequence/"+str(result_order)+".2.Sequence.README.txt"))
	assert not os.system('cp /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Readme/name.xls %s' % (str(result_order)+".2.Sequence/lncRNA.score.xls"))
	assert not os.system('cp /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Readme/name.xls %s' % (str(result_order)+".2.Sequence/mRNA.score.xls"))
	assert not os.system("sed -e '1d' %s|cut -f1,2,3,4,6 >>%s" % (root_dir+"/5.coding_noncoding_filt/step7_compare/conserve/lncRNA.score",str(result_order)+".2.Sequence/lncRNA.score.xls"))
	assert not os.system("sed -e '1d' %s|cut -f1,2,3,4,6 >>%s" % (root_dir+"/5.coding_noncoding_filt/step7_compare/conserve/mRNA.score",str(result_order)+".2.Sequence/mRNA.score.xls"))
	assert not os.system('cp %s %s' % (root_dir+"/5.coding_noncoding_filt/step7_compare/conserve/mRNAvslncRNA.cum.png",str(result_order)+".2.Sequence/"))
	assert not os.system('cp %s %s' % (root_dir+"/5.coding_noncoding_filt/step7_compare/conserve/mRNAvslncRNA.cum.pdf",str(result_order)+".2.Sequence/"))
'''
 
code+='''
print project+' report generating...'

import linecache
import glob

report_dir=root_dir+"/"+project+'_lncRNA_result/'+project+'_report'

os.chdir(root_dir)
assert not os.system('mkdir '+project+'_lncRNA_result/'+project+'_report')
# generate the noref html report   (separate the processes from above for clearness)
os.chdir(root_dir)
os.chdir(project+'_lncRNA_result/'+project+'_report')
assert not os.system('cp -r /PUBLIC/source/RNA/lncRNA/lncRNA_v1/report/new/src/ .')
assert not os.system('cp /PUBLIC/source/RNA/lncRNA/lncRNA_v1/report/Novofinder.exe ..')
assert not os.system('cp /PUBLIC/source/RNA/lncRNA/lncRNA_v1/report/Novofinder_manual.pdf ..')
os.chdir(root_dir)
os.chdir(project+'_lncRNA_result/'+project+'_results')
os.system('tree -d -v -L 3 -o ../'+project+'_report/DirectoryTree.html -H ../../../'+project+'_results')
os.system('perl /PUBLIC/source/RNA/lncRNA/lncRNA_v1/report/result_tree.pl ../'+project+'_report/DirectoryTree.html ../'+project+'_report/src/html/DirectoryTree.html')
os.system('rm ../'+project+'_report/DirectoryTree.html')

render={}
render['title']=project
render['flag_qc']=False
render['flag_tophat']=False
render['flag_as']=False
render['flag_snp']=False
render['flag_deu']=False
render['flag_joint']=False
render['flag_phylocsf']=False
render['flag_dupdiffexp']= False
render['flag_uniform']=False
render['flag_mrnanetwork'] = False
render['flag_lncrnaenrich']=False
render['flag_trans']=False
render['flag_comparestr']=False
render['flag_compareexp']=False
render['flag_lnccons']=False
render['flag_network']=False
render['flag_non_dupdiffexp']=False
render['flag_dupdiffexp'] = flag_repeat
if flag_repeat == False:
        render['flag_non_dupdiffexp']= True
##------------------------------QC-----------------------------
html_order=0
render['flag_qc']=True
html_order+=1  #for raw data stat
html_order+=1  #for qc
os.chdir(root_dir)
os.chdir(project+'_lncRNA_result/'+project+'_report')
assert not os.system('mkdir ./src/pictures')
#assert not os.system('mkdir '/result_dir+'/test')
render['figure_qc_1']=[]
render['figure_qc_2']=[]
render['table_qc']=[]

for eachsample in samples:
        assert not os.system('convert -resize 600 %s ./src/pictures/%s.error_rate_distribution.png' % (result_dir+'/1.QC/1.1.ErrorRate/'+eachsample+'.error_rate_distribution.png',eachsample))
        render['figure_qc_1'].append(["'"+'../pictures/'+eachsample+'.error_rate_distribution.png'+"'"])
        assert not os.system('convert -resize 600 %s ./src/pictures/%s.raw_reads_classification.png' % (result_dir+'/1.QC/1.3.ReadsClassification/'+eachsample+'.raw_reads_classification.png',eachsample))
        render['figure_qc_2'].append(["'"+'../pictures/'+eachsample+'.raw_reads_classification.png'+"'"])
#for eachLine in open(result_dir+'/1.QC/1.4.DataTable/datatable.xls'):
#	if (not eachLine.startswith('Sample')) and eachLine.strip() != '':		
#		render['table_qc'].append(eachLine.strip().split())

DataTable = open(result_dir+'/1.QC/1.4.DataTable/datatable.xls','r')
DataTable.next()
for eachLine in DataTable:
	if eachLine: render['table_qc'].append(eachLine.strip().split())


##-----------------------------------------tophat---------------------
render['flag_tophat']=True
render['flag_flow']=True
html_order+=1
render['html_order_tophat']=str(html_order)

temp = linecache.getline(result_dir+'/2.Mapping/2.1.MapStat/MapStat.xls',1)
tophat_sample = temp.strip().split('\\t')[1:]
render['sample_head']='</th><th>'.join(tophat_sample)
render['table_tophat_1']=[]
cont=open(result_dir+'/2.Mapping/2.1.MapStat/MapStat.xls').readlines()
cont.pop(0)
for eachLine in cont:
        if (not eachLine.startswith('Sample')) and eachLine.strip() != '':
                temp=eachLine.split('\\t')
                tophat='</td><td>'.join(temp)
                render['table_tophat_1'].append(tophat)

mappedfeature=linecache.getline(result_dir+'/2.Mapping/2.2.MapReg/all.mappedfeature.xls',1)
render['table_tophat_2']=[]
mappedfeature_sample = mappedfeature.strip().split('\\t')[1:]
render['mappedfeature_head']='</th><th>'.join(mappedfeature_sample)
cont=open(result_dir+'/2.Mapping/2.2.MapReg/all.mappedfeature.xls').readlines()
cont.pop(0)
for eachLine in cont:
        if (not eachLine.startswith('Sample')) and eachLine.strip() != '':
                temp=eachLine.split('\\t')
                tophat='</td><td>'.join(temp)
                render['table_tophat_2'].append(tophat)

render['figure_tophat_density']=[]
render['figure_tophat_classification']=[]
for eachsample in samples:
	assert not os.system('convert -resize 600 %s ./src/pictures/%s.density.png' % (result_dir+'/2.Mapping/2.3.ChrDen/'+eachsample+'.density.png',eachsample))
        render['figure_tophat_density'].append(["'"+'../pictures/'+eachsample+'.density.png'+"'"])
        assert not os.system('convert -resize 600 %s ./src/pictures/%s.Mapped_Region.png' % (result_dir+'/2.Mapping/2.2.MapReg/'+eachsample+'.Mapped_Region.png',eachsample))
        render['figure_tophat_classification'].append(["'"+'../pictures/'+eachsample+'.Mapped_Region.png'+"'"])



###-----------------------RNA-seq advanced QC----------			
render['flag_rnaseqqc']=True
html_order+=1
render['html_order_rnaseqqc']=str(html_order)			
render['figure_rnaseqxg']=[]
assert not os.system('convert -resize 600 %s ./src/pictures/cor_pearson.png' % (result_dir+'/*.AdvancedQC/*.1.DupCorr/cor_pearson.png'))
render['figure_rnaseqxg'].append(["'"+'../pictures/cor_pearson.png'+"'"])

for png in glob.glob(result_dir+'/*AdvancedQC/*.1.DupCorr/*.png'):
        if png !=glob.glob(result_dir+'/*.AdvancedQC/*.1.DupCorr/cor_pearson.png'):
                name=re.search(r'(/.*/)(.*)\.png',png).group(2)
                assert not os.system('convert -resize 600 %s %s' % (png,'./src/pictures/'+name+'.png'))
                render['figure_rnaseqxg'].append(["'"+'../pictures/'+os.path.basename(png)+"'"])

'''
if flag_uniform:
	code+='''
render['flag_uniform']=True
render['figure_RNASeq_mean_coverage']=[]
for eachsample in samples:
        assert not os.system('convert -resize 600 %s ./src/pictures/%s.Mean_coverage_distribution.png' % (result_dir+'/*.AdvancedQC/*MeanCov/'+eachsample+'.Mean_coverage_distribution.png',eachsample))
        render['figure_RNASeq_mean_coverage'].append(["'"+'../pictures/'+eachsample+'.Mean_coverage_distribution.png'+"'"])
'''
code+='''
##------------------------joint cufflinks scripture-----
render['flag_joint']=True
html_order+=1
render['html_order_joint']=str(html_order)
render['table_cufflinks']=[]
render['table_scripture']=[]
k=6
tmpstr=glob.glob(result_dir+'/*.lncRNA_Assembly*')
#print tmpstr
for eachLine in open(tmpstr[0]+'/cufflinks/'+samples[0]+'_cufflinks.gtf'):
        if (not eachLine.startswith('Chr Num')) and eachLine.strip() != '':
                temp=eachLine.split('\\t')
                render['table_cufflinks'].append(temp)
                k-=1
        if k == 0:
                break
k=6
for eachLine in open(glob.glob(result_dir+'/*.lncRNA_Assembly')[0]+'/scripture/'+samples[0]+'_scripture.bed'):
        if (not eachLine.startswith('Chr Num')) and eachLine.strip() != '':
                temp=eachLine.split('\\t')
                render['table_scripture'].append(temp)
                k-=1
        if k == 0:
                break
##-----------------------lncRNA filter---------------------------------
render['flag_filter']=True
html_order+=1
render['html_order_filter']=str(html_order)
assert not os.system('cp %s %s' % (result_dir+'/*.lncRNA_Filter/*Basic_Filter/step_number.png',report_dir+'/src/pictures/step_number.png'))
assert not os.system('cp %s %s' % (result_dir+'/*.lncRNA_Filter/*Basic_Filter/class_code.png',report_dir+'/src/pictures/class_code.png'))

#html_order+=1
#render['html_order_coding_filter']=str(html_order)

render['table_cpc']=[]
k=6
tmpdir= glob.glob(result_dir+'/*.lncRNA_Filter/*Coding_Potentiality_Filter')[0]

cont=open(tmpdir+'/CPC/CPC.result.xls').readlines()
cont.pop(0)
for eachLine in cont:
        if (not eachLine.startswith('transcript_id')) and eachLine.strip() != '':
                temp=eachLine.split('\\t')
                render['table_cpc'].append(temp)
                k-=1
        if k == 0:
                break
				
render['table_cnci']=[]
k=6
cont=open(tmpdir+'/CNCI/CNCI.result.xls').readlines()
cont.pop(0)
for eachLine in cont:
        if (not eachLine.startswith('transcript_id')) and eachLine.strip() != '':
                temp=eachLine.split('\\t')
                render['table_cnci'].append(temp)
                k-=1
        if k == 0:
                break
				
render['table_pfam']=[]
k=6
cont=open(tmpdir+'/PFAM/PFAM.result.xls').readlines()
cont.pop(0)
for eachLine in cont:
        if (not eachLine.startswith('seq id')) and eachLine.strip() != '':
                temp=eachLine.split('\\t')
                render['table_pfam'].append(temp)
                k-=1
        if k == 0:
                break
'''
if phyloCSF=='yes' or phyloCSF=='y':
    code+='''
render['flag_phylocsf'] = True
render['table_phylocsf'] = []
k=6
cont=open(tmpdir+'/phyloCSF/phyloCSF.result.xls').readlines()
cont.pop(0)
for eachLine in cont:
        if (not eachLine.startswith('transcript_id')) and eachLine.strip() != '':
                temp=eachLine.split('\\t')
                render['table_phylocsf'].append(temp)
                k-=1
        if k == 0:
                break
'''
				
code+='''
assert not os.system('cp %s %s' % (result_dir+'/*.lncRNA_Filter/*Coding_Potentiality_Filter/Venn/union.venn.png',report_dir+'/src/pictures/union.venn.png'))

##-----------------------------------lncRNA cuffdiff----------------------------------

html_order+=1
render['html_order_lncRNA_cuffdiff']=str(html_order)

render['table_lncrnafpkm']=[]
lncrnafpkm_temp=linecache.getline(glob.glob(result_dir+'/*.Quatification')[0]+'/lncRNA_FPKM.xls',1)
lncrnafpkm_sample = lncrnafpkm_temp.strip().split('\\t')[1:]
render['lncrnafpkm_head']='</th><th>'.join(lncrnafpkm_sample)
k=6
cont=open(glob.glob(result_dir+'/*.Quatification')[0]+'/lncRNA_FPKM.xls').readlines()
cont.pop(0)
for eachLine in cont:
        if (not eachLine.startswith('transcript_id')) and eachLine.strip() != '':
                temp=eachLine.split('\\t')
		temp2='</td><td>'.join(temp)
                render['table_lncrnafpkm'].append(temp2)
                k-=1
        if k == 0:
                break


##--------------lncRNA TargetGene---------------------------------------------------
render['flag_targetgene']=True
html_order+=1
render['html_order_targetgene']=str(html_order)
render['table_cis']=[]
k=6
cont=open(glob.glob(result_dir+'/*.Target/*.cis/100k.cis')[0]).readlines()
cont.pop(0)
for eachLine in cont:
        if (not eachLine.startswith('lncRNA_geneid')) and eachLine.strip() != '':
                temp=eachLine.split('\\t')
                render['table_cis'].append(temp)
                k-=1
        if k == 0:
                break	
'''

if len(samples)>4:
    code+='''
render['flag_trans'] = True
render['table_trans']=[]
k=6
cont=open(glob.glob(result_dir+'/*.Target/*.trans/lncRNA_mRNA.trans.xls')[0]).readlines()
cont.pop(0)
for eachLine in cont:
        if (not eachLine.startswith('lncRNA_geneid')) and eachLine.strip() != '':
                temp=eachLine.split('\\t')
                render['table_trans'].append(temp)
                k-=1
        if k == 0:
                break
'''
code+='''
## lncrna kegee pathway
render['flag_lncrnaenrich']=True

html_order+=1
render['html_order_lncrnaenrich']=str(html_order)

render['table_lncgo']=[]
k=6 
cont=open(glob.glob(result_dir+'/*.lncRNA_Enrichment/*cis/*.GOEnrichment/')[0]+compare_names[0]+'/'+compare_names[0]+'.GO_Enrichment_result.xls').readlines()
cont.pop(0)
for eachLine in cont:
	if (not eachLine.startswith('GO_accession')) and eachLine.strip() != '':
		temp=eachLine.split('\\t')
		render['table_lncgo'].append(temp)
		k-=1
	if k == 0:
		break

render['figure_lncgodag']=[]
for each in compare_names:
	os.system('convert -resize 600 %s ./src/pictures/%s.lnc_Enriched_GO_bp_DAG.png' % (result_dir+'/*.lncRNA_Enrichment/*cis/*GOEnrichment/'+each+'/'+each+'.Enriched_GO_bp_DAG.png',each))
	render['figure_lncgodag'].append(["'"+'../pictures/'+each+'.lnc_Enriched_GO_bp_DAG.png'+"'"])
	os.system('convert -resize 600 %s ./src/pictures/%s.lnc_Enriched_GO_mf_DAG.png' % (result_dir+'/*.lncRNA_Enrichment/*cis/*GOEnrichment/'+each+'/'+each+'.Enriched_GO_mf_DAG.png',each))
	render['figure_lncgodag'].append(["'"+'../pictures/'+each+'.lnc_Enriched_GO_mf_DAG.png'+"'"])
	os.system('convert -resize 600 %s ./src/pictures/%s.lnc_Enriched_GO_cc_DAG.png' % (result_dir+'/*.lncRNA_Enrichment/*cis/*GOEnrichment/'+each+'/'+each+'.Enriched_GO_cc_DAG.png',each))
	render['figure_lncgodag'].append(["'"+'../pictures/'+each+'.lnc_Enriched_GO_cc_DAG.png'+"'"])

render['figure_lncgobar']=[]
for each in compare_names:
	os.system('convert -resize 600 %s ./src/pictures/%s.lnc_Enriched_GO_classification.png' % (result_dir+'/*.lncRNA_Enrichment/*cis/*GOEnrichment/'+each+'/'+each+'.Enriched_GO_classification.png',each))
	render['figure_lncgobar'].append(["'"+'../pictures/'+each+'.lnc_Enriched_GO_classification.png'+"'"])

##table_lnckegg KEGGEnrichment
render['table_lnckegg']=[]
k=6 
cont=open(glob.glob(result_dir+'/*.lncRNA_Enrichment/*cis/*.KEGGEnrichment/')[0]+compare_names[0]+'/'+compare_names[0]+'.DEG_KEGG_pathway_enrichment_result.xls').readlines()
cont.pop(0)
for eachLine in cont:
    #Illegal= re.search(r'^[-#]',eachLine)
	if not (re.search(r'^[-#]',eachLine))  and eachLine.strip() != '':
		temp=eachLine.split('\\t')
		render['table_lnckegg'].append(temp)
		k-=1
	if k == 0:
		break
# Enriched_KEGG_Pathway_Scatterplot.png		

render['figure_lnckeggscat']=[]
for each in compare_names:
	assert not os.system('convert -resize 600 %s ./src/pictures/%s.cis_Enriched_KEGG_Pathway_Scatterplot.png' % (result_dir+'/*lncRNA_Enrichment/*cis/*.KEGGEnrichment/'+each+'/'+each+'.Enriched_KEGG_Pathway_Scatterplot.png',each))
	render['figure_lnckeggscat'].append(["'"+'../pictures/'+each+'.cis_Enriched_KEGG_Pathway_Scatterplot.png'+"'"])


k=6
kegg_png=[]
sig_map = ''

for eachLine in open(glob.glob(result_dir+'/*lncRNA_Enrichment/*cis/*KEGGEnrichment/')[0]+compare_names[0]+'/'+compare_names[0]+'.DEG_KEGG_pathway_enrichment_result.xls'):
	if not (re.search(r'^[-#]',eachLine))  and eachLine.strip() != '':
		temp=eachLine.split('\\t')
		if sig_map == '':
				sig_map=temp[2].strip()
				sig_map_png=glob.glob(result_dir+'/*lncRNA_Enrichment/*cis/*KEGGEnrichment/')[0]+compare_names[0]+'/src/'+sig_map+'.png'
				if os.path.isfile(sig_map_png):
						sig_map=temp[2].strip()
						sig_map_png=glob.glob(result_dir+'/*lncRNA_Enrichment/*cis/*KEGGEnrichment/')[0]+compare_names[0]+'/src/'+sig_map+'.png'
						kegg_png.append(sig_map_png)
						sig_map = ''
				else:
						sig_map = ''
		k-=1
	if k == 0:
		break


render['figure_lnckeggpathway']=[] 
for png in kegg_png:
        name=re.search(r'(/.*/)(.*)\.png',png).group(2)
        assert not os.system('convert -resize 600 %s %s' % (png,'./src/pictures/'+name+'.png'))
        render['figure_lnckeggpathway'].append(["'"+'../pictures/'+os.path.basename(png)+"'"])
'''
if phyloCSF=='yes' or phyloCSF=='y':
	code+='''
render['flag_lnccons'] = True
html_order+=1
render['html_order_lnccons']=str(html_order)
assert not os.system('cp %s %s' % (result_dir+'/*.Conversation/*.Sequence/mRNAvslncRNA.cum.png',report_dir+'/src/pictures/mRNAvslncRNA.cum.png'))

		
##-----------------------------------alternative splicing----------------------------------

render['flag_as']=True

html_order+= 1
render['html_order_as']=str(html_order)

assert not os.system('cp %s %s' % (result_dir+'/3.AS/AS.png',report_dir+'/src/pictures/AS.png'))

render['table_as_fpkm']=[]
k=4
cont=open(result_dir+'/3.AS/'+samples[0]+'.fpkm.xls').readlines()
cont.pop(0)
for eachLine in cont:
        if (not eachLine.startswith('event_id')) and eachLine.strip() != '':
                temp=eachLine.split('\\t')
                render['table_as_fpkm'].append(temp)
                k-=1
        if k == 0:
                break


##----------------------snp Indel------------------------


render['flag_snp']=True
html_order+=1
render['html_order_snp']=str(html_order)

SNP_temp=linecache.getline(result_dir+'/4.SNP_InDel/SNPs.xls',1)
SNP_sample = SNP_temp.strip().split('\\t')[4:]
render['SNP_head']='</th><th>'.join(SNP_sample)
render['table_snp']=[]
k=4
cont=open(result_dir+'/4.SNP_InDel/SNPs.xls').readlines()
cont.pop(0)
for eachLine in cont:
        if (not eachLine.startswith('#')) and eachLine.strip() != '':
                temp=eachLine.split('\\t')
                snp='</td><td>'.join(temp)
                render['table_snp'].append(snp)
                k-=1
        if k == 0:
                break
'''
##------------------------------DEU-------------------------
if DEXseq=='yes' or DEXseq=='y':
    code+='''
render['flag_deu'] = True
com=compare_names[0].split('vs')
render['log2fold_name']=[]
temp = linecache.getline(result_dir+'/5.DEU/'+DEXseq_compare_names[0]+'.DEU/'+ DEXseq_compare_names[0]+'.DEUdiff.xls',1)
render['log2fold_name'] = temp.strip().split('\\t')[-1]

html_order+=1
render['html_order_deu']=str(html_order)
render['table_deu']=[]
k=4
for eachLine in open(result_dir+'/5.DEU/'+DEXseq_compare_names[0]+'.DEU/'+ DEXseq_compare_names[0]+'.DEUdiff.xls'):
        if (not eachLine.startswith('geneID')) and eachLine.strip() != '':
                temp=eachLine.split('\\t')
                render['table_deu'].append(temp)
                k-=1
        if k == 0:
                break

assert not os.system('ls %s |while read a;do echo "$a">>%s;exit;done' % (result_dir+'/5.DEU/'+DEXseq_compare_names[0]+'.DEU/files/*expression.svg',report_dir+'/src/pic.log'))
assert not os.system('ls %s |while read a;do echo "$a">>%s;exit;done' % (result_dir+'/5.DEU/'+DEXseq_compare_names[0]+'.DEU/files/*counts.svg',report_dir+'/src/pic.log'))
assert not os.system('ls %s |while read a;do echo "$a">>%s;exit;done' % (result_dir+'/5.DEU/'+DEXseq_compare_names[0]+'.DEU/files/*splicing.svg',report_dir+'/src/pic.log'))
assert not os.system('ls %s |while read a;do echo "$a">>%s;exit;done' % (result_dir+'/5.DEU/'+DEXseq_compare_names[0]+'.DEU/files/*transcripts.svg',report_dir+'/src/pic.log'))

render['figure_deu'] = []
for each_svg in open(report_dir+'/src/pic.log'):
    each_svg = each_svg.rstrip()
    svg = os.path.basename(each_svg)
    svg=svg.strip()
    svg=svg.replace('(','_')
    svg=svg.replace(')','_')
    each_svg=each_svg.replace(')','\)')
    each_svg=each_svg.replace('(','\(')
    m = re.search(r'\w+',svg)
    if m: newSvgName = m.group()
    if '('in svg:
        assert not os.system("convert -resize 600 '%s' %s" % (each_svg,'./src/pictures/'+ newSvgName + '.png'))
    else:
        assert not os.system("convert -resize 600 %s %s" % (each_svg,'./src/pictures/'+ newSvgName + '.png'))
    render['figure_deu'].append(["'"+'../pictures/'+newSvgName+'.png'+"'"])
'''
code+='''

###------------------------cuffdiff----------------
				
render['flag_cuffdiff']=True

html_order+=1
mranfpkm_temp=linecache.getline(glob.glob(result_dir+'/*.Quatification')[0]+'/mRNA_FPKM.xls',1)
mranfpkm_sample = mranfpkm_temp.strip().split('\\t')[1:]
render['mranfpkm_head']='</th><th>'.join(mranfpkm_sample)
render['html_order_mRNA_cuffdiff']=str(html_order)				
render['table_mranfpkm']=[]
k=6
cont=open(glob.glob(result_dir+'/*.Quatification')[0]+'/mRNA_FPKM.xls').readlines()
cont.pop(0)
for eachLine in cont:
        if (not eachLine.startswith('transcript_id')) and eachLine.strip() != '':       
                temp=eachLine.split('\\t')
                temp2='</td><td>'.join(temp)
                render['table_mranfpkm'].append(temp2)
                k-=1
        if k == 0:
                break
				

render['figure_cuffdiff']=[]

assert not os.system('convert -resize 600 %s ./src/pictures/boxplot.png' % (result_dir+'/*.AdvancedQC/*ExpLev/boxplot.png'))
render['figure_cuffdiff'].append(["'"+'../pictures/boxplot.png'+"'"])

assert not os.system('convert -resize 600 %s ./src/pictures/density.png' % (result_dir+'/*.AdvancedQC/*ExpLev/density.png'))
render['figure_cuffdiff'].append(["'"+'../pictures/density.png'+"'"])




###---------------------------DifferentExpression-------------				
render['flag_diffexp']=True
#html_order+=1
#render['html_order_diffexp']=str(html_order)
render['table_diffexp']=[]
com=compare_names[0].split('_vs_')
render['com1']=com[0]
render['com2']=com[1]
k=6 
cont=open(glob.glob(result_dir+'/*.DiffExprAnalysis/*.DEGsList/')[0]+compare_names[0] +'.Differential_analysis_results.xls').readlines()
cont.pop(0)
for eachLine in cont:
        if (not eachLine.startswith('transcript_id')) and eachLine.strip() != '':
                temp=eachLine.split('\\t')
                if flag_repeat == True:
                        render['table_diffexp'].append([temp[0],temp[4],temp[5],temp[6],temp[7],temp[8]])
                else:
                        render['table_diffexp'].append([temp[0],temp[3],temp[4],temp[5],temp[6],temp[7]])
                k-=1
        if k == 0:
                break

render['figure_volcano']=[]
for each in compare_names:
        assert not os.system('convert -resize 600 %s ./src/pictures/%s.Volcanoplot.png' % (result_dir+'/*.DiffExprAnalysis/*DEGsFilter/'+each+'.Volcanoplot.png',each))
        render['figure_volcano'].append(["'"+'../pictures/'+each+'.Volcanoplot.png'+"'"])


assert not os.system('cp %s %s' % (result_dir+'/*DiffExprAnalysis/*DEGcluster/Hcluster_heatmap.png',report_dir+'/src/pictures/Hcluster_heatmap.png'))

if flag_venn:
        render['figure_venn']=[]
        render['flag_venn']=True
        for each in venn_cluster_vs_names:
                assert not os.system('convert -resize 600 %s ./src/pictures/%s.DEG_Venn_diagram.png' % (root_dir+'/8.Diff_analysis/Diff/'+each+'/'+each+'.DEG_Venn_diagram.png',each))
                render['figure_venn'].append(["'"+'../pictures/'+each+'.DEG_Venn_diagram.png'+"'"])
			
##--------------------------------mRNA_Enrichment
render['flag_mrnaenrich']=True
html_order+=1
render['html_order_mrnaenrich']=str(html_order)

render['table_mrnago']=[]
k=6 
cont=open(glob.glob(result_dir+'/*.mRNA_Enrichment/*.GOEnrichment/')[0]+compare_names[0]+'/'+compare_names[0]+'.GO_Enrichment_result.xls').readlines()
cont.pop(0)
for eachLine in cont:
	if (not eachLine.startswith('GO_accession')) and eachLine.strip() != '':
		temp=eachLine.split('\\t')
		render['table_mrnago'].append(temp[:7])
		k-=1
	if k == 0:
		break

render['figure_mrnagodag']=[]
for each in compare_names:
	os.system('convert -resize 600 %s ./src/pictures/%s.mRNA_Enriched_GO_bp_DAG.png' % (result_dir+'/*mRNA_Enrichment/*.GOEnrichment/'+each+'/'+each+'.Enriched_GO_bp_DAG.png',each))
	render['figure_mrnagodag'].append(["'"+'../pictures/'+each+'.mRNA_Enriched_GO_bp_DAG.png'+"'"])
	os.system('convert -resize 600 %s ./src/pictures/%s.mRNA_Enriched_GO_mf_DAG.png' % (result_dir+'/*mRNA_Enrichment/*.GOEnrichment/'+each+'/'+each+'.Enriched_GO_mf_DAG.png',each))
	render['figure_mrnagodag'].append(["'"+'../pictures/'+each+'.mRNA_Enriched_GO_mf_DAG.png'+"'"])
	os.system('convert -resize 600 %s ./src/pictures/%s.mRNA_Enriched_GO_cc_DAG.png' % (result_dir+'/*mRNA_Enrichment/*.GOEnrichment/'+each+'/'+each+'.Enriched_GO_cc_DAG.png',each))
	render['figure_mrnagodag'].append(["'"+'../pictures/'+each+'.mRNA_Enriched_GO_cc_DAG.png'+"'"])

render['figure_mrnagobar']=[]
for each in compare_names:
	os.system('convert -resize 600 %s ./src/pictures/%s.Enriched_GO_classification.png' % (result_dir+'/*mRNA_Enrichment/*.GOEnrichment/'+each+'/'+each+'.Enriched_GO_classification.png',each))
	render['figure_mrnagobar'].append(["'"+'../pictures/'+each+'.Enriched_GO_classification.png'+"'"])


render['table_mrnakegg']=[]
k=6
cont=open(glob.glob(result_dir+'/*.mRNA_Enrichment/*.KEGGEnrichment/')[0]+compare_names[0]+'/'+compare_names[0]+'.DEG_KEGG_pathway_enrichment_result.xls').readlines()
cont.pop(0)
for eachLine in cont:
        if (not eachLine.startswith('#')) and eachLine.strip() != '':
                temp=eachLine.split('\\t')
                render['table_mrnakegg'].append(temp[:7])
        k-=1
        if k == 0:
                break


# Enriched_KEGG_Pathway_Scatterplot.png		
render['figure_mrnakeggscat']=[]
for each in compare_names:
	assert not os.system('convert -resize 600 %s ./src/pictures/%s.mRNA_Enriched_KEGG_Pathway_Scatterplot.png' % (result_dir+'/*mRNA_Enrichment/*.KEGGEnrichment/'+each+'/'+each+'.Enriched_KEGG_Pathway_Scatterplot.png',each))
	render['figure_mrnakeggscat'].append(["'"+'../pictures/'+each+'.mRNA_Enriched_KEGG_Pathway_Scatterplot.png'+"'"])
		

k=6
kegg_png=[]
sig_map = ''

for eachLine in open(glob.glob(result_dir+'/*mRNA_Enrichment/*KEGGEnrichment/')[0]+compare_names[0]+'/'+compare_names[0]+'.DEG_KEGG_pathway_enrichment_result.xls'):
	if not (re.search(r'^[-#]',eachLine))  and eachLine.strip() != '':
		temp=eachLine.split('\\t')
		if sig_map == '':
				sig_map=temp[2].strip()
				sig_map_png=glob.glob(result_dir+'/*mRNA_Enrichment/*KEGGEnrichment/')[0]+compare_names[0]+'/src/'+sig_map+'.png'
				if os.path.isfile(sig_map_png):
						sig_map=temp[2].strip()
						sig_map_png=glob.glob(result_dir+'/*mRNA_Enrichment/*KEGGEnrichment/')[0]+compare_names[0]+'/src/'+sig_map+'.png'
						kegg_png.append(sig_map_png)
						sig_map = ''
				else:
						sig_map = ''
		k-=1
	if k == 0:
		break

render['figure_mrnakeggpathway']=[]
for png in kegg_png:
        name=re.search(r'(/.*/)(.*)\.png',png).group(2)
        assert not os.system('convert -resize 600 %s %s' % (png,'./src/pictures/'+name+'.png'))
        render['figure_mrnakeggpathway'].append(["'"+'../pictures/'+os.path.basename(png)+"'"])
	
## ppi
render['flag_mrnanetwork'] = True
html_order+=1
render['html_order_mRNAppi'] = str(html_order)



####################lncRNA VS mRNA###############################################

render['flag_comparestr']=True
html_order+=1
render['html_order_comparestr']=str(html_order)
assert not os.system('cp %s %s' % (result_dir+'/*.Constructure_Compare/length_compare.png',report_dir+'/src/pictures/length_compare.png'))
assert not os.system('cp %s %s' % (result_dir+'/*.Constructure_Compare/exon_compare.png',report_dir+'/src/pictures/exon_compare.png'))
assert not os.system('cp %s %s' % (result_dir+'/*.Constructure_Compare/orf_compare.png',report_dir+'/src/pictures/orf_compare.png'))

render['flag_compareexp'] = True
html_order+=1
render['html_order_compareexp']=str(html_order)
assert not os.system('cp %s %s' % (result_dir+'/*Expression_Compare/mRNAvslncRNA.violin.png',report_dir+'/src/pictures/mRNAvslncRNA.violin.png'))

render['flag_network']=True
html_order+=1
render['html_order_network']=str(html_order)


t=loader.get_template('src/html/right.html')
c=Context(render)
html=t.render(c)
open('src/html/right.html','w').write(html)

t=loader.get_template('index.html')
c=Context(render)
html=t.render(c)
open(project+'_Report.html','w').write(html)

t=loader.get_template('src/html/top.html')
c=Context(render)
html=t.render(c)
open('src/html/top.html','w').write(html)

assert not os.system("/PUBLIC/software/RNA/wkhtmltopdf/wkhtmltopdf --page-width 350mm --page-height 495mm -n --print-media-type --footer-center '[page] / [topage]' cover src/html/cover.html  toc  --toc-header-text  'content' --toc-text-size-shrink 1  src/html/right.html " +project+'_Report.pdf')
assert not os.system("sed -i 's/#00//g' " +project+'_Report.pdf')

'''

open('lnc_step6_Result_report.py','w').write(code)


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
''' % (fq,fa,gtf,root_dir+'/CAN_TR/CAN/NovelGene/GO/gene.go',sample,root_dir,project)

code+='''
assert not os.system('mkdir data_release')
data_release=open(root_dir+'/data_release/data_release.sh','w')
data_release.write('mkdir %s\\n' % (root_dir+'/data_release/data_give'))
dir=root_dir+'/data_release/data_give'
#data_release.write('echo "###############prepare SuppFiles#####################"\\n')
#data_release.write('mkdir %s\\n' % (root_dir+'/'+project+'_TR_result/'+project+'_results/0.SuppFiles'))
#SuppFiles_dir=root_dir+'/'+project+'_TR_result/'+project+'_results/0.SuppFiles'
#data_release.write('cp /PUBLIC/source/RNA/RefRNA/TransRef/data_release/SuppFiles.README.txt %s\\n' % (SuppFiles_dir))
#data_release.write('cp %s %s\\n' % (genenamefile,SuppFiles_dir+'/gene.description.xls'))
#data_release.write('perl /PUBLIC/source/RNA/RefRNA/TransRef/data_release/GO_class.pl %s /PUBLIC/software/RNA/GOseq/gene_ontology.1_2.obo.ab %s\\n' % (goann,SuppFiles_dir+'/gene.goannot.xls'))
#data_release.write('perl /PUBLIC/source/RNA/RefRNA/TransRef/data_release/extractInfo_from_GTF.pl -i %s -o %s\\n' % (gtf,SuppFiles_dir+'/gene.info.xls'))
#data_release.write('perl /PUBLIC/source/RNA/RefRNA/DGE/Enrichment/extractcDNAfromFA.pl %s %s %s\\n' % (gtf,fa,SuppFiles_dir+'/gene.fasta'))
#data_release.write('gzip %s\\n' % (SuppFiles_dir+'/gene.fasta'))
#data_release.write('echo "###############copy Novofinder#####################"\\n')
#data_release.write('cp /PUBLIC/source/RNA/RefRNA/TransRef/data_release/Novofinder.exe   %s\\n' % (root_dir+'/'+project+'_TR_result/Novofinder.exe'))
#data_release.write('cp /PUBLIC/source/RNA/RefRNA/TransRef/data_release/Novofinder_manual.pdf %s\\n' % (root_dir+'/'+project+'_TR_result'))
data_release.write('echo "###############tar Report#####################"\\n')
data_release.write('cd %s\\n' % (root_dir))
data_release.write("tar -cvzf %s %s\\n" % (project + '_lncRNA.tar.gz',project+'_lncRNA_result/'))
data_release.write('mv %s %s\\n' % (project + '_lncRNA.tar.gz',root_dir+'/data_release/data_give/'))
data_release.write('cd %s\\n' % (root_dir))
data_release.write('echo "###############prepare clean data#####################"\\n')
data_release.write('ln -s %s %s\\n' % (root_dir+'/1.QC/cleandata',dir+'/cleandata'))
#data_release.write('cp %s %s\\n' % (root_dir+'/1.QC/cleandata/cd_md5.txt',dir+'/cleandata/'))
data_release.write('echo "###############prepare raw data#####################"\\n')
data_release.write('mkdir %s\\n' % (dir+'/rawdata'))
data_release.write('cp %s %s\\n' % (root_dir+'/1.QC/rd_md5.txt',dir+'/rawdata/'))

samples=sample.split(',')
for eachsample in samples:
    data_release.write('ln -s %s %s\\n' % (root_dir+'/1.QC/'+eachsample+'/'+eachsample+'_1.fq.gz',dir+'/rawdata/'))
    data_release.write('ln -s %s %s\\n' % (root_dir+'/1.QC/'+eachsample+'/'+eachsample+'_2.fq.gz',dir+'/rawdata/'))

data_release.write('echo "###############prepare IGV_data#####################"\\n')
data_release.write('mkdir %s\\n' % (dir+'/IGV_data'))
data_release.write('cp /PUBLIC/source/RNA/RefRNA/TransRef/data_release/IGV_data.README.txt %s\\n' % (dir+'/IGV_data'))
data_release.write('cp %s %s\\n' % (fa,dir+'/IGV_data/genome.fasta'))
data_release.write('gzip %s\\n' % (dir+'/IGV_data/genome.fasta'))
data_release.write('cp %s %s\\n' % (gtf,dir+'/IGV_data/genome.gtf'))
data_release.write('gzip %s\\n' % (dir+'/IGV_data/genome.gtf'))
data_release.write('ln -s %s %s\\n' % (root_dir+'/1.QC/bam',dir+'/IGV_data/bam'))
data_release.write('echo "###############prepare Readme#####################"\\n')
data_release.write('cp /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Readme/data_release.README.txt %s\\n' % (dir+'/Readme.txt))
data_release.close()
os.chdir('data_release')
assert not os.system('qsub -V -cwd -l vf=1G -l p=1 data_release.sh')
os.chdir(root_dir)
'''

open('lnc_step7_data_release.py','w').write(code)



#################################################################################################
#for byebye
byebye=open(root_dir+'/byebye.sh','w')
byebye.write('sh /PUBLIC/source/RNA/lncRNA/lncRNA_v1/Pipeline/byebye_v1.sh -dir %s -sample %s \n' % (root_dir,sample))
byebye.close()
