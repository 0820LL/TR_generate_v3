import sys
import os
import os.path
import re
import argparse
root_dir=os.getcwd()

##########################################################################
##########################################################################

#parse the arguments
parser = argparse.ArgumentParser(description="DGE pipline v1.0")
parser.add_argument('--project',help='project name, maybe same with the name of root dir, which will be displayed in the final report title, [REQUIRED]',required=True)
parser.add_argument('--sample',help="sample names(sample1,sample2,...)  warning: order sensitive!",required=True)
parser.add_argument('--ss',help="strand specific, [REQUIRED]",choices=['no','yes','reverse'],required=True)
parser.add_argument('--fa',help="the reference FASTA file, [REQUIRED for mapping]",required=True)
parser.add_argument('--gtf',help="the annotation GTF file, [REQUIRED for mapping]",required=True)
parser.add_argument('--group',help="sample classification, e.g. sample1/sample2,sample3, [REQUIRED]",required=True)
parser.add_argument('--groupname',help="group names, [group1,group2,... ]",required=True)
parser.add_argument('--compare',help="group comparison strategy 1:2,1:3,..., [REQUIRED]",required=True)
parser.add_argument('--venn',help="venn and cluster mode, defult is all groups, 1:2_1:3,1:3_2:3,... ",default=None)
parser.add_argument('--goann',help="gene to GO annotations file, [REQUIRED]",required=True)
parser.add_argument('--species',help="abbreviation for species, for kegg pathway. (note: animal is 000, plant is 001, all species: /PROJ/RNA/share/software/kobas2.0-data-20120208/seq_pep)[REQUIRED]",required=True)
parser.add_argument('--ppi_number',help="species code, (ref file: /PROJ/RNA/share/database/string_ppi/species.v9.0.txt)",default=None)
parser.add_argument('--ppi_blast',help="whether to run blast to get the protein-protein interactions",choices=['y','n'],default=None)
parser.add_argument('--genenamefile',help="genenamefile, 1st col is geneID, 2nd col is genename",default=None)

# extract, check the parameters
argv = vars(parser.parse_args())
project=argv['project'].strip()
sample=argv['sample'].strip()
samples=sample.split(',')
samples_tmp=list(set(samples))
assert len(samples)==len(samples_tmp)
#-----------------------------------------------------------------------

ss = argv['ss'].strip()
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
assert len(groupnames)==len(list(set(groupnames)))
for i,each in enumerate(groupnames):
        group_tmp=groups[i]
        group_tmp=group_tmp.replace('/',',')
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
assert os.path.isfile(goann)
species = argv['species'].strip()
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


#################################################################################################
#for Diff analysis by edgeR
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
root_dir='%s'
''' % (sample,ss,gtf,group,groupname,compare,venn,fa,root_dir)

if argv['genenamefile']:
        code+='''
genenamefile='%s'
''' % (genenamefile)

code+='''
##Diff
sam=root_dir+'/QC_DGE/sam'
assert not os.system('mkdir Diff_DGE_edgeR')
os.chdir('Diff_DGE_edgeR')
f=open(root_dir+'/Diff_DGE_edgeR/generate_Diff.sh','w')
f.write('perl /PUBLIC/source/RNA/RefRNA/DGE/DiffExpression/runDiff_analysis_edgeR.ref.pl  -fa %s -sam %s -g %s -o %s -group %s -groupname %s -compare %s -venn %s -spe %s ' % (fa,sam,gtf,root_dir+'/Diff_DGE_edgeR',group,groupname,compare,venn,ss))
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
open('DGE_Diff_edgeR.py','w').write(code)
#################################################################################################
#for Enrichment and PPI
code='''
import os

fa='%s'
goann='%s'
species='%s'
groupname='%s'
compare='%s'
root_dir='%s'
''' % (fa,goann,species,groupname,compare,root_dir)

if argv['ppi_number']:
        code+='''
ppi_number='%s'
ppi_blast='%s'
''' % (ppi_number,ppi_blast)

code+='''
groupnames=groupname.split(',')
compares=compare.split(',')

##GOSeq
assert not os.system('mkdir GOSeq_DGE_edgeR')
os.chdir('GOSeq_DGE_edgeR')
go=open(root_dir+'/GOSeq_DGE_edgeR/runGOSeq.sh','w')
for each in compares:
        temp=each.split(':')
        dir=groupnames[int(temp[0])-1]+'vs'+groupnames[int(temp[1])-1]
        id=root_dir+'/Diff_DGE_edgeR/Diff/'+dir+'/'+dir+'.diffgeneID'
        out=root_dir+'/GOSeq_DGE_edgeR/'+dir
        length=root_dir+'/Diff_DGE/Diff/genelength'
        result=root_dir+'/GOSeq_DGE_edgeR/'+dir+'/'+dir+'.GO_enrichment_result.xls'
        diff=root_dir+'/Diff_DGE_edgeR/Diff/'+dir+'/'+dir+'.diffgene.xls'
        go.write('###############%s#####################\\n' % (dir))
        go.write('mkdir %s\\n' % (root_dir+'/GOSeq_DGE_edgeR/'+dir))
        go.write('perl /PUBLIC/source/RNA/RefRNA/DGE/Enrichment/goseq_graph_v3.pl  -i %s -goann %s -n %s -o %s -length %s\\n' % (id,goann,dir,out,length))
        go.write('perl /PUBLIC/source/RNA/RefRNA/DGE/Enrichment/changeGO_up_down.pl   %s %s %s\\n' % (result,diff,root_dir+'/GOSeq_DGE_edgeR/'+dir+'/'+dir+'.GO_enrichment_result_up_down.xls'))
        go.write('/PUBLIC/software/public/System/R-2.15.3/bin/Rscript /PUBLIC/source/RNA/RefRNA/DGE/Enrichment/goBar.R  %s %s %s\\n' % (result,out,dir))
        go.write('/PUBLIC/software/public/System/R-2.15.3/bin/Rscript /PUBLIC/source/RNA/RefRNA/DGE/Enrichment/goBar2.R  %s %s %s\\n' % (root_dir+'/GOSeq_DGE_edgeR/'+dir+'/'+dir+'.GO_enrichment_result_up_down.xls',out,dir))
go.close()
assert not os.system('qsub -V -cwd -l vf=5G -l p=1 runGOSeq.sh')
os.chdir(root_dir)
'''

if argv['genenamefile'] == None:
        code+='''
assert not os.system('mkdir Blast_DGE_edgeR')
os.chdir('Blast_DGE_edgeR')
f=open(root_dir+'/Blast_DGE_edgeR/runBlast_swissprot.sh','w')
query=root_dir+'/Diff_DGE_edgeR/Diff/diffgene_union.seq'
outdir1=root_dir+'/Blast_DGE_edgeR/Blast_Swissprot/'
out=root_dir+'/Blast_DGE_edgeR/Blast_Swissprot/diffgene_union.seq.blastout'
f.write('echo start blastx\\ndate\\n')
f.write('mkdir %s\\n' % (outdir1))
f.write('/PUBLIC/software/public/Alignment/ncbi-blast-2.2.28+/bin/blastx  -query %s -db /PUBLIC/database/Common/SwissProt/uniprot_sprot.fasta  -evalue 1e-5 -outfmt 5 -max_target_seqs 1 -num_threads 10 -out %s\\n' % (query,out))
f.write('echo blastx end\\ndate\\n')
f.write('perl /PUBLIC/source/RNA/RefRNA/DGE/scriptdir/extractIDsEVxml.pl   %s %s\\n' % (out,outdir1+'/diffgene_union.genenames'))
compare=root_dir+'/Diff_DGE/Diff/compare.txt'
indir=root_dir+'/Diff_DGE_edgeR/Diff/'
outdir2=root_dir+'/Diff_DGE_edgeR/Diff/DiffGeneList'
f.write('mkdir %s\\n' % (outdir2))
f.write('perl /PUBLIC/source/RNA/RefRNA/DGE/scriptdir/getdiffGN.up_down_v2.pl  %s %s %s %s\\n' % (indir,compare,outdir1+'/diffgene_union.genenames',outdir2))
f.close()
assert not os.system('qsub -V -cwd -l vf=5G -l p=10 runBlast_swissprot.sh')
os.chdir(root_dir)
'''

code+='''
##KOBAS
assert not os.system('mkdir KOBAS_DGE_edgeR')
os.chdir('KOBAS_DGE_edgeR')
ko=open(root_dir+'/KOBAS_DGE_edgeR/runKOBAS.sh','w')
path=open(root_dir+'/KOBAS_DGE_edgeR/runpathway.sh','w')
query=root_dir+'/Diff_DGE_edgeR/Diff/diffgene_union.seq'
blastout=root_dir+'/Blast_DGE_edgeR/KOBAS_blast.xml'
ko.write('###############Run Blast#####################\\n')
ko.write('mkdir %s\\n' % (root_dir+'/Blast_DGE_edgeR/'))
ko.write('perl /PUBLIC/source/RNA/RefRNA/DGE/KEGG_Enrichment/bin/KEGG_step1_blast.pl   %s %s %s %s\\n' % (query,species,blastout,root_dir+'/Blast_DGE_edgeR/KOBAS_blast.sh'))
ko.write('sh %s\\n' % (root_dir+'/Blast_DGE_edgeR/KOBAS_blast.sh'))
for each in compares:
        temp=each.split(':')
        dir=groupnames[int(temp[0])-1]+'vs'+groupnames[int(temp[1])-1]
        id=root_dir+'/Diff_DGE_edgeR/Diff/'+dir+'/'+dir+'.diffgeneID'
        out=root_dir+'/KOBAS_DGE_edgeR/'+dir
        ko.write('###############%s#####################\\n' % (dir))
        path.write('echo "###############%s#####################"\\n' % (dir))
        path.write('cd %s\\n' % (out))
        result=root_dir+'/KOBAS_DGE_edgeR/'+dir+'/'+'add.'+dir+'.identify.xls'
        diff=root_dir+'/Diff_DGE_edgeR/Diff/'+dir+'/'+dir+'.diffgene.xls'
        path.write('python /PUBLIC/source/RNA/RefRNA/DGE/KEGG_Enrichment/bin/pathway_annotation_flow_parallel_annotationfault_tolerant3.pyc   --table %s --diff %s\\n' % (result,diff))
        path.write('mv %s %s\\n' % ('add.'+dir+'.indentify.xls_rendered_html_detail.html',dir+'.html'))
        ko.write('mkdir %s\\n' % (root_dir+'/KOBAS_DGE_edgeR/'+dir))
        script=root_dir+'/KOBAS_DGE_edgeR/'+dir+'/run.sh'
        ko.write('perl /PUBLIC/source/RNA/RefRNA/DGE/KEGG_Enrichment/bin/KEGG_step2_enrich.pl -id %s -out-dir %s -species %s -blast-result %s -sample-names %s>%s\\n' % (id,out,species,blastout,dir,script))
        ko.write('sh %s\\n' % (script))
path.close()
ko.close()
assert not os.system('qsub -V -cwd -l vf=5G -l p=6 runKOBAS.sh')
os.chdir(root_dir)
'''

if argv['ppi_number']:
        code+='''
##PPI
assert not os.system('mkdir PPI_DGE_edgeR')
os.chdir('PPI_DGE_edgeR')
ppi=open(root_dir+'/PPI_DGE_edgeR/runPPI.sh','w')
for each in compares:
        temp=each.split(':')
        dir=groupnames[int(temp[0])-1]+'vs'+groupnames[int(temp[1])-1]
        id=root_dir+'/Diff_DGE_edgeR/Diff/'+dir+'/'+dir+'.diffgeneID'
        seq=root_dir+'/Diff_DGE_edgeR/Diff/Diff_Gene_Seq/'+dir+'.diffgene.seq'
        out=root_dir+'/PPI_DGE_edgeR/'+dir
        ppi.write('mkdir %s\\n' % (out))
'''
        if argv['ppi_blast'] == 'y':
                code+='''
        ppi.write('python /PUBLIC/source/RNA/RefRNA/DGE/Enrichment/BLASTX_TO_PPI_v4.py  --species %s --fa %s --outdir %s --name %s\\n' % (ppi_number,seq,out,dir))
'''
        else:
                code+='''
        ppi_dir='/PUBLIC/database/RNA/PPI_ref/PPI_99'
        ppi.write('python /PUBLIC/source/RNA/RefRNA/DGE/Enrichment/get_my_PPI.py  -p %s -g %s -o %s\\n' % (ppi_dir+'/PPI_'+ppi_number+'.txt',id,out+'/'+dir+'.ppi.txt'))
'''
        code+='''
ppi.close()
assert not os.system('qsub -V -cwd -l vf=5G -l p=4 runPPI.sh')
os.chdir(root_dir)
'''
open('DGE_Enrichment_edgeR.py','w').write(code)
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

##path
path=root_dir+'/KOBAS_DGE_edgeR/runpathway.sh'
assert not os.system('sh %s' % (path))

##Report
groupnames=groupname.split(',')
compares=compare.split(',')
compare_names=[]
for each in compares:
        temp=each.split(':')
        compare_names.append(groupnames[int(temp[0])-1]+'vs'+groupnames[int(temp[1])-1])
compare_names=','.join(compare_names)
assert not os.system('mkdir DGE_edgeR')
os.system('sh /PROJ/RNA/share/noRef/Pipeline_noRef/DGE_result_edgeR.sh -dir %s -sample %s -compare %s -title %s -results %s' % (root_dir,sample,compare_names,project,root_dir+'/DGE_edgeR'))
'''
open('DGE_result_edgeR.py','w').write(code)

