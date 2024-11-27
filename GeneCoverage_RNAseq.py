#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 07:47:30 2019

@author: gdiminin
"""

import HTSeq
import pickle
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd


def designGeneCoverage(aQualValue,gene,genesAnnotation_Path,library1,library2,saveLoc):
    def generateExperimentArray(library_Path,aQualValue):
    #open SAM/BAM Files
        aligned_file = HTSeq.SAM_Reader(library_Path)
        
        
        count_aligned = 0
        count_GoodQualityAlignment = 0
        count_total = 0
                
        for algnt in aligned_file: #To verify chromosome nomenclature
            if algnt.aligned:
                if algnt.iv.chrom.startswith('chr'):
                    chromosome_style = ''
                else:
                    chromosome_style = 'chr'
                break
                
        #Create a GenomicArray with all read intervals
        readsLib = HTSeq.GenomicArray('auto', stranded = False)
                
        for algnt in aligned_file:
            if algnt.aligned:
                if algnt.aQual >= aQualValue:
                    for block in algnt.cigar:
                        if block.type == 'M':
                            read = HTSeq.GenomicInterval('%s%s' %(chromosome_style,str(algnt.iv.chrom)),block.ref_iv.start,block.ref_iv.end)
                            readsLib[read]+=1
                    count_GoodQualityAlignment +=1                    
                count_aligned +=1
            count_total +=1
                
        string = '\t-Total reads: %i\n\t-Aligned reads: %i\n\t-Aligned Reads trusted: %i\n' %(count_total,count_aligned,count_GoodQualityAlignment)
        print string
                
        return readsLib
        
    
    
    
    
#Open an annotation file and get the gene that should be plotted

    with open (genesAnnotation_Path,'rb') as load:
        genesAnnotation = pickle.load(load)
    geneInfo = genesAnnotation.loc[gene]
    del genesAnnotation
    toRevert = False
    if geneInfo['genomic_interval'].strand == '-':
        toRevert = True
        
    cond1 = generateExperimentArray(library1[1],aQualValue)
    list1 = [cond1[HTSeq.GenomicPosition(geneInfo['genomic_interval'].chrom,x)] for x in range(geneInfo['genomic_interval'].start,geneInfo['genomic_interval'].end)]
    if toRevert:
        list1 = list1[:-1]
    del cond1

    cond2 = generateExperimentArray(library2[1],aQualValue)
    list2 = [cond2[HTSeq.GenomicPosition(geneInfo['genomic_interval'].chrom,x)] for x in range(geneInfo['genomic_interval'].start,geneInfo['genomic_interval'].end)]
    if toRevert:
        list2 = list2[:-1]
    del cond2

##    with open ('/media/BACKUP_DISK/Asun_FastQ/DesignGeneCoverage/Temporary.pkl','wb') as write:
##        pickle.dump((list1,list2),write)

##    with open ('/media/BACKUP_DISK/Asun_FastQ/DesignGeneCoverage/Temporary.pkl','rb') as load:
##        list1,list2 = pickle.load(load)
        
#Generate the figure panel  

    fig1 = plt.figure()
    fig1.suptitle('%s (%s)'%(gene,geneInfo['genomic_interval']))
    ax1 = fig1.add_subplot(311)    
    ax2 = fig1.add_subplot(312)
    ax3 = fig1.add_subplot(313) 
                
                
        #To define max read value in the two libraries for get ymax position
    xmin = geneInfo['genomic_interval'].start
    xmax = geneInfo['genomic_interval'].end
    ymax = max(list1+list2)
                
        #Design columns for 1st experiment    
    ax1.axis(xmin = xmin, xmax = xmax, ymin = 0, ymax = ymax) #First scatter plot
    ax1.set_ylabel(library1[0]) #Name of the experiment to show in ylabel
    ax1.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        labelbottom='off') # labels along the bottom edge are off
##    if ymax > 500:
##        ax1.set_yscale('log')
##        ax1.axis(xmin =xmin, xmax = xmax, ymin = 1, ymax = ymax)
   
    ax1.plot([x for x in range(xmin,xmax)],list1)



#Design columns for 2nd experiment    
    ax3.axis(xmin =xmin, xmax = xmax, ymin = 0, ymax = ymax) #First scatter plot
    ax3.invert_yaxis()
    ax3.set_ylabel(library2[0]) #Name of the experiment to show in ylabel
    ax3.tick_params(
       axis='x',          # changes apply to the x-axis
       which='both',      # both major and minor ticks are affected
       bottom='off',      # ticks along the bottom edge are off
       top='off',         # ticks along the top edge are off
       labelbottom='off') # labels along the bottom edge are off
#    if ymax > 500:
#        ax3.set_yscale('log')
#        ax3.axis(xmin =xmin, xmax = xmax, ymin = 1, ymax = ymax)
        
    ax3.plot([x for x in range(xmin,xmax)],list2)


#Design gene models
    
    transcripts = geneInfo['variants']
    
    ax2.axis([xmin,xmax,0,len(transcripts)+1])
    ax2.axis('off')
    
    y_value = 0 #location of gene in y axes, acccording to transcripts number
    for transcript in transcripts:
        y_value +=1 #move 1 up
        
        ax2.text((xmax), # Transcript name starting position x-axis (end of gene)
           (y_value-0.2), # Transcript name starting position y-axis
           ('   ' + transcript),fontsize = 10)
                
        #line rapresenting all the transcript length        
        ax2.plot([min([exon.start_d for exon in transcripts[transcript]]), max([exon.end_d for exon in transcripts[transcript]])],[y_value,y_value],'k',linewidth = 2./len(transcripts))
        
        for exon in transcripts[transcript]:
           ax2.add_patch(patches.Rectangle(
                (exon.start,(y_value-0.2)), exon.length,0.4,linewidth = 0.1,facecolor='black'))


    fig1.savefig(saveLoc)


#INFORMATIONS#

# library1= ('HATX8','/media/BACKUP_DISK/Asun_FastQ/HATX8_aligned_R.sa')
# library2= ('HATX3','/media/BACKUP_DISK/Asun_FastQ/HATX3_aligned_R.sa')
# saveLoc = '/media/BACKUP_DISK/Asun_FastQ/DesignGeneCoverage/Hira_HATX8.svg'
#genesAnnotation_Path = '/media/BACKUP_DISK/Asun_FastQ/DesignGeneCoverage/mm10_GenesRef_190412.pkl' # pickel for Ubuntu 14.04
# gene = 'Hira'
aQualValue = 0

#####################

# a = designGeneCoverage(aQualValue,gene,genesAnnotation_Path,library1,library2,saveLoc)

# Automate generation of multiple figures:
libs = list()
"""
libs.append((('HATX8','/media/BACKUP_DISK/Asun_FastQ/HATX8_aligned_R.sa'),'/media/BACKUP_DISK/Asun_FastQ/DesignGeneCoverage/Hira_HATX8.svg','Hira'))
libs.append((('H23','/media/BACKUP_DISK/Asun_FastQ/H23_aligned_R.sa'),'/media/BACKUP_DISK/Asun_FastQ/DesignGeneCoverage/Hira_H23.svg','Hira'))

libs.append((('HATX8','/media/BACKUP_DISK/Asun_FastQ/HATX8_aligned_R.sa'),'/media/BACKUP_DISK/Asun_FastQ/DesignGeneCoverage/Ubn1_HATX8.svg','Ubn1'))
libs.append( (('U9','/media/BACKUP_DISK/Asun_FastQ/U9_aligned_R.sa'),'/media/BACKUP_DISK/Asun_FastQ/DesignGeneCoverage/Ubn1_U9.svg','Ubn1') )
libs.append( (('U14','/media/BACKUP_DISK/Asun_FastQ/U14_aligned_R.sa'),'/media/BACKUP_DISK/Asun_FastQ/DesignGeneCoverage/Ubn1_U14.svg','Ubn1') )
libs.append( (('U1_36','/media/BACKUP_DISK/Asun_FastQ/U1_36_aligned_R.sa'),'/media/BACKUP_DISK/Asun_FastQ/DesignGeneCoverage/Ubn1_U1_36.svg','Ubn1') )
libs.append( (('U12_8','/media/BACKUP_DISK/Asun_FastQ/U12_8_aligned_R.sa'),'/media/BACKUP_DISK/Asun_FastQ/DesignGeneCoverage/Ubn1_U12_8.svg','Ubn1') )
libs.append( (('U12_10','/media/BACKUP_DISK/Asun_FastQ/U12_10_aligned_R.sa'),'/media/BACKUP_DISK/Asun_FastQ/DesignGeneCoverage/Ubn1_U12_10.svg','Ubn1') )

libs.append((('HATX8','/media/BACKUP_DISK/Asun_FastQ/HATX8_aligned_R.sa'),'/media/BACKUP_DISK/Asun_FastQ/DesignGeneCoverage/Ubn2_HATX8.svg','Ubn2'))
libs.append( (('U1_36','/media/BACKUP_DISK/Asun_FastQ/U1_36_aligned_R.sa'),'/media/BACKUP_DISK/Asun_FastQ/DesignGeneCoverage/Ubn2_U1_36.svg','Ubn2') )
libs.append( (('U12_8','/media/BACKUP_DISK/Asun_FastQ/U12_8_aligned_R.sa'),'/media/BACKUP_DISK/Asun_FastQ/DesignGeneCoverage/Ubn2_U12_8.svg','Ubn2') )
libs.append( (('U12_10','/media/BACKUP_DISK/Asun_FastQ/U12_10_aligned_R.sa'),'/media/BACKUP_DISK/Asun_FastQ/DesignGeneCoverage/Ubn2_U12_10.svg','Ubn2') )
# next round
libs.append( (('U9','/media/BACKUP_DISK/Asun_FastQ/U9_aligned_R.sa'),'/media/BACKUP_DISK/Asun_FastQ/DesignGeneCoverage/Hira_U9.svg','Hira') )
libs.append( (('U14','/media/BACKUP_DISK/Asun_FastQ/U14_aligned_R.sa'),'/media/BACKUP_DISK/Asun_FastQ/DesignGeneCoverage/Hira_U14.svg','Hira') )
libs.append( (('U1_36','/media/BACKUP_DISK/Asun_FastQ/U1_36_aligned_R.sa'),'/media/BACKUP_DISK/Asun_FastQ/DesignGeneCoverage/Hira_U1_36.svg','Hira') )
libs.append( (('U12_8','/media/BACKUP_DISK/Asun_FastQ/U12_8_aligned_R.sa'),'/media/BACKUP_DISK/Asun_FastQ/DesignGeneCoverage/Hira_U12_8.svg','Hira') )
libs.append( (('U12_10','/media/BACKUP_DISK/Asun_FastQ/U12_10_aligned_R.sa'),'/media/BACKUP_DISK/Asun_FastQ/DesignGeneCoverage/Hira_U12_10.svg','Hira') )

libs.append((('H23','/media/BACKUP_DISK/Asun_FastQ/H23_aligned_R.sa'),'/media/BACKUP_DISK/Asun_FastQ/DesignGeneCoverage/Ubn1_H23.svg','Ubn1'))
libs.append((('H23','/media/BACKUP_DISK/Asun_FastQ/H23_aligned_R.sa'),'/media/BACKUP_DISK/Asun_FastQ/DesignGeneCoverage/Ubn2_H23.svg','Ubn2'))
libs.append((('HE46','/media/BACKUP_DISK/Asun_FastQ/HE46_aligned_R.sa'),'/media/BACKUP_DISK/Asun_FastQ/DesignGeneCoverage/Ubn1_HE46.svg','Ubn1'))
libs.append((('HE46','/media/BACKUP_DISK/Asun_FastQ/HE46_aligned_R.sa'),'/media/BACKUP_DISK/Asun_FastQ/DesignGeneCoverage/Ubn2_HE46.svg','Ubn2'))

libs.append( (('U9','/media/BACKUP_DISK/Asun_FastQ/U9_aligned_R.sa'),'/media/BACKUP_DISK/Asun_FastQ/DesignGeneCoverage/Ubn1_U9.svg','Ubn1') )
"""
#my_datasets = ("HATX8", "H23", "HE46", "U9", "U14", "U1_36", "U12_8", "U12_10")
#my_gene_list = ("Hira", "Ubn1", "Ubn2")
DATAPATH="/home/linux/Asuns_Data/NEW_RNAseq/"
SVGOUTPUTPATH="/home/linux/Asuns_Data/Asun_FastQ/DesignGeneCoverage/"
genesAnnotation_Path="/home/linux/Asuns_Data/Asun_FastQ/DesignGeneCoverage/mm10_GenesRef_191009.pkl"  # pickel for Ubuntu 16.04

""" Generate a list of SAM files named <sample name>_aligned.sa with: ls *_aligned.sa > sa_files
H18_aligned_R.sa
H18_DOX_aligned_R.sa
HATX12_aligned_R.sa
HATX12_DOX_aligned_R.sa
U1_10_aligned_R.sa
U1_10_DOX_aligned_R.sa
U1_13_aligned_R.sa
U1_13_DOX_aligned_R.sa
U12_27_aligned_R.sa
U12_27_DOX_aligned_R.sa
U12_82_aligned_R.sa
U12_82_DOX_aligned_R.sa
Ue38_aligned_R.sa
Ue38_DOX_aligned_R.sa
The sa_files will then be read from the data folder and the file and sample names constructed.""""

with open(DATAPATH+"sa_files","r") as inFile:
	samples=list()
	for line in inFile:
		samples.append(line.rstrip())
	print samples
my_datasets = [i.split("_aligned")[0] for i in samples]
my_gene_list = ("Hira", "Ubn1", "Ubn2")
for gene_name in my_gene_list:
    for dataset_name in my_datasets:
        sa_file_name = DATAPATH+dataset_name+"_aligned_R.sa"
        svg_file_name = SVGOUTPUTPATH+gene_name+"_"+dataset_name+".svg"
        libs.append( ((dataset_name, sa_file_name), svg_file_name, gene_name) )

i=0
library2= ('HATX3','/home/linux/Asuns_Data/Asun_FastQ/HATX3_aligned_R.sa')
for library1,saveLoc,gene in libs:
    print i, "Library :", library1[0],"; file :", library1[1]
    print "\tGene :", gene, "; output :", saveLoc
    print
    a = designGeneCoverage(aQualValue,gene,genesAnnotation_Path,library1,library2,saveLoc)
    i+=1
print "ended successfully!"









