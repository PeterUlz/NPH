#! /usr/bin/python

# Analyze and plot coverage of TSS of genes from gene list
import os
from subprocess import call
import argparse
import sys
import scipy.stats
import numpy
import matplotlib.pyplot as plt
import pandas
##########################################################################################################
# Normalize values ############################################################################################
def normalizeValues(value_list):
    norm_list = list()
    minimum = numpy.amin(value_list)
    maxmimum = numpy.amax(value_list)
    for i in value_list:
        norm_list.append( (i - minimum) / (maxmimum - minimum))
    return norm_list

##########################################################################################################    
# Load Values per Position fro BedGraph file ###################################################################
def loadPositionListFromBedGraphFile(filename,start_count,end_count):
    BEDGRAPH_FILE = open(filename,"r")
    value_list = list()
    wig_content = BEDGRAPH_FILE.readlines()
    bam_value_line_list = list()
    former_end = None
    for line in wig_content:
        if line[:8] == "#":
            continue
        chrom,start,end,signal = line.split("\t")
        for i in range(int(start),int(end)):
            if i >= start_count and i < end_count:
                value_list.append(float(signal))
        if former_end and former_end != start:
            for i in range(int(former_end),int(start)):
                value_list.append(float(0.0))
        former_end = end
    return value_list

##########################################################################################################
# Create Wiggle file from BEDFILE ###################################################
def getCoveragePerPosition(gene,chrom,position,bam_samplename,args):
    # Create wiggle from BamToBigWig
    bam_bigwig_file = "intermediate/single_gene_plots/"+bam_samplename+"."+gene+".bigwig"
    bam_wig_file = "intermediate/single_gene_plots/"+bam_samplename+"."+gene+".wig"
    call(["./scripts/bam_to_wiggle.py","--outfile="+bam_bigwig_file,"--chrom="+chrom, "--start="+str(int(position)-1100),"--end="+str(int(position)+1100),args.bam_file])
    call(["./scripts/bigWigToBedGraph","-chrom="+chrom,"-start="+str(int(position)-1100),"-end="+str(int(position)+1100),bam_bigwig_file,bam_wig_file])

    #parse wiggle file from BAM and normalize
    bam_values=loadPositionListFromBedGraphFile(bam_wig_file,int(position)-1000,int(position)+1000)
    print gene,"Length:",len(bam_values)
    if len(bam_values) != 2000:
        return None
    if args.normalize:
        bam_norm_values = normalizeValues(bam_values)
    else:
        bam_norm_values = bam_values
    return bam_norm_values

##########################################################################################################
# Plot data ###################################################
def plotNormValues(gene_name,norm_values,bam_samplename,outdir):
    bam_norm_values_series = pandas.Series(norm_values)
    bam_norm_values_series.plot()
    output=plt.gcf()
    if outdir[:-1] != "/":
        outdir = outdir+"/"
    plt.show()
    output.savefig(outdir+gene_name+"."+bam_samplename+".png")

##########################################################################################################
# Parse command line arguments ###################################################################################
parser = argparse.ArgumentParser(description='Analyze read depth in comparison to transcription end')
parser.add_argument('-rg','--ref-gene', dest='refgene_file', 
                   help='RefGene file with RefSeq Gene definitons, transcription start should be stored in column 1-3',required=True)
parser.add_argument('-b','--bam', dest='bam_file',
                   help='BAM file',required=True)
parser.add_argument('-gl','--gene-list', dest='genelist',
                   help='Gene list to use for plotting',required=True)
parser.add_argument('-s','--smoothing', dest='genelist',
                   help='Specify window size for sliding window smoothing',default=0,type=int)
parser.add_argument('-o','--output-dir', dest='output_dir',
                   help='Specify output directory for images',default=".")
parser.add_argument('-n','--normalize', dest='normalize',
                   help='Normalize [default:False]',action="store_true")

args = parser.parse_args()
sys.stderr.write("\nGenelist: "+args.genelist)
sys.stderr.write("\nBam file: "+args.bam_file)
sys.stderr.write("\nRefGene file: "+args.refgene_file)
sys.stderr.write("\nPic dir: "+args.output_dir+"\n")
sys.stderr.flush()

#Load gene definitions from RefGene file
refgene_starts = dict()
REFGENE=open(args.refgene_file,"r")
lines = REFGENE.readlines()
for line in lines:
    if line[:1] == "#":
        continue
    info = line.split("\t")
    if "_" in info[0]:
        continue
    if info[1] == "+":
        refgene_starts[info[4][:-1]] = info[0]+":"+info[2]
    elif info[1] == "-":
        refgene_starts[info[4][:-1]] = info[0]+":"+info[3]
    else:
        print "Something went horribly wrong"
        sys.exit(1)

GENELIST=open(args.genelist,"r")
for gene in GENELIST.readlines():
    gene_name = gene.rstrip()
    if gene_name not in refgene_starts.keys():
        print gene_name+" not found in RefGene file"
    else:
        chrom,position = refgene_starts[gene_name].split(":")
        bam_samplename = os.path.basename(args.bam_file)
        bam_norm_values = getCoveragePerPosition(gene_name,chrom,position,bam_samplename,args)
        if not bam_norm_values:
            continue
        plotNormValues(gene_name, bam_norm_values,bam_samplename,args.output_dir)


