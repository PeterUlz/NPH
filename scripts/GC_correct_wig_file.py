#! /usr/bin/python

#GC correct WIG file

import numpy
import statsmodels.api
from subprocess import call
from scipy import interpolate
import argparse
import sys
import os.path

def calculateGcContent(chrom,position,window_size,fasta,wig):
    half_window = int(window_size/2)
    basefasta = os.path.basename(fasta)
    basewig = os.path.basename(wig)
    region = chrom+":"+str(int(position)-half_window)+"-"+str(int(position)+half_window)
    temp_file = "intermediate/GC_correction/"+basewig+"_"+basefasta+"_"+chrom+"_"+position
    TEMP = open(temp_file,"w")
    call(["samtools","faidx",fasta,region],stdout=TEMP)
    TEMP.close()
    GC_bases = 0
    total_bases = 0
    IN_FASTA = open(temp_file,"r")
    for line in IN_FASTA.readlines():
        if line[:1] == ">":
            continue
        tmp_line = line.rstrip().lower()
        gc = tmp_line.count("g")+tmp_line.count("c")
        at = tmp_line.count("a")+tmp_line.count("t")
        GC_bases += gc
        total_bases += (gc+at)
    IN_FASTA.close()
    call(["rm",temp_file])
    return float(GC_bases)/float(total_bases)

# Parse command line arguments ###################################################################################
parser = argparse.ArgumentParser(description='Analyze read depth in comparison to transcription start')
parser.add_argument('-w','--wig', dest='wig_file',
                   help='BAM file',required=True)
parser.add_argument('-o','--output-file', dest='output_file',
                   help='BAM file',required=True)
parser.add_argument('-ws','--window-size', dest='window_size',
                   help='Use this window size for GC correction',default=100,type=int)
parser.add_argument('-gcf','--gc-file', dest='gc_file',
                   help='Optionally specify a file to save GC content and raw/smoothed signal data')
parser.add_argument('-fa','--fasta-file', dest='fasta',
                   help='Use this FastA file for GC calculation',required=True)

args = parser.parse_args()
sys.stderr.write("WIG file: "+args.wig_file+"\n")
sys.stderr.write("FASTA: "+args.fasta+"\n")
sys.stderr.write("Window Size: "+str(args.window_size)+"\n")

try:
    INPUT = open(args.wig_file,"r")
    OUTPUT = open(args.output_file,"w")
except:
    print "Can't open input file"
    sys.exit(1)

if (args.gc_file != ""):
    GC_FILE = open(args.gc_file,"w")

chrom = None
signal_list=list()
position_list=list()
signal_gc = list()
header = ""
varheader = ""
last_position = None
for line in INPUT.readlines():
    if line[:5] == "track":
        header = line
        continue
    elif line[:12] == "variableStep":
        info = line.split()
        chrom = info[1][6:]
        varheader = line
    else:
        info = line.split()
        position = info[0]
        if last_position and int(position) != last_position+1:
            for i in range(last_position+1,int(position)):
                signal_list.append(0)
                gc_content = calculateGcContent(chrom,str(i),args.window_size,args.fasta,args.wig_file)
                signal_gc.append(gc_content)
                position_list.append(str(i))
        signal = info[1].rstrip()
        gc_content = calculateGcContent(chrom,position,args.window_size,args.fasta,args.wig_file)
        position_list.append(position)
        signal_list.append(float(signal))
        signal_gc.append(gc_content)
        last_position = int(position)

lowess_signal = statsmodels.api.nonparametric.lowess(numpy.array(signal_list),numpy.log(numpy.array(signal_gc)),frac=0.05,return_sorted=False)
interpolated_signal = statsmodels.api.nonparametric.lowess(numpy.array(signal_list),numpy.log(numpy.array(signal_gc)),frac=0.05,return_sorted=False)

OUTPUT.write(header.rstrip()+"gc_corrected_"+str(args.window_size)+"\n")
OUTPUT.write(varheader)
for i in range(0,len(signal_list)):
    if (args.gc_file != ""):
        GC_FILE.write(position_list[i]+"\t"+str(signal_list[i])+"\t"+str(signal_gc[i])+"\t"+str(lowess_signal[i])+"\n")
    OUTPUT.write(position_list[i]+"\t"+str(lowess_signal[i])+"\n")

if (args.gc_file != ""):
    GC_FILE.close()
