#! /usr/bin/python

#GC correct WIG file

import numpy
import statsmodels.api
from subprocess import call

def calculateGcContent(chrom,position,window_size,fasta,wig):
    half_window = int(window_size/2)
    region = chrom+":"+str(position-half_window)+"-"+str(position-half_window)
    temp_file = "intermediate/GC_correction/"wig+"_"+fasta+"_"+chrom+"_"+position
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
parser.add_argument('-fa','--fasta-file', dest='fasta',
                   help='Use this FastA file for GC calculation',required=True)

args = parser.parse_args()
sys.stderr.write("WIG file: "+args.wig_file+"\n")
sys.stderr.write("FASTA: "+args.fasta+"\n")
sys.stderr.write("Window Size: "+args.window_size+"\n")

try:
    INPUT = open(args.wig_file,"r")
    OUTPUT = open(args.output_file,"w")
except:
    print "Can't open input file"
    sys.exit(1)

chrom = None
for line in INPUT.readlines():
    if line[:5] == "track":
        OUTPUT.write(line)
    elif line[:12] == "variableStep"
        info = line.split()
        chrom = info[1][6:]
        OUTPUT.write(line)
    else:
        info = line.split()
        position = info[0]
        signal = info[1].rstrip()
        gc_content = calculateGcContent(chrom,position,args.window_size,args.fasta,args.wig_file)
        OUTPUT.write(position+"\t"+signal+"\t"+str(gc_content)+"\n")
