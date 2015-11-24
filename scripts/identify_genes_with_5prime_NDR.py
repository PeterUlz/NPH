#! /usr/bin/python

# Identify genes/transcription starts having the 5' nucleosome depleted region (NDR)
#  by analyzing MNase Seq data from Eoncode Project (Gm12878)

from subprocess import call
import sys
import numpy

refgene_file = "ref/refSeq_extended_names_strand.bed"
bigwig = "ref/wgEncodeSydhNsomeGm12878Sig.bigWig"
temp_wig = "intermediate/correlation/wgEncodeSydhNsomeGm12878Sig.tmp.wig"
output_file = "ref/refSeq_TSS_with_5prime_NDR.txt"

# Normalize values ############################################################################################
def normalizeValues(value_list):
    norm_list = list()
    minimum = numpy.amin(value_list)
    maxmimum = numpy.amax(value_list)
    for i in value_list:
        norm_list.append( (i - minimum) / (maxmimum - minimum))
    return norm_list

# Load Values per Position from BedGraph file ###################################################################
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


################################################################################################################
tss_visited = list()

REFGENE = open(refgene_file,"r")
lines = REFGENE.readlines()
REFGENE.close()

OUTPUT = open(output_file,"w")

tss_count = 0
tss_5ndr_count = 0
total_tss = len(lines)
print "Analysis started"

for line in lines:

    if line[:1] == "#":
        continue
    info = line.split("\t")
    if info[1] == "+":
        chrom = info[0]
        tss_pos = info[2]
    if info[1] == "-":
        chrom = info[0]
        tss_pos = info[3]
    tss = chrom+":"+tss_pos
    if tss in tss_visited:
        continue
    else:
        tss_visited.append(tss)

    tss_count += 1
    call(["./scripts/bigWigToBedGraph","-chrom="+chrom,"-start="+str(int(tss_pos)-1100),"-end="+str(int(tss_pos)+1100),bigwig,temp_wig]) 
    values = loadPositionListFromBedGraphFile(temp_wig,int(tss_pos)-1000,int(tss_pos)+1000)
    call(["rm",temp_wig])
    if len(values) < 2000:
        continue
    norm_values = normalizeValues(values)
    cum_values = 0
    for index in range(900,1100):
        cum_values += norm_values[index]
    if (cum_values / 200.0) < 0.2:
        tss_5ndr_count += 1
        OUTPUT.write(info[4])
 
    sys.stdout.write("\r "+str(tss_count)+" TSS visited; "+str(tss_5ndr_count)+" containing 5' NDR; Total TSS: "+str(total_tss)+"; "+str(float(tss_count)/float(total_tss))+" %                   ")
    sys.stdout.flush()

OUTPUT.close()
