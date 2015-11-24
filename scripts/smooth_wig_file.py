#! /usr/bin/python

# window smooth WIG file

import sys 
import argparse
import numpy

# Parse command line arguments ###################################################################################
parser = argparse.ArgumentParser(description='Apply window smoothing to WIG file data')
parser.add_argument('-wig','--wig-file', dest='wig_file', 
                   help='Wig file should contain data for every position',required=True)
parser.add_argument('-o','--out-file', dest='out_file',
                   help='Outputfile to write smoothed wig file to',required=True)
parser.add_argument('-w','--window-size', dest='windowsize',
                   help='Window size to use for smoothing',default=10,type=int)

args = parser.parse_args()

try:
    INPUT_WIG  = open(args.wig_file,"r")
    OUTPUT_WIG = open(args.out_file,"w")
except:
    print "Error in opening specified files. Aborting!"
    sys.exit(1)

#Assume first line is track info for UCSC Browser
line = INPUT_WIG.readline()
OUTPUT_WIG.write(line.rstrip()+"_smoothed_"+str(args.windowsize)+"\n")

#Assume second line is wig region header
line = INPUT_WIG.readline()
OUTPUT_WIG.write(line)

lines = INPUT_WIG.readlines()
raw_values=list()
raw_positions=list()
smoothed_values=list()
smoothed_positions=list()

for line in lines:
    info = line.split("\t")
    if not info[0].isdigit():
        continue
    else:
        if len(raw_positions) > 1 and int(info[0]) != raw_positions[-1]+1:
            start = raw_positions[-1]+1
            for i in range(start,int(info[0])):
                raw_positions.append(i)
                raw_values.append(0)
                if len(raw_values) >= args.windowsize:
                    OUTPUT_WIG.write( str(raw_positions[-int(args.windowsize/2)])+"\t")
                    OUTPUT_WIG.write( str(numpy.mean(raw_values[-args.windowsize:]))+"\n")
        raw_values.append(int(info[1].rstrip()))
        raw_positions.append(int(info[0]))
        if len(raw_values) >= args.windowsize:
            OUTPUT_WIG.write( str(raw_positions[-int(args.windowsize/2)])+"\t")
            OUTPUT_WIG.write( str(numpy.mean(raw_values[-args.windowsize:]))+"\n")

            
