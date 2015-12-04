#! /usr/bin/python

import sys
import argparse

# Parse command line arguments ###################################################################################
parser = argparse.ArgumentParser(description='Combine prediction parameters')
parser.add_argument('-o','--output-file', dest='output_file',
                   help='BAM file',required=True)
parser.add_argument('-bcov','--broad-coverage-file', dest='broad_coverage_file',
                   help='Broad Coverage File',required=True)
parser.add_argument('-slope','--slope-file', dest='slope_file',
                   help='Slope File',required=True)
parser.add_argument('-scov','--small-coverage-file', dest='small_coverage_file',
                   help='Small Coverage File',required=True)
args = parser.parse_args()

try:
    BROAD = open(args.broad_coverage_file,"r")
    SMALL = open(args.small_coverage_file,"r")
    SLOPE = open(args.slope_file,"r")
except:
    print "Can't open one of the input files"
    sys.exit(1)

OUTPUT = open(args.output_file,"w")

broad_cov_dict = dict()
small_cov_dict = dict()
slope_dict = dict()
expressed_dict = dict()

for line in BROAD.readlines():
    if line == "\n" or "Gene" in line:
        continue
    info = line.split("\t")
    expressed_dict[info[0]] = info[1]
    broad_cov_dict[info[0]] = info[2].rstrip()
BROAD.close()

for line in SMALL.readlines():
    if line == "\n" or "Gene" in line:
        continue
    info = line.split("\t")
    small_cov_dict[info[0]] = info[2].rstrip()
SMALL.close()

for line in SLOPE.readlines():
    if line == "\n" or "Gene" in line:
        continue
    info = line.split("\t")
    slope_dict[info[0]] = info[2]+"\t"+info[3].rstrip()
SLOPE.close()

OUTPUT.write("Gene\tType\t5'Slope\t3'Slope\tBroad TSS Coverage\tSmall TSS coverage\n")
for gene in expressed_dict.keys():
    OUTPUT.write(gene+"\t"+expressed_dict[gene]+"\t"+slope_dict[gene]+"\t"+ broad_cov_dict[gene]+"\t"+small_cov_dict[gene]+"\n")
OUTPUT.close()
