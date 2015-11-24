#! /usr/bin/python

#pipe in BAM file via stdin

import fileinput
import sys

for line in fileinput.input():

    #print header to output if present
    if line[:1] == "@":
        sys.stdout.write(line)
        continue
    info = line.split("\t")
    name = info[0]
    flag = info[1]
    chrom = info[2]
    position = int(info[3])
    mapq = info[4]
    cigar = info[5]
    chrom_paired = info[6]
    pos_paired = info[7]
    mapq_paired = info[8]
    seq = info[9]
    qual = info[10]
    rest = info[11:]
   
    if chrom == "*":
        continue
    if len(seq) < 140:
        continue
    if "N" in cigar or "S" in cigar or "I" in cigar:
        continue

    trim_seq = seq[53:113]
    trim_qual = qual[53:113]

    if int(flag) & 16:
       trim_pos = position + (len(seq) - 113)
    else:
       trim_pos = position + 53

    trim_cigar = "60M"
    sys.stdout.write(name+"\t"+flag+"\t"+chrom+"\t"+str(trim_pos)+"\t"+mapq+"\t"+trim_cigar+"\t"+chrom_paired+
                   "\t"+pos_paired+"\t"+mapq_paired+"\t"+trim_seq+"\t"+trim_qual)
    for part in rest:
        sys.stdout.write("\t"+part)

sys.stdout.flush()
