#! /usr/bin/python

#Analyze read depth in comparison to transcription start

import sys
import argparse
from subprocess import call
import numpy
import scipy
import scipy.stats
import os.path
import multiprocessing

# sliding window smoothing ###################################################################################
def sliding_window(data, window_size = 100):
    data_smoothed = list()
    for i in range(0,len(data)-window_size):
        data_smoothed.append(numpy.mean(numpy.array(data[i:i+window_size])))
    return data_smoothed

# Calculate mean and confidence intervals ###################################################################################
def mean_confidence_interval(data, confidence=0.95):
    a = 1.0*numpy.array(data)
    n = len(a)
    m, se = numpy.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1+confidence)/2., n-1)
    return m, m-h, m+h

# Run this for each thread ###################################################################################
def thread_proc(q,thread_number,transcript_list):
    sys.stderr.write("Thread "+str(thread_number)+" started\n")
    sys.stderr.flush()
    coverage_dict=dict()

    #create a list of visited TSS not to count some more than once
    tss_visited = list()

    line_count = 0
    skipped = 0

    #iterate through transcript list from UCSC genome browser
    for transcript in transcript_list:
        coverage_list=dict()
        gene_name = transcript.split("\t")[4].rstrip()
        for i in range(-args.wide_start,0):
            coverage_list[i] = 0
        for i in range(0,args.wide_end+1):
            coverage_list[i] = 0
        line_count += 1
        if (line_count  % 100 == 0):
            sys.stderr.write("\rThread "+str(thread_number)+"\t"+str(line_count)+" genes analyzed")
            sys.stderr.flush()
        #transcription starts are marked at txEnd Field of RefGene Txt file for reverse transcribed genes
        if transcript.split()[1] == '+':
            forward = True
            chrom = transcript.split()[0]
            pos = int(transcript.split()[2])
        else:
            forward = False
            chrom = transcript.split()[0]
            pos = int(transcript.split()[3])
        if chrom+"_"+str(pos) in tss_visited:
            continue
        tss_visited.append(chrom+"_"+str(pos))

        TMP_COVERAGE_BED = open("intermediate/"+os.path.basename(args.bam_file)+str(thread_number)+"tmp_coverage.bed","w")
        call(["samtools","depth","-r",chrom+":"+str(pos-args.wide_start)+"-"+str(pos+args.wide_end),args.bam_file],stdout=TMP_COVERAGE_BED)
        TMP_COVERAGE_BED.close()
    
        TMP_COVERAGE_BED_OUTPUT = open("intermediate/"+os.path.basename(args.bam_file)+str(thread_number)+"tmp_coverage.bed","r")
        content = TMP_COVERAGE_BED_OUTPUT.readlines()
        for i in range(pos-args.wide_start,pos+args.wide_start):
            found = False
            for line in content:
                chrom_found = line.split()[0]
                pos_found = int(line.split()[1])
                if pos_found == i:
                    found = True
                    coverage = int(line.split()[2])
                    if forward:
                        coverage_list[i-pos] = coverage
                    elif not forward:
                        coverage_list[-(i-pos)] = coverage
                    continue
            if not found:
                if forward:
                    coverage_list[i-pos] = 0
                elif not forward:      
                    coverage_list[-(i-pos)] = 0
        TMP_COVERAGE_BED_OUTPUT.close()
        list_values_in_order = list()
        small_values_list = list()
        for i in range(-args.wide_start,0):
            list_values_in_order.append(coverage_list[i])
            if i in range(-args.small_start-50,0):
                small_values_list.append(coverage_list[i])
        for i in range(0,args.wide_end+1):
            list_values_in_order.append(coverage_list[i])
            if i in range(0,args.small_end+50):
                small_values_list.append(coverage_list[i])
        smoothed_data = sliding_window(small_values_list)
        call(["rm","intermediate/"+os.path.basename(args.bam_file)+str(thread_number)+"tmp_coverage.bed"])    
        mean,lower,upper = mean_confidence_interval(coverage_list.values())
        mean_small,lower,upper = mean_confidence_interval(smoothed_data)
        small_ndr = mean_small / mean
        coverage_dict[gene_name] = small_ndr
    
    sys.stderr.write("\rThread "+str(thread_number)+"\t finished\n")
    sys.stderr.flush() 
    q.put(coverage_dict)
#######################################################################################################

# Parse command line arguments ###################################################################################
parser = argparse.ArgumentParser(description='Analyze read depth in comparison to transcription start')
parser.add_argument('-rg','--ref-gene', dest='refgene_file', 
                   help='RefGene file, transcription start should be stored in column 1-3',required=True)
parser.add_argument('-b','--bam', dest='bam_file',
                   help='BAM file',required=True)
parser.add_argument('-ws','--wide-start', dest='wide_start',
                   help='Start analyzing coverage at this point before TSS [default:1000]',default=1000,type=int)
parser.add_argument('-we','--wide-end', dest='wide_end',
                   help='Stop analyzing coverage at this point after TSS [default:1000]',default=1000,type=int)
parser.add_argument('-ss','--small-start', dest='small_start',
                   help='Start analyzing coverage at this point before TSS [default:1000]',default=1000,type=int)
parser.add_argument('-se','--small-end', dest='small_end',
                   help='Stop analyzing coverage at this point after TSS [default:1000]',default=1000,type=int)
parser.add_argument('-t','--threads', dest='threads',
                   help='Threads to use for computation [default:1]',default=1,type=int)
parser.add_argument('-egl','--expr-gene-list', dest='expr_gene_list',
                   help='List of gene names believed to be expressed',required=True)
parser.add_argument('-uegl','--unexpr-gene-list', dest='unexpr_gene_list',
                   help='List of gene names believed to be unexpressed',required=True)


args = parser.parse_args()
sys.stderr.write("Bam file: "+args.bam_file+"\n")
sys.stderr.write("RefGene file: "+args.refgene_file+"\n")
sys.stderr.write("Expressed Genes: "+str(args.expr_gene_list)+"\n")
sys.stderr.write("Unexpressed Genes: "+str(args.unexpr_gene_list)+"\n")
sys.stderr.write("Threads: "+str(args.threads)+"\n")

###############################################################################################
# Analyze data ###################################################################################

try:
    REFGENE = open(args.refgene_file,"r")
except:
    print "Fail to open files specified"
    sys.exit(1)
target_genes = dict()
try:
    EXPR_GENELIST_H = open(args.expr_gene_list,"r")
    for item in EXPR_GENELIST_H.readlines():
        target_genes[item.rstrip()]="expressed"
    EXPR_GENELIST_H.close()
    UNEXPR_GENELIST_H = open(args.unexpr_gene_list,"r")
    for item in UNEXPR_GENELIST_H.readlines():
        target_genes[item.rstrip()]="unexpressed"
    UNEXPR_GENELIST_H.close()
except:
    print "Failed to open genelist"
    sys.exit(1)
    
#filter genes from genelist if specified
header = REFGENE.readline()
refgene_content = REFGENE.readlines()
target_genes_count = 0
target_content = list()
for i in refgene_content:
    chrom = i.split()[0]
    if chrom.find("_") != -1:
        continue
    if chrom.find("X") != -1:
        continue
    if chrom.find("Y") != -1:
        continue
    if i.split()[4].rstrip() in target_genes.keys():
        target_content.append(i)

#initialize input data
gene_count = 0
sys.stderr.write("\n")
sys.stderr.flush()
thread_input_list = dict()
thread_coverage_list = dict()

for thread_number in range(0,args.threads):
    thread_input_list[thread_number] = list()
    thread_coverage_list[thread_number] = dict()

#dispatch input data
max_gene = len(target_content)
partition_point = max_gene/ args.threads
for thread_number in range(0,args.threads):
    thread_input_list[thread_number] = target_content[(thread_number*partition_point):(thread_number+1)*partition_point]

sys.stderr.write(str(len(thread_input_list[0]))+" genes per thread\n")
sys.stderr.write("--------------------------------------------------\n")
sys.stderr.flush()

#start multiple processes
processes = dict()
queues = dict()
for thread in range(0,args.threads):
    queues[thread] = multiprocessing.Queue()
    processes[thread] = multiprocessing.Process(target=thread_proc,args=(queues[thread],thread,thread_input_list[thread]))
    processes[thread].start()

#wait for processes to finish
for thread in range(0,args.threads):
    thread_coverage_list[thread] = queues[thread].get()
    processes[thread].join()

#collect all data
coverage_dict_all=dict()
for thread in range(0,args.threads):
    for i in  thread_coverage_list[thread].keys():
        coverage_dict_all[i] = thread_coverage_list[thread][i]

expr_list = list()
unexpr_list = list()
print "Gene\tType\tCoverage\n"
for gene in coverage_dict_all.keys():
    if gene in target_genes.keys() and target_genes[gene] == "expressed":
        print gene+"\tExpressed\t"+str(coverage_dict_all[gene])
        expr_list.append(coverage_dict_all[gene])
    if gene in target_genes.keys() and target_genes[gene] == "unexpressed":
        print gene+"\tUnexpressed\t"+str(coverage_dict_all[gene])
        unexpr_list.append(coverage_dict_all[gene])

mean_expr,lower_expr,upper_expr = mean_confidence_interval(expr_list)
mean_unexpr,lower_unexpr,upper_unexpr = mean_confidence_interval(unexpr_list)

sys.stderr.write("\nExpressed mean: "+str(mean_expr)+"\t( "+str(lower_expr)+"-"+str(upper_expr)+")")
sys.stderr.write("\nUnexpressed mean: "+str(mean_unexpr)+"\t( "+str(lower_unexpr)+"-"+str(upper_unexpr)+")")
sys.stderr.flush()



 
