#! /usr/bin/python

# Check for high expressed genes in breast cancer

import sys
import numpy

input_file = "ref/TCGA_BRCA_RNASeq/BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt"
output_file = "ref/TCGA_BRCA_RNASeq/BRCA_mean_expr_values.txt"
genes_expr=dict()
INPUT = open(input_file,"r")
OUTPUT = open(output_file,"w")

header1 = INPUT.readline()
header2 = INPUT.readline()
print ""
#read in data line by line to save memory
line = INPUT.readline()
count = 0
while line:
    count += 1
    info = line.split("\t") 
    values = list()
    gene_name = info[0].split("|")[0]
    for value in info[1:]:
        values.append(float(value))
    genes_expr[gene_name] = numpy.mean(values)
  
    if count % 1000 == 0:
        sys.stdout.write("\r"+str(count)+" lines analyzed")
        sys.stdout.flush()
    line = INPUT.readline()

for gene in genes_expr.keys():
    OUTPUT.write(gene+"\t"+str(genes_expr[gene])+"\n")
