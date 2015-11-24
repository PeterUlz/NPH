#! /usr/bin/python

import sys
from subprocess import call
import os.path

run_data_file = sys.argv[1]
offset = sys.argv[2]
threads = sys.argv[3]

ascp_binary_loc = "/home/peter/.aspera/connect/bin/ascp"
ascp_certificate = "/home/peter/.aspera/connect/etc/asperaweb_id_dsa.openssh"

bad_experiment_types=["AMPLICON","Bisulfite-Seq","POOLCLONE","miRNA-Seq"]

RUNDATA = open(run_data_file)

header = RUNDATA.readline()
line = RUNDATA.readline()

#start at certain offset
count = 0
while (count < int(offset)):
   line = RUNDATA.readline()
   count += 1

for i in range(0,10000):
   info = line.split(",")
   run_name = info[0]
   prefix1 = run_name[:3]
   prefix2 = run_name[:6]

   #check if already analyzed
   if os.path.isfile("./MitoGT/"+run_name+".mito"):
      line = RUNDATA.readline()
      continue 

   #skip if experiment type is eother POOLCLONE, AMPLICON, Bisulfite
   bad_type = False
   for experiment in bad_experiment_types:
      if experiment in line:
         bad_type = True
   
   if bad_type:
      print "Skipped run: "+run_name
      line = RUNDATA.readline()
      continue 

   address="anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/"+prefix1+"/"+prefix2+"/"+run_name+"/"+run_name+".sra"
   call([ascp_binary_loc,"-i",ascp_certificate,"-k","1","-T","-l640m",address,"./"])
   call(["./Software/fastq-dump.2.3.4","-X","300000000","-B","--gzip","./"+run_name+".sra"])
   call(["rm","./"+run_name+".sra"])

	# Check whether Download and FastQ Conversion succeeded
   if not os.path.isfile(run_name+".fastq.gz"):
      print "Download or FastQ conversion did not succeed for run: "+run_name
      line = RUNDATA.readline()
      continue 

	# Align to rCRs at first
   SAM = open(run_name+".sam","w")
   call(["./Software/bwa","mem","-t",threads,"./rCRS/rCRS",run_name+".fastq.gz"],stdout=SAM)
   SAM.close()
   call(["rm",run_name+".fastq.gz"])

	# Filter low quality alignments 
   call(["./Software/samtools", "view","-S","-q","15","-h","-o",run_name+".filtered.sam",run_name+".sam"])
   call(["rm",run_name+".sam"])

	# convert filtered SAM fiole to FastQ File and align against hg19 + rCRS
   call(["./Software/sam2fastq.py",run_name+".filtered.sam",run_name+".filtered.fastq"])
   call(["rm",run_name+".filtered.sam"])
   HG19SAM = open(run_name+".hg19.sam","w")
   call(["./Software/bwa","mem","-t",threads,"./rCRS/hg19_rCRS",run_name+".filtered.fastq"],stdout=HG19SAM)
   HG19SAM.close()
   call(["rm",run_name+".filtered.fastq"])

	# Filter only hits to rCRS
   NUMTSAM = open(run_name+".numt_rm.sam","w")
   call(["awk","$3 == \"gi|251831106|ref|NC_012920.1|\"",run_name+".hg19.sam"],stdout=NUMTSAM)
   NUMTSAM.close()
   call(["rm",run_name+".hg19.sam"])
   NUMTHEADERSAM = open(run_name+".numt_rm.header.sam","w")
   call(["cat","./Software/header",run_name+".numt_rm.sam"],stdout=NUMTHEADERSAM)
   NUMTHEADERSAM.close()
   call(["rm",run_name+".numt_rm.sam"])
   call(["./Software/samtools", "view","-S","-h","-b","-o",run_name+".bam",run_name+".numt_rm.header.sam"])
   call(["rm",run_name+".numt_rm.header.sam"])
   call(["./Software/samtools", "sort", run_name+".bam",run_name+".sorted"])
   call(["rm",run_name+".bam"])
   call(["./Software/samtools","index",run_name+".sorted.bam"])

   #create usual pileup file
   MPILEUP = open(run_name+".pileup","w")
   call(["./Software/samtools", "mpileup","-f","./rCRS/rCRS.fa",run_name+".sorted.bam"],stdout=MPILEUP)
   MPILEUP.close()

   #create pileup file and do bcf variant calling
   MPILEUP = open(run_name+".mpileup","w")
   call(["./Software/samtools", "mpileup", "-uf","./rCRS/rCRS.fa",run_name+".sorted.bam"],stdout=MPILEUP)
   MPILEUP.close()
   VARIANT = open("./VCF/"+run_name+".vcf","w")
   call(["./Software/bcftools","view","-vcg",run_name+".mpileup"],stdout=VARIANT)
   VARIANT.close()

   #create consensus sequence file
   BCF = open(run_name+".all.vcf","w")
   call(["./Software/bcftools", "view","-cg", run_name+".mpileup"], stdout=BCF)
   BCF.close()
   CONSENSUS = open("./Seq/"+run_name+".fq","w")
   BCF = open(run_name+".all.vcf","r")
   call(["./Software/vcfutils.pl","vcf2fq"],stdin=BCF,stdout=CONSENSUS)
   BCF.close()
   CONSENSUS.close()


   call(["rm",run_name+".mpileup"])
   call(["rm",run_name+".all.vcf"])

   #process pileupfile and "call" genotypes
   call(["./Software/process_pileup.py",run_name+".pileup","./MitoGT/"+run_name+".mito"])
   call(["rm",run_name+".sorted.bam"])
   call(["rm",run_name+".sorted.bam.bai"])
   call(["rm",run_name+".pileup"])
   line = RUNDATA.readline()
