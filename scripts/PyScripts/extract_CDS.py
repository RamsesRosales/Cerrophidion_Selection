#!/bin/env python
import argparse 

parser = argparse.ArgumentParser(description="""
Script to extract cds from a gen bank file or a fasta file

if you use a genbank file it will take the start and end of the CDS from the file
it only works if there is just one range

if you use a fasta file you can use the flags cs and ce to delimit the CDS

""")

###############################################

parser.add_argument("-i","--input",
					type=str,
					default=None,
					help="input file with")
parser.add_argument("-f","--format",
					type=str,
					default="gb",
					help="format of the file, might be fasta or genebank")
parser.add_argument("-cs","--cds_start",
					type=int,
					default=None,
					help="start of the cds, if several split with a ,")
parser.add_argument("-ce","--cds_end",
					type=int,
					default=None,
					help="end of the cds, if several split with a ,")
parser.add_argument("-o","--output",
					type=int,
					default=None,
					help="name of the output file if none the name of the original file will be used")
args=parser.parse_args()

input = args.input
format = args.format
out = args.output

from Bio import SeqIO
import os
import sys
import subprocess

if not args.input:
    print("""

    ##########################################################################################################
         missing input file, use flag -i to give the path to the input file
    ##########################################################################################################

    """)
    quit()

print("""

##########################################################################################################""")
print("extracting cds from {} ".format(input))
print("""##########################################################################################################

""")



if format == "fasta":
    if not args.cds_start | args.cds_end :
        print("""

        ##########################################################################################################
             missing CDS ranges
        ##########################################################################################################

        """)
        quit()
    cs = args.cds_start -1
    ce = args.cds_end -1
    seq = SeqIO.read(input,"fasta")
    seq.seq = seq.seq[cs:ce]
    if args.output == None:
        out = str(input.split(".fasta")[0])
    out = str(out+"_CDS.fasta")
    with open(out,"w") as output:
        SeqIO.write(seq,output,"fasta")

if format == "gb":
    cmd = str("grep 'CDS' " + input + " | cut -d ' ' -f 19 | perl -pi -e 's/[>,<]//g'")
    print(cmd)
    capture = subprocess.run(cmd,shell=True,capture_output=True)
    capture = str(capture.stdout).split("'")[1].split("\\")[0]
    print(capture)
    cs = int(capture.split(".")[0]) - 1
    ce = int(capture.split(".")[2]) - 1
    seq = SeqIO.read(input,"gb")
    seq.seq = seq.seq[cs:ce]
    if args.output == None:
        out = str(input.split(".gb")[0])
    out = str(out+"_CDS.fasta")
    print(str("writing file " + out))
    with open(out,"w") as output:
        SeqIO.write(seq,output,"fasta")
