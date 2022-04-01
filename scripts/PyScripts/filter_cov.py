#!/usr/bin/env python

# ChimeraSorter was made as a followup to ChimeraChecker. It uses several filters to identify and sort filters into good and bad. 
# If you have questions about this script, send a postcard with your information and question addressed to My Ass, You're Screwed Blvd., Imalone, WI 54848

# modify version of this script, I already run bwa picard samtools and bedtools 
# I am only interested in the last part to filter contigs with a coverage of 0 in more that 5% of the lenght

# Additional software necessary to run this:
# (1) bwa 
# (2) gatk
# (3) samtools
# (4) bedtools
# (5) picard
# (6) bam-read from transrate
import time
start = time.time()

import argparse
import sys, os, shutil
import subprocess
import csv
import numpy as np
import pandas as pd

# Set up and extract the command line options.
parser = argparse.ArgumentParser(description='Filter contings based in the bedtools coverage output(not genomecov). Filter by a minimun coverage in a percentage of the lenght')
parser.add_argument("-i","--input",
                    type=str,
                    default="Cgodm-CLP2359_coverage.txt",
                    help="Coverage file, output of bedtools coverage function.")
parser.add_argument("-o","--outfile",
                    type=str,
                    help="name of the output .csv file")
parser.add_argument("-c","--cov",
                    type=int,
                    default=1,
                    help="coverage threshold, default is 1 which mean 0 coverage sites will be discarted")
parser.add_argument("-p","--percentage",
                    type=int,
                    default=0.05,
                    help="percentage of the lenght wich we need to discard a contig")

args = parser.parse_args()

# Read in the coverage data
print("*"*100)
print("Importing coverage information.")
X = pd.read_csv(args.input,sep='\t',names=['transcript','start','end','pos','cov'])
print("Identified coverage data for " + str(len(set(X['transcript']))) + " transcripts.")
T = list(set(X['transcript']))

## check for sequences with less than 5X coverage in greater than 10 percent of sequences
transcript_list = []
for i in T :
	#length = len(list(X[X['transcript'] == T[0]]['cov']))
	length = len(list(X[X['transcript'] == i]['cov']))
	x = 0
	for j in range(0,args.cov):
	#for j in range(0,20):
		#thing = list(X[X['transcript'] == T[0]]['cov']).count(j)
		thing = list(X[X['transcript'] == i]['cov']).count(j)
		x = x + thing
	if x > (args.percentage * length):
		entry = [i, 'absent']
		transcript_list.append(entry)
	else:
		entry = [i, 'present']
		transcript_list.append(entry)
		

name = args.outfile + '.csv'

with open(name, 'w') as csv_file:
	csv_writer = csv.writer(csv_file, delimiter = ',')
	for row in transcript_list:
		csv_writer.writerow(row)
		
csv_file.close()

 
end = time.time()
print(end - start)
