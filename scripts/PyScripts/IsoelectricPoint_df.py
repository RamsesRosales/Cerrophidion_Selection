#!/bin/env python


import argparse
import sys,os
import subprocess
import shutil
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqIO import FastaIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

#Returns a csv with the isoelectric poin of each sequence in the fasta file
parser = argparse.ArgumentParser(description = """Takes a aminoacid fasta an returns a table\n
with the toxins and the  theoretical isoelectric points""")
parser.add_argument("-f","--fasta",
								type=str,
								help="Translated fasta")
parser.add_argument("-o","--output",
								type=str,
								default="No",
								help="output name")

args = parser.parse_args()
f = args.fasta
o = args.output

if f == None:
	print("Missing fasta")
	quit()
if o == "No":
	print("No output selected, the name of the input will be used for the output")



sequences = list(SeqIO.parse(f,'fasta'))

reading = str("reading "+ f + ": fasta file with " + str(len(sequences)) + " proteins")
print (reading)

ip = []
tox= []
for i in range (1,len(sequences)): #go through fasta
    tox.append(sequences[i].name)
    X = ProteinAnalysis(str(sequences[i].seq))
    ip.append(X.isoelectric_point())

d = {"Tox": tox, "isoelectric_point" : ip}

df = pd.DataFrame(d)

if o == "No":
    pd.DataFrame.to_csv(df,str(f.split(".")[-1] + "_IsoelectricPoint.csv"),index = False)
    print(str( " writing file: "+f.split(".")[-1] + "_IsoelectricPoint.csv"))

if o != "No":
    pd.DataFrame.to_csv(df,str(o + "_IsoelectricPoint.csv"),index = False)
    print(str( " writing file: "+ o + "_IsoelectricPoint.csv"))

