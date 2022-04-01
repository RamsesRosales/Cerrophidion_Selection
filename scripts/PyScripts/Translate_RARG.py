#!/bin/env python


import argparse

#Returns a csv with the isoelectric poin of each sequence in the fasta file
parser = argparse.ArgumentParser(description = "translates or a nucleotide fasta to a aminoacide fasta\n make sure your sequences are CDS")
parser.add_argument("-f","--fasta",
								type=str,
								help="Translated fasta")
#parser.add_argument("-o","--output",
#								type=str,
#								default="Protein",
#								help="translate to Protein or Nucleotid")
parser.add_argument("-n","--name",

								type=str,
								default="Translated",
								help="name for the output fasta file")

args = parser.parse_args()
f = args.fasta
o = args.name

if f == None:
	print("Missing fasta")
	quit()
#
#if o != "Protein" and o != "Nucleotide":
#    print("Choose the right otput, the options are Protein of Nucleotide")
#    quit()
#if o == "Protein":
#	print("Translate from nucleotide to protein")
#if o == "Nucleotide":
#    print("Translate from Protein to Nucleotide")


import sys,os
import subprocess
import shutil
from Bio import SeqIO
from Bio.SeqIO import FastaIO
#from Bio import Translate
#
####### functions

#standard_translator = Translate.unambiguous_dna_by_id[1]

def N2P(seq,list):
    x = seq
    x.seq = seq.translate().seq
    list.append(x)

#def P2N(seq,list):
#    x = standard_translator.back_translate(seq)
#    list.append(x)

sequences = list(SeqIO.parse(f,'fasta'))

reading = str("reading "+ f + ": fasta file with " + str(len(sequences)) + " sequences")
print (reading)


o_sequences = list()
for i in range(0,len(sequences)):
    seq = sequences[i]
    N2P(seq,o_sequences)

#o_sequences = list()
#if o == "Nucleotide":
#    for i in range(0,len(sequences)):
#        seq = sequences[i]
#        P2N(seq,o_sequences)

SeqIO.write(o_sequences,str(o  + ".fasta"),"fasta")
