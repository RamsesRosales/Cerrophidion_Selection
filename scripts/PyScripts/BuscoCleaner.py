#!/usr/bin/env python3

import sys, os, shutil, errno
import datetime as dt
import argparse
from argparse import RawTextHelpFormatter
import csv
import subprocess as sp
import multiprocessing as mp
import glob
import pyfastx
from dfply import *
from p_tqdm import p_map
try:
	from Bio import SeqIO
	from Bio.Seq import Seq
	from Bio.SeqRecord import SeqRecord
	from Bio.SeqIO import FastaIO
except:
	print("Error: biopython module is not properly installed.")
	quit()

########################################
############### ARGUMENTS ##############
########################################

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="""

BuscoCleaner
1. If no fna files, extract sequences from original assembly
1. Check multi-copy busco sequences for 100% duplicates
	- Because they are 100% duplicates, they are not true multi-copy genes
	- This will result in more single-copy genes
2. Rename sequences and concatenate

Results will have the format:
>locus_name|sample_name

Delimiter | can be changed with the -d flag.

""")

parser.add_argument("-f","--folder",
					type=str,
					help="BUSCO output folder")
parser.add_argument("-n","--name",
					type=str,
					help="Name to apply to sequence")
parser.add_argument("-o","--output",
					type=str,
					help="Optional argument. For concatenating multiple busco runs, specify path here")
parser.add_argument("-d","--delim",
					type=str,
					default='|',
					help="Delimiter used to separate locus names from sample names (default: %(default)s)")
#parser.add_argument("-a","--assembly",
#					type=str,
#					help="Original assembly on which BUSCO was performed \n Only required if fna files were not produced (default for BUSCO v5)")
parser.add_argument("-c","--cpu",
					type=int,
					default=mp.cpu_count(),
					help="Number of threads to be used in each step. (default: %(default)s)")
args=parser.parse_args()

########################################
############### FUNCTIONS ##############
########################################

def mkdir_p(path):
	try:
		os.makedirs(path)
	except OSError as exc: # Python >2.5
		if exc.errno == errno.EEXIST and os.path.isdir(path):
			pass
		else: raise

def Distributor(seq):
	locus=seq.id.split("_")[0]
	contig=seq.id.split("|")[1]
	seq.id = seq.name = seq.description = contig
	if locus in sc_buscos:
		ref_faa = [x.id.split(":")[0] for x in list(SeqIO.parse(sc_folder + "/" + locus + ".faa","fasta"))]
		if contig in ref_faa:
			with open(sc_folder + "/" + locus + ".fna", "a") as output_handle:
				SeqIO.write(seq, output_handle, "fasta")
	if locus in mc_buscos:
		ref_faa = [x.id.split(":")[0] for x in list(SeqIO.parse(mc_folder + "/" + locus + ".faa","fasta"))]
		if contig in ref_faa:
			with open(mc_folder + "/" + locus + ".fna", "a") as output_handle:
				SeqIO.write(seq, output_handle, "fasta")
	if locus in frag_buscos:
		ref_faa = [x.id.split(":")[0] for x in list(SeqIO.parse(frag_folder + "/" + locus + ".faa","fasta"))]
		if contig in ref_faa:
			with open(frag_folder + "/" + locus + ".fna", "a") as output_handle:
				SeqIO.write(seq, output_handle, "fasta")

def mcChecker(locus):
	sp.call("cd-hit-est -i " + locus + " -o " + locus + ".clust.fasta -d 0 -c 1 -M 20000 &> " + locus + ".clust.log", shell=True)
	sequences = list(SeqIO.parse(locus + ".clust.fasta","fasta"))
	if len(sequences)==1:
		handle=open(locus.replace("multi_copy_busco_sequences","single_copy_busco_sequences"), "w")
		writer = FastaIO.FastaWriter(handle, wrap=None)
		writer.write_file(sequences)
		handle.close()
		
		shutil.move(locus, mc_folder+"/Duplicates/")
		shutil.move(locus.replace("fna","faa"), mc_folder+"/Duplicates/")
	
	sp.call("rm " + locus + ".clust*", shell=True)

def Renamer(locus):
	locus_name=locus.split('/')[-1]
	locus_name=locus_name.split(".fna")[0]
	sp.call("sed 's/>.*/>" + locus_name + "\\" + delim + name + "/g' " + locus + " > " + os.path.join(folder,"BuscoCleaner",locus_name) + ".fasta", shell=True)
	sp.call("cat " + os.path.join(args.folder,"BuscoCleaner",locus_name) + ".fasta >> " + args.folder + "/" + name + ".busco.fasta", shell=True)

########################################
################# SETUP ################
########################################

if args.folder==None:
	print("No busco output folder specified. Use -f/--folder flag to specify input")
	quit()
else:
	folder = os.path.abspath(args.folder)

if args.name==None:
	print("No sample name specified. Use -n/--name flag to specify input")
	quit()
else:
	name = args.name

sc_folder = glob.glob(args.folder + "/run_*/busco_sequences/single_copy_busco_sequences")[0]
mc_folder = glob.glob(args.folder + "/run_*/busco_sequences/multi_copy_busco_sequences")[0]
frag_folder = glob.glob(args.folder + "/run_*/busco_sequences/fragmented_busco_sequences")[0]

tmp_count = glob.glob(sc_folder + "/*.fna")

delim = args.delim
cpus = args.cpu

print("""

Clean BUSCO results
1. Check multi-copy busco sequences for 100% duplicates
	- This will result in more single-copy genes
2. Rename sequences and concatenate

Results will have the format:
>locus_name|sample_name

Delimiter | can be changed with the -d flag.

""")

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: starting BuscoCleaner...")
print("\tBUSCO output folder -> "+ args.folder)
print("\tSample name -> "+ args.name)
print("\tDelimiter -> " + delim)
print("\tNumber of CPU -> "+ str(cpus))

########################################
################# CODE #################
########################################

#### CREATING FNA FILES
if len(tmp_count)==0:
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Creating fna files :::")
	 
	metaeuk_output = glob.glob(args.folder + "/run*/metaeuk_output/rerun_results/*.codon.fas")[0]
	metaeuk_output = list(SeqIO.parse(metaeuk_output,"fasta"))
	
	sc_buscos = [x.split("/")[-1].split(".")[0] for x in glob.glob(sc_folder + "/*.faa")]
	mc_buscos = [x.split("/")[-1].split(".")[0] for x in glob.glob(mc_folder + "/*.faa")]
	frag_buscos = [x.split("/")[-1].split(".")[0] for x in glob.glob(frag_folder + "/*.faa")]
	
	p_map(Distributor, metaeuk_output, num_cpus=cpus)

#### CHECKING FOR DUPLICATES IN MULTI-COPY FNA FILES
print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Checking for 100% Duplicates in Multi-Copy Sequence Folder :::")
mc_fnas1 = glob.glob(mc_folder+"/*.fna")
sc_fnas1 = glob.glob(sc_folder+"/*.fna")

mkdir_p(mc_folder + "/Duplicates")
results=p_map(mcChecker, mc_fnas1, num_cpus=cpus)

sc_fnas2 = glob.glob(sc_folder+"/*.fna")
mc_fnas2 = glob.glob(mc_folder+"/*.fna")
dup_count = glob.glob(mc_folder+"/Duplicates/*.fna")

print("Single Copy Busco Loci: "+ str(len(sc_fnas1)))
print("Multi Copy Busco Loci: "+ str(len(mc_fnas1)))
print("::::::::::::::")
print("Duplicates in Multi Copy Busco Loci: "+ str(len(dup_count)))
print("::::::::::::::")
print("New Single Copy Busco Loci: "+ str(len(sc_fnas2)))
print("New Multi Copy Busco Loci: "+ str(len(mc_fnas2)))


#### RENAMING AND CONCATENATING
print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Renaming and Concatenating Loci into " + args.folder + "/BuscoCleaner :::")
mkdir_p(args.folder+"/BuscoCleaner")
results=p_map(Renamer, sc_fnas2, num_cpus=cpus)

