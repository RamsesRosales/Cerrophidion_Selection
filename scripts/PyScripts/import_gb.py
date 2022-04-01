#!/bin/env python

import argparse

parser = argparse.ArgumentParser(description="Script to download gene bank files")

###############################################

parser.add_argument("-a","--accesion",
					type=str,
					default=None,
					help="genebank accession number")
parser.add_argument("-d","--database",
					type=str,
					default="nucleotide",
					help="genebank data base")
parser.add_argument("-e","--email",
					type=str,
					default="ramsesr@g.clemson.edu",
					help="email")
parser.add_argument("-f","--format",
					type=str,
					default="gb",
					help="output format")
args=parser.parse_args()


if not args.accesion:
    print("""

    ##########################################################################################################
             missing accession number, use flag -a to give a accession number
    ##########################################################################################################

    """)
    quit()

query = args.accesion
database = args.database
email = args.email
format = args.format

print("""

##########################################################################################################""")
print("importing {} file from genebank record with accession number {}".format(format,query))
print("""##########################################################################################################

""")


from Bio import Entrez
Entrez.email = email
handle_w = Entrez.efetch(db= database, id= query, rettype= format, retmode="text")
output_w = handle_w.read()
with open (str( query + "." + format ), "w") as output:
    print(output_w, file=output)

print("file imported")
