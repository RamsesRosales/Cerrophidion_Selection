from Bio import Phylo
import argparse
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="""

nx2nw
convert nexus files to newick (build to run Astral with IQtree output)

""")

parser.add_argument("-i","--input",
					type=str,
					help="input nexus file")

parser.add_argument("-f","--format",
					type=str,
					default="newick",
					help="output format considering that the input is the oposit (nexus to newick for default")
args=parser.parse_args()

if args.input==None:
	print("No nexus tree specified. Use -i/--input flag to specify input")
	quit()

tree = args.input

if args.format=='newick':
	with open(tree, "r") as input_handle:
		with open(tree.split('.')[0] + "_nw.tre", "w") as output_handle:
			sequences = Phylo.parse(input_handle, 'nexus')
			count = Phylo.NewickIO.write(sequences, output_handle)
	print("Converted %i records" % count)


if args.format=='nexus':
	with open(tree, "r") as input_handle:
		with open(tree.split('.')[0] + "_nx.tre", "w") as output_handle:
			sequences = Phylo.parse(input_handle, 'newick')
			count = Phylo.NexusIO.write(sequences, output_handle)
	print("Converted %i records" % count)

if args.format != 'nexus' and args.format != 'newick':
	print('format name is not correct use -f newick of -f nexus')
