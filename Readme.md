Sequence divergence in venom genes within and between montane pitvipers
(*Viperidae* :*Crotalinae*  :*Cerrophidion*) species is driven by
mutation-drift equilibrium.
================
Ramses A. Rosales-Garcia, Rhett M. Rautsaw , Erich P. Hofmann, Christoph
I. Grunwald, Jason M. Jones, Hector Franz-Chavez, Ivan T.
Ahumada-Carrillo, Ricardo Ramirez-Chaparro, Miguel Angel De la
Torre-Loranca, Jason L. Strickland, Andrew J. Mason, Matthew L. Holding,
Miguel Borja, Gamaliel Castaneda-Gaytan, Darin R. Rokyta, Tristan D.
Schramer, N. Jade Mellor, Edward A. Myers, Christopher Parkinson
2023 May 04

- <a href="#assembly-expression-seq-evolution"
  id="toc-assembly-expression-seq-evolution">Assembly Expression Seq
  Evolution</a>
- <a href="#second_part" id="toc-second_part">second_part</a>
- <a href="#supplementary" id="toc-supplementary">Supplementary</a>
- <a href="#data" id="toc-data">Data</a>

<figure>
<img src="pictures/C_tzotzilorum1.jpg" style="width:500.0%"
alt="Cerrophidion tzotzilorum" />
<figcaption aria-hidden="true"><em>Cerrophidion
tzotzilorum</em></figcaption>
</figure>

This Github page contains the code and scripts used in the development
of the paper.

# Assembly Expression Seq Evolution

The
[assembly_expression_seq_evolution](assembly_expression_seq_evolution/Scripts_Paper.md)
contains the shell code used to de novo assembly the transcriptomes,
annotation of toxins and nontoxin gene, gene expression, expression
ploting, differential expression, phylogenetics, variant calling and
selection analisys. There are several
[scripts](scripts/Scripts_Description.md) used that can be found in the
[scripts](scripts) folder. Some scripts were written by a third party,
in those cases we mentioned the location of the original script or the
contact were the script might be requested.

# second_part

The
[sequence_evolution_statistics](sequence_evolution_statistics/Cgodm_selection_BMC.md)
contains the R scripts with the statistics analysis that compare
selection signals in toxins vs nontoxins, for Tajima’s D, Tajima’s D in
synonymous SNPS, Tajima’D in nonsynonymous SNPs, Fst, and BUSTED. The
necessary data to run the analysis might be found withing the folder.

# Supplementary

The
[supplementary](https://github.com/RamsesRosales/Cerrophidion_Selection/tree/main/supplementary)
directory contains the suplementary material submitted to the table.

- Additional_Figures.pdf
- Additional_Tables.xlsx
- Annotation_Scripts.html
- Selection_Scripts.html

# Data

The
[data](https://github.com/RamsesRosales/Cerrophidion_Selection/tree/main/data)
directory contains the [consensus
transcriptomes](https://github.com/RamsesRosales/Cerrophidion_Selection/tree/main/data/Transcriptomes),
the [RSEM
results](https://github.com/RamsesRosales/Cerrophidion_Selection/tree/main/data/Expression)
and other files generated and analised, during this work.
