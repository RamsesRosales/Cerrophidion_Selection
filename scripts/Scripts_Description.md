Scripts\_Description
================
Ramses Alejandro Rosales Garcia
2022 April 11

-   [R Scripts](#r-scripts)
    -   [RSEM Formating and Ploting
        Pipeline](#rsem-formating-and-ploting-pipeline)
    -   [Differential Expression and Heat Map
        Pipeline](#differential-expression-and-heat-map-pipeline)
    -   [Phylogenetics Formating](#phylogenetics-formating)
    -   [Selection Table Formating](#selection-table-formating)
    -   [Other Scripts](#other-scripts)
-   [Python Scripts](#python-scripts)
    -   [Manual Annotation pipeline](#manual-annotation-pipeline)
    -   [Phylogenetic formating](#phylogenetic-formating)
-   [Shell Scripts](#shell-scripts)

These directory contain several small scripts used in the paper.

# R Scripts

Script on R language

## RSEM Formating and Ploting Pipeline

The scripts are designed to be used in the following order.

-   data\_frame\_consense.R
-   consense.R
-   FancyPlots3.R
-   FancyPlotsAverage.R
-   FancyPlotAverageNS.R
-   TransConPlot\_Species.R
-   stacked\_toxin\_plot.R
-   PlottingFunctions1.R (sourced by FancyPlots scripts)

## Differential Expression and Heat Map Pipeline

The script are design to be run in the following order

-   consense\_expected\_count.R
-   ToxNames.R
-   Dif\_Exp\_Cerrophidion.R
-   Dif\_Exp\_Tab.R
-   Colors\_HM\_Cerrophidion.R (list of colors used by heatmap script)
-   Cerrophidion\_HeatMap\_AverageExpressionOrder.R
-   Cerrophidion\_HeatMap\_order\_difExp.R (specific for Cgodmani data)

## Phylogenetics Formating

Scripts used in different steps of the phylogenetic analysis mostly to
format files.

-   Rename\_Fasta\_Headers.R
-   Clip\_trees.R
-   Change\_tree\_names.R

## Selection Table Formating

Scripts used in other steps of the project, or not used.

-   make\_FD\_1.R
-   make\_FD\_2\_final.R

## Other Scripts

-   tox\_numbers.R
-   mixing\_df\_tx.R

# Python Scripts

Script on Python language

## Manual Annotation pipeline

The manual annotation pipeline scripts were writen by [Darin Rokyta](https://drokyta.com) and [Andrew Mason](https://eeob.osu.edu/people/mason.501)
email them or to the correspondence authors to ask for those scripts

- Blast_parse_v6.py
- AutoAnnotator.py
- AutoAnnotatorFormat.py
- NextAnnotate_v0.3.py
- XML_Barber_v3.py

## Phylogenetic formating

Scripts used at different steps of the phylogenetic analysis for busco
genes and for PLA2s

-   BuscoCleaner.py (modified from Rhett Rautsau)
-   Translate\_RARG.py
-   filter\_cov.py (Modified from Andrew Script)
-   nx2nw.py
-   IsoelectricPoint\_df.py
-   extract\_CDS.py
-   import\_gb.py

# Shell Scripts

Some of the shell script used, other shell scripts might be seen in the
first section directory markdown.

-   add\_shebang.sh
-   busted\_maketable
