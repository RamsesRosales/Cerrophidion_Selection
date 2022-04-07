Scripts\_Description
================
Ramses Alejandro Rosales Garcia
4/7/2022

-   [R Scripts](#r-scripts)
    -   [RSEM Formating and Ploting
        Pipeline](#rsem-formating-and-ploting-pipeline)
    -   [Differential Expression and Heat Map
        Pipeline](#differential-expression-and-heat-map-pipeline)
    -   [Phylogenetics Formating](#phylogenetics-formating)
    -   [Selection Table Formating](#selection-table-formating)
    -   [Other Scripts](#other-scripts)
-   [Python Scripts](#python-scripts)
    -   [Phylogenetic formating](#phylogenetic-formating)

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

*make\_FD\_1.R *make\_FD\_2\_final.R

## Other Scripts

-   tox\_numbers.R
-   mixing\_df\_tx.R

# Python Scripts

Script on Python language

## Phylogenetic formating

Scripts used at different steps of the phylogenetic analysis for busco
genes and for PLA2s
