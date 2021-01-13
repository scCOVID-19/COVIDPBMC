# Cydar Analysis

This folder contains scripts that were used to perform the Cydar analysis on the protein expression data.
The scripts are in order of how they are run to process the data and take the background corrected ADT expression SCE object as input.

- 01_fastMNN.R performs batch correction on the background corrected protein counts
- 02_CreateHyper.R assigns cells into hyperspheres
- 03_Cydar_Figure.Rmd then test for differential abundance of cells in each sphere and produces the plots from the manuscript
