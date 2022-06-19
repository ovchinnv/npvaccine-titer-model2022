
## Description

 This repository contains the code and use cases to reproduce the data
shown in the figures of the paper **Phenomenological modeling of antibody
response from vaccine strain composition**, 2022 by *V. Ovchinnikov and
M. Karplus* (submitted).

## Code Summary

The modeling of antibody titers performed here is implemented as a
collection of scripts for [Matlab](http://www.mathworks.com). The scripts
have been tested using Matlab Release 2016a with the Bioinformatics
Toolbox on ArchLinux X86_84 with kernel version 5.10.61-1-lts; however,
other versions of Matlab should work as well. Furthermore, the some
scripts will also work with [GNU
Octave](http://www.gnu.org/software/octave/index). However, those that
involve sequence alignments will not work with Octave.


## Required software:
 Name | Publisher | Purpose
----------|-------|---------
 *MATLAB* | The Mathworks, Inc. | Simulation & Analysis
 *mex-sqlite3 extension* | https://github.com/rmartinjak/mex-sqlite3 | Regenerate flu HA alignment from sqlite3 database (optional)


## How to Run

The main code for the antibody titer model resides in the root directory
of the repository. The files are:


File Name | Description
----------|------------
*aln2coor.m | transform a multiple sequence alignment (MSA) to a matrix with numerical values using the Grantham, 1974 encoding or the Atchley 2005 encoding
*ntaa.m* | functions to convert between 1- and 3-letter aa and nt code representations
*pdbout.m* | helper function to write a subset of coordinates, occupancy and temperature factors to a file in PDB format
*flu/dist2ave.m* | distance-to-average strain model (#1) for the flu hemagglutinin (HA)
*flu/avedist.m* | average distance to strain model (#2) for the flu hemagglutinin (HA)
*flu/3lzg.pdb* | a processed PDB file based on structure 3LZG for visualizing weights
*flu/mkpdbwgt.m* | transfer residue weights to beta column of structure 3LZG for visualization, as in Fig. 3D.
*flu/flu_strains.db* | our curated sqlite3 database of influenza HA strains compiled from several sources (GISAID, fludb, HIN); up to Oct. 2019.
*flu/strains.m* | misc. script to define strain name abbreviations for HA
*flu/mkmsa.m* | use flu_strains.db to create a MSA and store in file msaheadstem.mat
*flu/mkmsa* | example of an executable bash script that calls Matlab to generate the MSA
*flu/mkcoor.m* | load MSA for flu HA and convert to numerical vectors using either Grantham or Atchley encoding
*flu/mkvac.m* | define vaccine composition using strain codes in 'strains' and read experimental titer data of Cohen et al 2021
*flu/getind.m* | misc. function to return the array index corresponding to a flu virus strain abbreviation
*flu/showfit.m* | print comparison of model and experimental titers, e.g. Figs. 3A,B in the paper.
*flu/mkfig3ab* | executable bash script that calls Matlab to generate paper Fig. 3A,B
*flu/mkfig4ab* | executable bash script that calls Matlab to generate paper Fig. 4A,B
*flu/mkfigS1ab* | executable bash script that calls Matlab to generate paper Fig. S1A,B
*flu/mkfigS2ab* | executable bash script that calls Matlab to generate paper Fig. S2A,B
*flu/mkfigS3ab* | executable bash script that calls Matlab to generate paper Fig. S3A,B
*flu/train-test/allvacs.m* | compute correlation coefficients for all possible train/test experimental titer data splits
*flu/train-test/test.m* | compute test-set correlation coefficients for a single train/test data split
*flu/train-test/show.m* | generate correlation coefficient scatter plot
*flu/train-test/mkvac.m* | generate a vaccination cocktail based on the current permutation in allvacs.m
*flu/train-test/mkfig3c* | executable bash script that calls Matlab to generate paper Fig. 3C
*flu/train-test/mkfig4c* | executable bash script that calls Matlab to generate paper Fig. 4C
*flu/train-test/mkfigS1c* | executable bash script that calls Matlab to generate paper Fig. S1C
*flu/train-test/mkfigS2c* | executable bash script that calls Matlab to generate paper Fig. S2C
*flu/train-test/mkfigS3c* | executable bash script that calls Matlab to generate paper Fig. S3C
*flu/scan/mkscan.m* | script to perform the parameter scan shown in Fig. 2
*flu/scan/show.m* | generate a 2D contour plot of parameter scan (Fig. 2A)
*flu/scan/show1d.m* | generate a 1D contour plot of parameter scan (Fig. 2B)
*flu/scan/mkscan* | executable bash script to generate in Fig. 2
*cov/dist2ave.m* | distance-to-average strain model (#1) for the coronavirus (COV) RBD
*cov/avedist.m* | average distance to strain model (#2) for the COV-RBD
*cov/getind.m* | misc. function to return the array index corresponding to a COV strain abbreviation
*cov/mkcoor.m* | create MSA for COV and convert to numerical vectors using either Grantham or Atchley encoding
*cov/strains.m* | misc. script to define strain name abbreviations for COV
*cov/mkvac.m* | define vaccine composition using strain codes in 'strains' and read experimental titer data of Cohen et al 2021
*cov/seqwrite.m* | function to write a sequence of characters in FASTA or SELEX format
*cov/seqs/\*fasta* | sequence files for generating alignments
*cov/showfit.m* | print comparison of model and experimental titers, e.g. Figs. 5-6A,B in the paper.
*cov/mkmsa* | example of an executable bash script that calls Matlab to generate the alignment for COV-RBD
*cov/mkfig5ab* | executable bash script that calls Matlab to generate paper Fig. 5A,B
*cov/mkfig6ab* | executable bash script that calls Matlab to generate paper Fig. 6A,B
*cov/train-test/allvacs.m* | compute correlation coefficients for all possible train/test experimental titer data splits
*cov/train-test/test.m* | compute test-set correlation coefficients for a single train/test data split
*cov/train-test/show.m* | generate correlation coefficient scatter plot
*cov/train-test/mkvac.m* | generate a vaccination cocktail based on the current permutation in allvacs.m
*cov/train-test/mkfig5c* | executable bash script that calls Matlab to generate paper Fig. 5C
*cov/train-test/mkfig6c* | executable bash script that calls Matlab to generate paper Fig. 6C
*cov/6vxx.pdb* | a processed PDB file based on structure 6VXX for visualizing COV-RBD weights
*cov/mkpdbwgt.m* | transfer residue weights to beta column of structure 6VXX for visualization, as in Fig. 5D.

## NOTES

The directory structure of the scripts must be preserved, e.g. the script
files in flu/scan depend on the files in the parent directory (../), etc.
Note also that to regenerating the sequence alignment used for the flu
model requires the Matlab/Octave plugin mex-sqlite3. However, the
provided alignment file msaheadstem.mat can be used instead.

To generate the data for the figures, look for executable bash files
within the directory structure. For example, plots for figure Fig 3A,B
can be generated as follows :

 `cd flu`

 `./mkfig3ab`

Which runs the appropriate scripts in Matlab, in this case dist2ave.m
followed by showfit.m

