# Nautilus 

This repository contains the nucleic acid model building software Nautilus, and is the site where any further developments on the software will be held. The original Nautilus software was written by Kevin Cowtan and this work is part of a PhD studentship funded by the Biotechnology and Biological Sciences Research Council (BBSRC)  awarded to Jordan Dialpuri.  
### Development

#### Prerequisites

You must have 
 - Clipper
 - MMDB2

installed, which come included in the CCP4. To ensure they are in your path, you must source the appropriate script. To do this run: 

    source /opt/xtal/ccp4-X.X/bin/ccp4.setup-sh 
where X.X is your CCP4 version. 
 
#### Development

To compile this, just simply run: 

    make
and the executable 'cnautilus' should be created in the root directory of the project. 

### Running

To run the executable run: 

    ./cnautilus <PARAMS>

The available parameters are
 
	-mtzin <filename>               
    -seqin <filename>
    -pdbin <filename>
    -pdblistin <filename> 
    -pdblistdir <filepath> 
    -pdbout <filename>
    -predicted-phos-map <filename>
    -xmlout <filename>
    -colin-fo <colpath>
    -colin-hl <colpath> or -colin-phifom <colpath>
    -colin-fc <colpath>
    -colin-free <colpath>
    -cycles <number>
    -anisotropy-correction
    -fragments <number>
    -resolution <resolution/A>
    -pdbin-ref <filename>
    -cif           
    -search-step <float>  
    -verbose <number>
    
#### Descriptions
mtzin - the input reflections data from the experiment e.g. 1hr2.mtz

seqin - the sequence of the nucleic acid e.g. 1hr2.fasta

pdbin - the filepath of an input model 

pdblistin - the list of PDB files for use in the library

pdblistdir - the directory containing all PDB files in the pdblistin 

pdbout - the filepath of the output model e.g. 1hr2_output.pdb

predicted-phos-map - the filepath of a predicted phosphate map 

xmlout - 

colin-fo - column headings for F<sub>obs</sub> e.g. FP, SIGFP

coln-hl - column headings for Hendrickson-Lattman coefficients (ABCD) e.g. sfcalc.A, sfcalc.B, sfcalc.C, sfcalc.D

colin-fc - column headings for F<sub>calc</sub> e.g FWT, PHWT

colin-free - column heading for Free R Flag

cycles - number of cycles to run

anisotropy-correction - apply anisotropy correction

-fragments - number of hits to be found in the first step 

-resolution - the resolution of the data

-pdbin-ref  - single PDB file for use as the library 

-cif - output model in CIF format only

-search-step - the step used in the fragment fitting angle search

-verbose - level of logging output
