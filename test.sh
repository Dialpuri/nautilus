#!/usr/bin/env bash

export LD_LIBRARY_PATH=$CLIB

./cnautilus \
-seqin test_data/1hr2_final.fasta \
-mtzin test_data/1hr2_final.mtz \
-colin-fo FP,SIGFP \
-colin-fc FWT,PHWT \
-colin-free FREE \
-cycles 3 \
-anisotropy-correction \
-pdbout test_data/nautilus.pdb
