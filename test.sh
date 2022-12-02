#!/usr/bin/env bash

export LD_LIBRARY_PATH=$CLIB

./cnautilus \
-seqin data/1hr2/1hr2_final.fasta \
-mtzin data/1hr2/1hr2_final.mtz \
-colin-fo FP,SIGFP \
-colin-fc FWT,PHWT \
-colin-free FREE \
-cycles 3 \
-anisotropy-correction \
-pdbout data/1hr2/nautilus_output.pdb
