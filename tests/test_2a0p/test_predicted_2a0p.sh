#!/usr/bin/env bash

source "/Applications/ccp4-8.0/bin/ccp4.setup-sh"
export PATH=/usr/bin:$PATH
make -j
./fix_library.sh
export LD_LIBRARY_PATH=$CLIB:/home/jordan/dev/nautilus/lib

rm nautilus_output*

time ./cnautilus \
-seqin tests/test_2a0p/2a0p.seq \
-mtzin tests/test_2a0p/2a0p.mtz \
-colin-fo FP,SIGFP \
-colin-fc sfcalc.F_phi.F,sfcalc.F_phi.phi \
-colin-free FREE \
-cycles 1 \
-anisotropy-correction \
-pdbout tests/test_2a0p/nautilus_output.pdb \
-predicted-phos-map tests/test_2a0p/2a0p_phos.map


#csymmatch -pdbin-ref tests/test_structures/pdb2a0p.ent -pdbin debug/2a0p/final.pdb -pdbout debug/2a0p/final_symmatch.pdb -connectivity-radius 2 -origin-hand