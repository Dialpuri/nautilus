#!/usr/bin/env bash

source "/Applications/ccp4-8.0/bin/ccp4.setup-sh"
export PATH=/usr/bin:$PATH
make -j
./fix_library.sh
export LD_LIBRARY_PATH=$CLIB:/home/jordan/dev/nautilus/lib

rm nautilus_output*

time ./cnautilus \
-seqin tests/test_1u9s/1u9s.seq \
-mtzin tests/test_1u9s/1u9s.mtz \
-colin-fo FP,SIGFP \
-colin-fc sfcalc.F_phi.F,sfcalc.F_phi.phi \
-colin-free FREE \
-cycles 1 \
-anisotropy-correction \
-pdbout tests/test_1u9s/nautilus_output.pdb \
-predicted-phos-map tests/test_1u9s/1u9s_phos.map


#csymmatch -pdbin-ref tests/test_structures/pdb1u9s.ent -pdbin debug/1u9s/final.pdb -pdbout debug/1u9s/final_symmatch.pdb -connectivity-radius 2 -origin-hand