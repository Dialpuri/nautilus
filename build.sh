source /Applications/ccp4-8.0/bin/ccp4.setup-sh 
source ~/Development/privateer/emsdk/emsdk_env.sh

emcmake cmake .
emmake make -j

mv pnautilus.js webserver/wasm/pnautilus.js
mv pnautilus.wasm webserver/wasm/pnautilus.wasm
mv pnautilus.data webserver/public/pnautilus.data

