# Nautilus Webserver

![Workflow](https://github.com/Dialpuri/nautilus/actions/workflows/main.yml/badge.svg)

This repository contains the nucleic acid model building software Nautilus, and is the site where any further developments on the software will be held. The original Nautilus software was written by Kevin Cowtan and this work is part of a PhD studentship funded by the Biotechnology and Biological Sciences Research Council (BBSRC)  awarded to Jordan Dialpuri.  
### Development

#### Installation

To build backend:
```
git clone https://github.com/emscripten-core/emsdk.git
cd emsdk
git pull
./emsdk install latest
./emsdk activate latest
cd ..
./get_sources
source emsdk/emsdk_env.sh
emcmake cmake .
emmake make -j
mv pnautilus.js webserver/wasm/pnautilus.js
mv pnautilus.wasm webserver/wasm/pnautilus.wasm
mv pnautilus.data webserver/public/pnautilus.data
```

To build frontend: 
```
cd webserver
wget https://github.com/moorhen-coot/Moorhen/releases/download/v0.5.1/moorhen-0.5.1.tgz
npm ci
cp -r node_modules/moorhen/baby-gru public/
npm run dev 
```
