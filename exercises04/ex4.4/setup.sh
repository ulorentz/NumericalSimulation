#!/bin/bash
cd ../code
make
cd ../ex4.4

for d in */; do cp ../code/MolDyn "$d"; done
for d in */; do cp ../code/config.0 "$d"; done
