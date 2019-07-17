#!/bin/bash
cd ../code
make
cd ../ex07.1-4

for d in */; do cp ../code/MolecularMC "$d"; done
for d in */; do cp ../code/config.0 "$d"; done
