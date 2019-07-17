#!/bin/bash
cd ../code_mol
make
cd ../ex07.4_mol

for d in */; do cp ../code_mol/MolDyn "$d"; done
for d in */; do cp ../code_mol/config.0 "$d"; done
