#!/bin/bash
(cd solid && exec ./MolecularMC --new-sim --istant)
(cd gas && exec ./MolecularMC --new-sim --istant)
(cd liquid && exec ./MolecularMC --new-sim --istant)
