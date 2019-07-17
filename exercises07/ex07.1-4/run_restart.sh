#!/bin/bash
(cd solid && exec ./MolecularMC --restart --istant)
(cd gas && exec ./MolecularMC --restart --istant)
(cd liquid && exec ./MolecularMC --restart --istant)
