#!/bin/bash
mpirun -np 4 --use-hwthread-cpus rsMPI 1 1 100000 squared L1

