#!/bin/bash
(cd exercises01 && exec ./run.sh)
(cd exercises02 && exec ./run.sh)
(cd exercises03 && exec ./run.sh)
echo -e "\n... Running ex5.1 ..."
(cd exercises05 && exec ./hydro>/dev/null)

echo -e "\nExercise 4 not run due to being computationally intense."
echo "If you wish to run it, go into the folder, read the README and execute manually."
