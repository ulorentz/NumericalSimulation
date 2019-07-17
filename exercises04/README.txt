INSTRUCTIONS
The code for this exercise is in the folder 'code' and there could be inspected.
The actual executables are placed in folder ex4.4, there you can find three 
scripts:
    - setup.sh: compile all and copy executables in the required subfolders.
    - run_one_config.sh: run simulations for solid, liquid and gas starting from
      the fcc configuration.
    - run_restart.sh: run all the simulations restarting from previous results.


NOTE: all the simulations in ex4.1, 4.2,4.3 and 4.4 are compiled with '-Ofast'
which reduces execution time on a modern processor and g++ 8.x of about 50% 
compared to a standard '-O3'. 
If you don't like non standard compilation flags, or on different compiler
you get no benefits or different results, just change in "code" folder the 
'CXXFLAGS' in the makefile. (but note that, since the code makes a massive use 
of std::vector and STL, using no optimization at all may have a huge impact). 
