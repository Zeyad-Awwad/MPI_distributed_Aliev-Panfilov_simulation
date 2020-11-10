# Aliev-Panfilov Cardiac Simulation with MPI

The software in this repository uses MPI to perform the Aliev-Panfilov cardiac simulation on UCSD's Comet supercomputer. The code demonstrates strong scaling and exceeded 1 tflop/s in our largest tests (400 cores).

The code was collaboratively developed with my project partner, Mohammad Samragh, as part of a parallel computation course. It was built upon a substantial starter code shared by Prof Bryan Chin, with our work primarily contained in the "solve.cpp" file. 

Our task was to integrate MPI into an existing model to properly distribute a finite-difference calculation among an arbitrary rectangular mesh of processors. This required (approximately) even load distribution among cores and efficient ghost cell management among the nearest neighbors to minimize communication bottlenecks. 
