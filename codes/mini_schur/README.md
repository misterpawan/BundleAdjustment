Idea of the mini schur complement method is to approximate global Schur complement with local Schur complements
Main test file: test_nssolve*.m

For the moment, it is for navier stokes matrices. Change it for bundle adjustment matrix. 
Right now, the code uses independent set ordering to make domain decomposition. 
But since bundle adjustment matrix is already in DD form, this is not required. That is the blocks D, E, F, and G can be 
identified directly.

Moreover, you may use domain_decomposition.m from the ddm folder. It does the required partitioning suitable for bundle adjustment problem.
