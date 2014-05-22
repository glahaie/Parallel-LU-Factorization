Parallel LU Factorisation using MPI
===========

This project implements different parallel message-passing approaches for Matrix
LU Factorisation. The goal is to compare the execution time for different
agglomeration schemes and the cost of partial pivoting. We based our work
on Michael T. Heath's class 
[Parallel Numerical Algorithms](http://courses.engr.illinois.edu/cs554/fa2013/notes/)

We implemented 3 agglomeration schemes and a sequential version:

 - Agglomeration by adjacent lines
 - Agglomeration by cyclic lines
 - Agglomeration by 2D adjacent blocks

For the all versions except the 2D blocks, it is possible to specify if you want
to use partial pivoting. The default behavior of the applications is to generate
a random square matrix and calculate the factorization. The program does not check
if a factorization actually exists. The result is stocked in the same matrix, the 
lower-triangular matrix containing the L matrix and the main diagonal and upper 
triangular matrix containing the U matrix.

**Usage**

Arguments:

 - **-m** <matrix size>: **Required**. Size of the matrix. For the moment the
   applications only treat square matrixes, so only one size is required.
 - **-f** <file path>: Text file to specify the input matrix. If not specified
   the application will generate a random matrix.
 - **-s** <file path>: Solution file, to verify the result calculated by the
   application.
 - **-i**: print the result to screen.
 - **-w** <file path>: write the result to file.
 - **-p**: use partial pivoting. Not available for lu\_mpi\_cart
 - **-x, y** <number>: **Required for lu\_mpi\_cart**, specify the size of a block,
   the sizes must be divisors of the number of nodes used and must 

The makefile contains basic tests and a command to compare times for all applications
given a size of a matrix and a number of nodes to use.

This work was done for INF7235 at UQAM. Our original report was written in French and
is available in the latex folder.
