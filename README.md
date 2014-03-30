# Sieve of Eratosthenes in C using OpenMPI

sieve_block_decomp.c is derived from examples and ideas presented in Parallel Programming in C with MPI and OpenMP by Michael J. Quinn, McGraw Hill, 2003. It also draws from course work and lectures for Brooklyn College CISC 4335 Parallel and Distributed Computing.

## Compiling

make

or: 

mpicc -o sieve_block_decomp sieve_block_decomp.c

## To Run Tests

There is a bash script (runscript.sh) which can be used to run a series of trials on a given input size.  It can be used as follows:

./runscript.sh sieve_block_decomp [input_size] [max_processors] [num_test_repetitions]

For example, to run the program with an input size of 50,000,000 on a max of 32 processors 100 times for each worker pool size, run:

./runscript.sh sieve_block_decomp 50000000 32 100

### Authors
Peter Aleksa Jan 2014 -  
