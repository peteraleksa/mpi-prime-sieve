//
//  main.c
//
//  Created by Peter on 3/14/14.
//
//

#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "blockdecomp.h"

#define MIN(a,b) ((a)<(b)?(a):(b))

typedef int bool;
#define true 1
#define false 0

int main(int argc, char * argv[])
{
    /* Constant Declarations */
    //long const 	SET_SIZE = 7920;

    /* Variable Declarations */
    int		count = 0;				// local count
    double 	elapsed_time = 0.00;			// time elapsed
    int		first;					// index of first multiple
    int 	global_count = 1;			// global count
    int 	high_value;				// highest value on processor
    char 	hostname[MPI_MAX_PROCESSOR_NAME];	// host process is running on
    int	 	i;					// counter variable
    int 	id;					// process id number
    int		index;
    int 	init_status;			// initialization error status flag
    bool  	initialized = false;		// mpi initialized flag
    int 	len;				// hostname length
    int 	low_value;			// lowest value on the processor
    char*	marked;				// portion of 2 to n that is marked
    int		n;			// number of elements to sieve
    int		n_sqrt;			// square root of n
    int 	p;			// number of processes
    int		prime;
    int		proc0_size;		// size of process 0's subarray
    int		size;			// elements in marked
    int*	sqrt_primes;		// primes up to the square root
    int		sqrt_primes_index;	// index in the square root primes array
    char*	sqrt_primes_marked;	// numbers up to sqrt marked prime or not
    int		sqrt_primes_size;	// size of square root primes array

    /* Function Declarations */
    //int is_prime( int );

    /* Initialization */
    MPI_Initialized( &initialized );                     // set initialized flag
    if( !initialized )                                  // if MPI is not initialized
        init_status = MPI_Init( &argc, &argv );        // Initialize MPI
    else
        init_status = MPI_SUCCESS;   	               // otherwise set init_status to success
    if( init_status != MPI_SUCCESS ) {     	       // if not successfully initialized
        printf ("Error starting MPI program. Terminating.\n");      // print error message
        fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, init_status);                     // abort
    }
    MPI_Get_processor_name( hostname, &len );                       // set hostname

    MPI_Comm_rank( MPI_COMM_WORLD, &id );                           // set process rank
    MPI_Comm_size( MPI_COMM_WORLD, &p );                            // set size of comm group
    //printf("Process rank %d started on %s.\n", id, hostname);     // print start message
    //fflush(stdout);
    //MPI_Barrier(MPI_COMM_WORLD );

    /* Start Timer */
    MPI_Barrier( MPI_COMM_WORLD );                                  // synchronize
    elapsed_time = - MPI_Wtime();                                   // start time

    /* Check that a set size was passed into the program */
    if(argc != 2) {
        if(id==0) {
            printf("Command line: %s <m>\n", argv[0]);
            fflush(stdout);
	}
        MPI_Finalize();
        exit(1);
    }

    n = atoi(argv[1]);
    n_sqrt = ceil(sqrt((double)n));
    //if(id==0)
    //	printf("square root: %i\n", n_sqrt);
    // debug
    //if(id==0) {
	//printf("n sqrt: %i\n", n_sqrt);
   	//fflush(stdout);
    //}

    sqrt_primes_marked = (char *) malloc(n_sqrt + 1);
    sqrt_primes_marked[0] = 1;
    sqrt_primes_marked[1] =1;

    for(i = 2; i <= n_sqrt; ++i) {
	sqrt_primes_marked[i] = 0;
    }

    prime = 2;
    sqrt_primes_size = n_sqrt;
    //printf("sqrt primes size: %i\n", sqrt_primes_size);

    do {
	for(i = prime * prime; i < n_sqrt; i+=prime) {
	     sqrt_primes_marked[i] = 1;
	     //sqrt_primes_size--;
	}
	while(sqrt_primes_marked[++prime]);    
    } while (prime * prime <= n_sqrt);
    //printf("sqrt primes size: %i\n", sqrt_primes_size);
    sqrt_primes = (int *) malloc(sqrt_primes_size);
    sqrt_primes_index = 0;

    //sqrt_primes_size = 0;

    for(i = 3; i <= n_sqrt; ++i) {
	if(!sqrt_primes_marked[i]) {
	    
	    sqrt_primes[sqrt_primes_index] = i;
	   // printf("%i, ", sqrt_primes[sqrt_primes_index]);
	    sqrt_primes_index++;
                
        }
    }

    sqrt_primes_size = sqrt_primes_index;

    //printf("sqrt primes size: %i\n", sqrt_primes_size);
    //fflush(stdout);

    /* Set process's array share and first and last elements */
    low_value = 2 + BLOCK_LOW(id,p,n-1);
    high_value = 2 + BLOCK_HIGH(id,p,n-1);
    size = BLOCK_SIZE(id,p,n-1);

    //printf("Process %i block low: %i\n", id, low_value);
    //fflush(stdout);
    //printf("Process %i block high: %i\n", id, high_value);
    //fflush(stdout);
    //printf("Block size: %i\n", size);
    //fflush(stdout);

    if(low_value % 2 == 0) {
	if(high_value % 2 == 0) {
	     size = (int)floor((double)size / 2.0);
	     high_value--;
	}
	else {
	    size = size / 2;
	}
	low_value++;
    }
    else {
	if(high_value % 2 == 0) {
	     size = size / 2;
	     high_value--;
	}
	else {
	     size = (int)ceil((double)size / 2.0);
	}
    }

    //printf("Process %i block low: %i\n", id, low_value);
    //fflush(stdout);
    //printf("Process %i block high: %i\n", id, high_value);
    //fflush(stdout);
    //printf("Block size: %i\n", size);
    //fflush(stdout);

    //proc0_size = (n-1)/p;

    /* if process 0 doesn't have all the primes for sieving, then bail*/
    /*if((2+proc0_size) < (int)sqrt((double)n)) {
        if(id==0) {
            printf("Too many processes\n");
            fflush(stdout);
        }
        MPI_Finalize();
        exit(1);
    }
    */

    /* Allocate share of array */
    marked = (char *) malloc(size);

    if(marked == NULL) {
        printf("Cannot allocate enough memory\n");
        fflush(stdout);
	MPI_Finalize();
        exit(1);
    }

    /* Run Sieve */

    //printf("made it to sieve\n");
    //fflush(stdout);

    for(i = 0; i < size; i++)
	marked[i] = 0;

    if(id==0)
	first = 0;
    
    sqrt_primes_index = 0;
    prime = sqrt_primes[sqrt_primes_index];

    //printf("first prime: %i\n", prime);
    //fflush(stdout);

    //for(i = 0; i < sqrt_primes_size; i++) {

      //              printf("%i,", sqrt_primes[i]);
        //            fflush(stdout);

        //}
     

    do {
	if(prime >= low_value)
	    first = ((prime - low_value) / 2) + prime;
	else if(prime * prime > low_value) {
		first = (prime * prime - low_value) / 2;
	}
	else {
	    if(low_value % prime == 0)
		first = 0;
	    else {
		first = 1;
		while ((low_value + (2 * first)) % prime != 0)
			++first;
	    }
	}

	//printf("first: %i\n", first);
	//fflush(stdout);

	for(i = first; i < size; i += (prime))
		marked[i] = 1;

	//printf("made it to prime assignment\n");
	prime = sqrt_primes[++sqrt_primes_index];

	//printf("prime: %i\n", prime);
	//fflush(stdout);

    } while(prime * prime <= n && sqrt_primes_index < sqrt_primes_size);

    count = 0;

    for(i = 0; i < size; i++) {
	if(!marked[i])
	    count++;
    }

    //printf("size: %i\ncount: %i\n", size, count);

//    for( i=id; i<SET_SIZE; i+=p )                                                       // interleaved allocation
//        count += is_prime( i );                                                             // check if prime w/ sieve of eratosthenes

    /* Reduce Sum */
    MPI_Reduce( &count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );        // reduce the primes count, root: proces 0

    /* Stop Timer */
    elapsed_time += MPI_Wtime();                                                        // end time

    //printf("Process %i found %i primes.\n", id, count);
    //fflush(stdout);

    //printf("Process %d is done in %d, running on %s.\n", id, elapsed_time, hostname);   // print process done message
    if( id == 0 ) {                                                                     // rank 0 prints global count
        printf("There are %d primes in the first %i integers.\nExecution took %10.6f.\n",
               global_count, n, elapsed_time);
	fflush(stdout);
	
//	printf("Debug:\n");
//	fflush(stdout);
//	printf("sqrt primes size: %i\n", sqrt_primes_size);
//        fflush(stdout);
	for(i = 0; i < sqrt_primes_size; i++) {
		if(!sqrt_primes[i]){
		    printf("%i,", sqrt_primes[i]);
		    fflush(stdout);
		}
	}
    }

    MPI_Barrier(MPI_COMM_WORLD);

  //  printf("rank: %i\nlow value: %i\nhigh value: %i\ncount: %i\n", id, low_value, high_value, count);

    //fflush(stdout);
    MPI_Finalize();                                                                     // finalize
    return 0;
}


/**
 *
 *  function is_prime
 *
 *      input: int x  -  number to check
 *
 *      process: mark prime with sieve of eratosthenes
 *
 *      output: 1 for prime, 0 for not prime
 *


int is_prime( int x )
{
    int j;
    if( (x == 1) || (x % 2 == 0 && x != 2) ) {                                      // if x is 1 or an even number other than 2
        return 0;                                                                       // then its not prime
    }
    else {
        for( j=3; j*j <= x; j+=2 ) {                                                // else for j as odd ints 3 to the square root of x
            if( (x % j) == 0 )                                                          // if x is divisible by any j
                return 0;                                                                   // then its not prime
        }
        exit(1);                                                                   // otherwise its prime
    }
}

*/
