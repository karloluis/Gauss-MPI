Gaps between consecutive primes
	

by 
Cassandra Schaening-Burgos
	
and Karlo Martinez Martos



parallel_gauss.cpp
Uses gaussian elimination to solve a dense system of linear equations 

NOTE: In the same directory as the program, there must be a file named "coefficients.txt" which contains the system of linear equations represented as coefficients separated by spaces and newlines.




Compile with:

	$ mpic++ -fopenmp parallel_gauss.cpp -o parallel_gauss 

Run as:

	$ mpirun -np <numprocs> parallel_gauss
        Where <numprocs> is the number of processor to be used.
	This must be done on a cluster or computer with multiple processors and
 support for MPI.


