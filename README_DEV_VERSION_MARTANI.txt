This file contains a brief description of the different development files developed during the 
internship of MARTANI Fayssal at the polsys team in LIP6.

For information about how to build the LELA library, refer to the README file.

Code organization
=================

The different versions of the code are organized in directories as following:

new-faugere-lachartre: contains the very first implementation. This implementation is generic and 
			works over any field. Computes only the REDUCED row echelon form. [COMPLETE]
			
only-D:	This version is an amelioration of the 'new-faugere-lachartre' version. It computed echelon
			forms in addition to reduced row echelon forms. [COMPLETE]
			
FG-multiline: In this version, the multiline data structure is used allowing more performance 
			improvements. An implementation of the new method (the new ordering of operations) is 
			also provided starting from this version. [COMPLETE]
			
FG-multiline-generic: Contains a generic implementation of the multiline data structure. [INCOMPLETE] 	

FG-bloc*: In these directories, the block decomposition of matrices is used. In each directory, 
			a different data structure is used to represent the blocks. (the code is not necessarily
			 complete, only the performance tests are developed in these versions.) [INCOMPLETE]
			
final_FAUGERE_LACHARTRE: The final functional version using blocks and the multiline data structure.
			Both the standard Faugere-Lachartre and the new methods are implemented along with their
			parallel versions. 
			

Program options:
================
To get the list of options used in the program, after a make, invoke the executable with '-?', 
for example for the sequential program in the final version:
		final_FAUGERE_LACHARTRE/test-FGL-seq -?
		
		
The programs (sequential and parallel) are invoked as following:
		./test-FGL-[seq|parallel] - -f <MATRIX_FILE> <params>

The matrix file format, the multiline data structure and the block data structures are explained in 
the internship report [1]. Please contact MARTANI Fayssal 
martani.net@gmail.com if you need more information.

The main options are:
**Required**
-f FILE: The path to the file where the matrix is stored.

**optional**
-: (only a dash "-") shows progress and different information to the use (the rank, density etc.)
-r: to compute the REDUCED row echelon form. If this flag is not specified, only a row echelon form
	is computed.
-o: to use the standard Faugere-Lachartre method. If this flag is not present, the new method 
	proposed in [1] is used.
-s: validate results against structured Gaussian elimination method. By default, the results are 
	not validated. This could take very long time if used on non small matrices.
	

-p NUM_THREADS : The number of threads to used in case of the parallel version, default to 8.


example
-------

-Compute the REDUCED row echelon form using the standard Faugere-Lachartre
		./test-FGL-[seq|parallel] - -f <MATRIX_FILE> -r -o
		
Compute the a row echelon form (not reduced) using the new method and validating results
		./test-FGL-[seq|parallel] - -f <MATRIX_FILE> -s
		
Compute a REDUCED row echelon form using the new method, no result validation, parallel with 32 
threads
		./test-FGL-[seq|parallel] - -f <MATRIX_FILE> -r -p 32


Building test code
==================

1. Build LELA
-------------
in the LELA directory: ./autogen; ./configure; make

2.Build the new Faugere-Lachartre code:

In the directories that we have mentioned above, there is always a test file starting with "test-*"
to test the code in that directory.

For the final version in the directory 'final_FAUGERE_LACHARTRE', to make the sequential version:
	cd final_FAUGERE_LACHARTRE
	make test-FGL-seq
	
for the parallel version:
	cd final_FAUGERE_LACHARTRE
	make test-FGL-parallel
	
Note: clang doesn't support OpenMp, gcc or another OpenMP compliant compiler must be used.



Compiling options:
==================
Several options are available as compile time entries. They can be changed in the Makefile.am.
Every time Makefile.am is changed, ./autogen must be invoked so that the actual Makefiles are 
generated.

To change the size of the blocks to 128x128 for example, in final_FAUGERE_LACHARTRE/Makefile.am, 
edit -DDEFAULT_BLOC_SIZE like this:
		-DDEFAULT_BLOC_SIZE=128 

To enable time profiling (warning: this slows greately the code) enable:
		-DDETAILED_PROFILE_TIMERS
		
To show progress metrics enable (Do not used when saving output to files):
		-DSHOW_PROGRESS




[1] "Faugère-Lachartre Parallel Gaussian Elimination for Gröbner Bases Computations Over Finite 
Fields: Implementation and New Algorithms" MARTANI Fayssal 2012. 