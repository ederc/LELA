New Faugère-Lachartre implementation
===================================

The code for the new implementation of the Faugère-Lachartre method can be found 
in the directories `FAUGERE_LACHARTRE_final` and `EXPERIMENTAL_CODE/LEGACY_COMPLETE_VERSIONS`.

`FAUGERE_LACHARTRE_final` contains the final version, which is implemented using blocking, 
it contains the parallel and the sequential versions. Both the standard and the new method
of the algorithm are implemented in this version.

Building the test programs
--------------------------

In the `LELA` directory type:
<pre>
$ ./autogen
$ ./configure
$ make
</pre>

The test programs are compiled automatically in the directory `FAUGERE_LACHARTRE_final`.
### Note
clang doesn't support OpenMP for now, gcc or another OpenMP compliant compiler must be used.

Invoking the test programs
--------------------------

This is an examples of how to invoke the test programs:

* `./test-FGL-seq - -f path/to/matrix/file` 		performs a sequential Gaussian elimination using the new method.
* `./test-FGL-seq - -f path/to/matrix/file -r`		performs a sequential *reduced* Gaussian elimination using the new method.
* `./test-FGL-seq - -f path/to/matrix/file	-o -s`	performs a sequential Gaussian elimination using the standard method with result validation.

* `./test-FGL-parallel - -f path/to/matrix/file`		performs a parallel (default 8 threads) Gaussian elimination using the new method.
* `./test-FGL-parallel - -f path/to/matrix/file -o -r -p 16`	performs a parallel reduced Gaussian elimination using the standard method with 16 threads.

List of options in the test programs
------------------------------------

To get the list of the options available in the test programs type:

    ./test-FGL-seq -?


Mainly the options can be mixed together to get the desired computation;
for example, using the options `-r -o -s` means that one desires to compute:
* a reduced echelon form (the option `-r`)
* using the standard method (the option `-o` for old method)
* the results are to be validated against a more stable algorithm (option `-s`)

Compile time options
--------------------

In all the versions (in `FAUGERE_LACHARTRE_final` and `EXPERIMENTAL_CODE/LEGACY_COMPLETE_VERSIONS`):
* To show progress information (very verbose) enable the `-DSHOW_PROGRESS` flag in the `Makefile`.
* To enable time profiling (warning: this slows greately the code): `-DDETAILED_PROFILE_TIMERS`.

To change the block size in the `FAUGERE_LACHARTRE_final` version:
* Edit `-DDEFAULT_BLOC_SIZE` like this: `-DDEFAULT_BLOC_SIZE=128` for 128x128 elements blocks.



Note on the state of the code & earlier versions
------------------------------------------------

The code is in its very early stages of development and has a significant amount of debug/experimentation code;
the block version (`FAUGERE_LACHARTRE_final`) is pretty complicated (and lacks documentation/comments).
For earlier (and simpler) versions, refer to the directory `EXPERIMENTAL_CODE/LEGACY_COMPLETE_VERSIONS` where you can find:

* `FGL-sparse-vector-version`: the earliest and simplest implementation of the **standard algorithm**, using only sparse vectors (no blocks).
* `FGL-sparse-vector-new-method`: the same as the above version, but implements the **new method** (new ordering of operations).
* `FG-multiline: In this version`, the data structure **multiline** is used for both the standard and the new Faugère-Lachartre.