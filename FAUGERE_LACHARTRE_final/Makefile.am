INCLUDES=-I$(top_srcdir) -I$(top_builddir)
AM_CXXFLAGS = -std=c++0x -g -Wall -O3 -DDEBUG -Wno-sign-compare -DSHOW_PROGRESS_ -DDETAILED_PROFILE_TIMERS_ -DDEFAULT_BLOC_SIZE=256 -funroll-loops -finline-functions -msse2 -fopenmp

BENCHMARK_CXXFLAGS = -O2

SUBDIRS = 

AM_CPPFLAGS= $(GMP_CFLAGS) $(PNG_CFLAGS) $(M4RI_CFLAGS)
LDADD = $(GMP_LIBS)  $(PNG_LIBS) $(M4RI_LIBS) $(BLAS_LIBS) $(top_builddir)/lela/liblela.la -lstdc++ -lm -lpthread -lrt

AM_LDFLAGS = -fopenmp

noinst_PROGRAMS = \
	test-FGL-seq	\
	test-FGL-parallel

BASIC_TESTS =			\
		test-FGL-seq	\
		test-FGL-parallel	\
		test-structured-gauss	\
		test-LELA-FGL-Blas		
		
		
EXTRA_PROGRAMS = $(NON_COMPILING_TESTS) $(BENCHMARKS)

TESTS =				\
		$(BASIC_TESTS)

check_PROGRAMS = $(TESTS)

CLEANFILES = $(TESTS)

test_FGL_seq_SOURCES =						\
		test-FGL-seq.C						\
		../util/support.C
		
		
test_FGL_parallel_SOURCES =						\
		test-FGL-parallel.C						\
		../util/support.C										

		
test_structured_gauss_SOURCES =				\
		test-structured-gauss.C					\
		../util/support.C
		
test_LELA_FGL_Blas_SOURCES =				\
		test-LELA-FGL-Blas.C				\
		../util/support.C
		
		
noinst_HEADERS =	\
		../util/support.h		
