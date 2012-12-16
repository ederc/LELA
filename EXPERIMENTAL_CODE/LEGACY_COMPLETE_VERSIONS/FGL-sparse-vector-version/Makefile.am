INCLUDES=-I$(top_srcdir) -I$(top_builddir)
AM_CXXFLAGS = -g -Wall -DDEBUG -O3 -funroll-loops -Wno-sign-compare

BENCHMARK_CXXFLAGS = -O2

SUBDIRS = 

AM_CPPFLAGS= $(GMP_CFLAGS) $(PNG_CFLAGS) $(M4RI_CFLAGS)
LDADD = $(GMP_LIBS)  $(PNG_LIBS) $(M4RI_LIBS) $(BLAS_LIBS) $(top_builddir)/lela/liblela.la -lstdc++ -lm

BASIC_TESTS =					\
		test-faugere-lachartre		\
		test-structured-gauss		\
		test-splicer-indexer

EXTRA_PROGRAMS = $(NON_COMPILING_TESTS) $(BENCHMARKS)

TESTS =                               \
        $(BASIC_TESTS)

check_PROGRAMS = $(TESTS)

CLEANFILES = $(TESTS)

test_faugere_lachartre_SOURCES =		\
		test-faugere-lachartre.C		\
		matrix-util.h					\
		matrix-op.h						\
		indexer.h						\
		structured-gauss-lib.h			\
		../../../util/support.C

		
test_structured_gauss_SOURCES =		\
	test-structured-gauss.C 		\
	structured-gauss-lib.h			\
	matrix-util.h					\
	../../../util/support.C


test_splicer_indexer_SOURCES =		\
	test-splicer-indexer.C			\
	matrix-util.h					\
	indexer.h						\
	LELA-FG-util.C					\
	../../../util/support.C
	
	
noinst_HEADERS =	\
	../../../util/support.h