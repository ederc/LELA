INCLUDES=-I$(top_srcdir) -I$(top_builddir)
AM_CXXFLAGS = -g -Wall

BENCHMARK_CXXFLAGS = -O2

SUBDIRS = 

AM_CPPFLAGS= $(GMP_CFLAGS) $(PNG_CFLAGS) $(M4RI_CFLAGS)
LDADD = $(GMP_LIBS)  $(PNG_LIBS) $(M4RI_LIBS) $(BLAS_LIBS) $(top_builddir)/lela/liblela.la

BASIC_TESTS =			\
		test-FG-multiline

EXTRA_PROGRAMS = $(NON_COMPILING_TESTS) $(BENCHMARKS)

TESTS =				\
		$(BASIC_TESTS)

check_PROGRAMS = $(TESTS)

CLEANFILES = $(TESTS)

test_FG_multiline_SOURCES =						\
		test-FG-multiline.C						\
		../new-faugere-lachartre/matrix-util.h	\
		../new-faugere-lachartre/indexer.h		\
		multiline-indexer.h						\
		../util/support.C

	
noinst_HEADERS =	\
	../util/support.h