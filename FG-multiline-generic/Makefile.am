INCLUDES=-I$(top_srcdir) -I$(top_builddir)
AM_CXXFLAGS = -g -Wall -O3 -DSHOW_PROGRESS -DDETAILED_PROFILE_TIMERS

BENCHMARK_CXXFLAGS = -O2

SUBDIRS = 

AM_CPPFLAGS= $(GMP_CFLAGS) $(PNG_CFLAGS) $(M4RI_CFLAGS)
LDADD = $(GMP_LIBS)  $(PNG_LIBS) $(M4RI_LIBS) $(BLAS_LIBS) $(top_builddir)/lela/liblela.la -lstdc++ -lm

BASIC_TESTS =			\
		test-FG-multiline	\
		test-FG-multiline4

EXTRA_PROGRAMS = $(NON_COMPILING_TESTS) $(BENCHMARKS)

TESTS =				\
		$(BASIC_TESTS)

check_PROGRAMS = $(TESTS)

CLEANFILES = $(TESTS)

test_FG_multiline_SOURCES =						\
		test-FG-multiline.C						\
		../new-faugere-lachartre/matrix-util.h	\
		../new-faugere-lachartre/matrix-op.h	\
		../new-faugere-lachartre/indexer.h		\
		multiline-indexer.h						\
		matrix-util.C	\
		../util/support.C


test_FG_multiline4_SOURCES =						\
		test-FG-multiline4.C						\
		../new-faugere-lachartre/matrix-util.h	\
		../new-faugere-lachartre/matrix-op.h	\
		../new-faugere-lachartre/indexer.h		\
		multiline-indexer.h						\
		matrix-util.C	\
		../util/support.C

	
noinst_HEADERS =	\
	../util/support.h