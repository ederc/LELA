INCLUDES=-I$(top_srcdir) -I$(top_builddir)
AM_CXXFLAGS = -g -Wall -O3 -DDEBUG -Wno-sign-compare -DSHOW_PROGRESS_ -DDETAILED_PROFILE_TIMERS_    \
	 -funroll-loops -finline-functions	-msse2

BENCHMARK_CXXFLAGS = -O2

SUBDIRS = 

AM_CPPFLAGS= $(GMP_CFLAGS) $(PNG_CFLAGS) $(M4RI_CFLAGS)
LDADD = $(GMP_LIBS)  $(PNG_LIBS) $(M4RI_LIBS) $(BLAS_LIBS) $(top_builddir)/lela/liblela.la -lstdc++ -lm -lpthread

BASIC_TESTS =			\
		test-FG-multiline

EXTRA_PROGRAMS = $(NON_COMPILING_TESTS) $(BENCHMARKS)

TESTS =				\
		$(BASIC_TESTS)

check_PROGRAMS = $(TESTS)

CLEANFILES = $(TESTS)

test_FG_multiline_SOURCES =						\
		test-FG-multiline.C						\
		multiline-indexer.h						\
		matrix-util-m.C							\
		matrix-op-m.C							\
		../../../util/support.C
	
noinst_HEADERS =	\
	../../../util/support.h