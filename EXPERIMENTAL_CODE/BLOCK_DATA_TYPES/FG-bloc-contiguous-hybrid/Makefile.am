INCLUDES=-I$(top_srcdir) -I$(top_builddir)
AM_CXXFLAGS = -g -Wall -O3 -DDEBUG -DSHOW_PROGRESS -DDETAILED_PROFILE_TIMERS_ -DDEFAULT_BLOC_HEIGHT=256    \
	 -DDEFAULT_BLOC_WIDTH=256 -funroll-loops -finline-functions

BENCHMARK_CXXFLAGS = -O2

SUBDIRS = 

AM_CPPFLAGS= $(GMP_CFLAGS) $(PNG_CFLAGS) $(M4RI_CFLAGS)
LDADD = $(GMP_LIBS)  $(PNG_LIBS) $(M4RI_LIBS) $(BLAS_LIBS) $(top_builddir)/lela/liblela.la -lstdc++ -lm -lpthread

BASIC_TESTS =			\
		test-FG-bloc
		
		
EXTRA_PROGRAMS = $(NON_COMPILING_TESTS) $(BENCHMARKS)

TESTS =				\
		$(BASIC_TESTS)

check_PROGRAMS = $(TESTS)

CLEANFILES = $(TESTS)

test_FG_bloc_SOURCES =						\
		test-FG-bloc.C						\
		../../../util/support.C
		
		
		
noinst_HEADERS =	\
		../../../util/support.h