#ifndef _LELA_LELA_CONFIG_H
#define _LELA_LELA_CONFIG_H 1
 
/* lela/lela-config.h. Generated automatically at end of configure. */
/* lela/lela-config.h.  Generated from lela-config.h.in by configure.  */
/* lela/lela-config.h.in.  Generated from configure.in by autoheader.  */

/* Define if building universal (internal helper macro) */
/* #undef AC_APPLE_UNIVERSAL_BUILD */

/* Define if BLAS routines are available */
/* #undef BLAS_AVAILABLE */

/* Define if GMP is version 3.xxx */
/* #undef GMP_VERSION_3 */

/* Define that architecture uses big endian storage */
/* #undef HAVE_BIG_ENDIAN */

/* Define if BLAS is installed */
/* #undef HAVE_BLAS */

/* Define if C interface to BLAS is available */
/* #undef HAVE_CBLAS */

/* Define to 1 if you have the <dlfcn.h> header file. */
#ifndef __LELA_HAVE_DLFCN_H
#define __LELA_HAVE_DLFCN_H 1
#endif

/* Define if GMP is installed */
#ifndef __LELA_HAVE_GMP
#define __LELA_HAVE_GMP 1
#endif

/* Define to 1 if you have the <inttypes.h> header file. */
#ifndef __LELA_HAVE_INTTYPES_H
#define __LELA_HAVE_INTTYPES_H 1
#endif

/* Enable use of libpng */
#ifndef __LELA_HAVE_LIBPNG
#define __LELA_HAVE_LIBPNG 1
#endif

/* Define that architecture uses little endian storage */
#ifndef __LELA_HAVE_LITTLE_ENDIAN
#define __LELA_HAVE_LITTLE_ENDIAN 1
#endif

/* Define if M4RI is installed */
/* #undef HAVE_M4RI */

/* Version of M4RI is at least 20110601 */
/* #undef HAVE_M4RI_GE_20110601 */

/* Define to 1 if you have the <memory.h> header file. */
#ifndef __LELA_HAVE_MEMORY_H
#define __LELA_HAVE_MEMORY_H 1
#endif

/* Define to 1 if you have the <stdint.h> header file. */
#ifndef __LELA_HAVE_STDINT_H
#define __LELA_HAVE_STDINT_H 1
#endif

/* Define to 1 if you have the <stdlib.h> header file. */
#ifndef __LELA_HAVE_STDLIB_H
#define __LELA_HAVE_STDLIB_H 1
#endif

/* Define to 1 if you have the <strings.h> header file. */
#ifndef __LELA_HAVE_STRINGS_H
#define __LELA_HAVE_STRINGS_H 1
#endif

/* Define to 1 if you have the <string.h> header file. */
#ifndef __LELA_HAVE_STRING_H
#define __LELA_HAVE_STRING_H 1
#endif

/* Define to 1 if you have the <sys/stat.h> header file. */
#ifndef __LELA_HAVE_SYS_STAT_H
#define __LELA_HAVE_SYS_STAT_H 1
#endif

/* Define to 1 if you have the <sys/types.h> header file. */
#ifndef __LELA_HAVE_SYS_TYPES_H
#define __LELA_HAVE_SYS_TYPES_H 1
#endif

/* Define to 1 if you have the <unistd.h> header file. */
#ifndef __LELA_HAVE_UNISTD_H
#define __LELA_HAVE_UNISTD_H 1
#endif

/* Canonical 16-bit data type */
#ifndef __LELA_INT16
#define __LELA_INT16 short
#endif

/* Canonical 32-bit data type */
#ifndef __LELA_INT32
#define __LELA_INT32 int
#endif

/* Canonical 64-bit data type */
#ifndef __LELA_INT64
#define __LELA_INT64 long
#endif

/* Canonical 8-bit data type */
#ifndef __LELA_INT8
#define __LELA_INT8 char
#endif

/* Define to the sub-directory in which libtool stores uninstalled libraries.
   */
#ifndef __LELA_LT_OBJDIR
#define __LELA_LT_OBJDIR ".libs/"
#endif

/* Name of package */
#ifndef __LELA_PACKAGE
#define __LELA_PACKAGE "lela"
#endif

/* Define to the address where bug reports for this package should be sent. */
#ifndef __LELA_PACKAGE_BUGREPORT
#define __LELA_PACKAGE_BUGREPORT "lela-users@googlegroups.com"
#endif

/* Define to the full name of this package. */
#ifndef __LELA_PACKAGE_NAME
#define __LELA_PACKAGE_NAME "lela"
#endif

/* Define to the full name and version of this package. */
#ifndef __LELA_PACKAGE_STRING
#define __LELA_PACKAGE_STRING "lela 0.1.0"
#endif

/* Define to the one symbol short name of this package. */
#ifndef __LELA_PACKAGE_TARNAME
#define __LELA_PACKAGE_TARNAME "lela"
#endif

/* Define to the home page for this package. */
#ifndef __LELA_PACKAGE_URL
#define __LELA_PACKAGE_URL ""
#endif

/* Define to the version of this package. */
#ifndef __LELA_PACKAGE_VERSION
#define __LELA_PACKAGE_VERSION "0.1.0"
#endif

/* The size of `char', as computed by sizeof. */
#ifndef __LELA_SIZEOF_CHAR
#define __LELA_SIZEOF_CHAR 1
#endif

/* The size of `int', as computed by sizeof. */
#ifndef __LELA_SIZEOF_INT
#define __LELA_SIZEOF_INT 4
#endif

/* The size of `long', as computed by sizeof. */
#ifndef __LELA_SIZEOF_LONG
#define __LELA_SIZEOF_LONG 8
#endif

/* The size of `long long', as computed by sizeof. */
#ifndef __LELA_SIZEOF_LONG_LONG
#define __LELA_SIZEOF_LONG_LONG 8
#endif

/* The size of `short', as computed by sizeof. */
#ifndef __LELA_SIZEOF_SHORT
#define __LELA_SIZEOF_SHORT 2
#endif

/* The size of `__int64', as computed by sizeof. */
#ifndef __LELA_SIZEOF___INT64
#define __LELA_SIZEOF___INT64 0
#endif

/* The size of `__uint128_t', as computed by sizeof. */
#ifndef __LELA_SIZEOF___UINT128_T
#define __LELA_SIZEOF___UINT128_T 16
#endif

/* The size of `__uint256_t', as computed by sizeof. */
#ifndef __LELA_SIZEOF___UINT256_T
#define __LELA_SIZEOF___UINT256_T 0
#endif

/* Define to 1 if you have the ANSI C header files. */
#ifndef __LELA_STDC_HEADERS
#define __LELA_STDC_HEADERS 1
#endif

/* Canonical 128-bit data type */
#ifndef __LELA_UINT128
#define __LELA_UINT128 __uint128_t
#endif

/* Canonical 256-bit data type */
/* #undef UINT256 */

/* Version number of package */
#ifndef __LELA_VERSION
#define __LELA_VERSION "0.1.0"
#endif

/* Define WORDS_BIGENDIAN to 1 if your processor stores words with the most
   significant byte first (like Motorola and SPARC, unlike Intel). */
#if defined AC_APPLE_UNIVERSAL_BUILD
# if defined __BIG_ENDIAN__
#  define WORDS_BIGENDIAN 1
# endif
#else
# ifndef WORDS_BIGENDIAN
/* #  undef WORDS_BIGENDIAN */
# endif
#endif
 
/* once: _LELA_LELA_CONFIG_H */
#endif
