#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
AC_INIT([epi-emmax], [v1.01-r0], [slowkoni@gmail.com])
AC_CONFIG_SRCDIR([emmax-kin.c])
AM_INIT_AUTOMAKE

AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.
AC_CHECK_LIB([pthread], [pthread_create])
AC_CHECK_LIB([m], [exp])
AC_CHECK_LIB([z], [gzopen])
AC_SEARCH_LIBS([dgetri_],[lapack lapacke])
AC_SEARCH_LIBS([cblas_dgemv],[cblas gslcblas])

# Checks for header files.
AC_CHECK_HEADERS([float.h limits.h stdint.h stdlib.h string.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor gettimeofday memset sqrt strdup strerror strtol strtoul])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
