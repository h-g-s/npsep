#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.50)
AC_INIT(npsep, 1.0, haroldo.santos@gmail.com)
AC_CONFIG_SRCDIR(src/memory.c)
AM_INIT_AUTOMAKE
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_LIBTOOL
AM_PROG_CC_C_O

# Checks for libraries.
AC_CHECK_LIB([m], [fabs])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h limits.h stdlib.h string.h stdio.h math.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_ALLOCA
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([memset strstr])

PKG_CHECK_MODULES([cbc], [cbc])

AC_CONFIG_FILES([Makefile src/Makefile])

AC_OUTPUT


