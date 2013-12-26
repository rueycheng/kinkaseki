#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.65)
AC_INIT([kinkaseki],1.0,[rueycheng@turing.csie.ntu.edu.tw])
AM_INIT_AUTOMAKE([kinkaseki],1.0)
AC_CONFIG_SRCDIR([src/tokenize.cpp])
AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

AC_CXX_HEADER_UNORDERED_MAP
AC_CXX_HEADER_TR1_UNORDERED_MAP
AM_CONDITIONAL([USE_CXX0X], [test $ac_cv_cxx_unordered_map = yes])

# Checks for library functions.
AC_CONFIG_FILES([Makefile
                 src/Makefile
		 script/Makefile])
AC_OUTPUT
