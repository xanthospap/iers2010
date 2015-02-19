#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([libiers10++], [1.0], [xanthos@mail.ntua.gr])
AM_INIT_AUTOMAKE

AC_LANG(C++)

AC_CONFIG_SRCDIR([src/gpt.cpp])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_RANLIB
AC_PROG_SED
AC_PROG_MKDIR_P
AC_CHECK_PROGS([DOXYGEN], [doxygen])
if test -z "$DOXYGEN"
then
    AC_MSG_WARN([Doxygen not found - continuing without Doxygen])
fi
AM_CONDITIONAL([HAVE_DOXYGEN],
               [test -n "$DOXYGEN"])AM_COND_IF([HAVE_DOXYGEN],
                                               [AC_CONFIG_FILES([doc/Doxyfile])])

# Checks for optional programs.

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE

# Checks for library functions.
AC_CHECK_FUNCS([floor modf pow sqrt])

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 test/Makefile
                 doc/Makefile])
AC_OUTPUT
