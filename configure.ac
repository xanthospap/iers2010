#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([libiers10++], [1.0], [xanthos@mail.ntua.gr])
# AM_INIT_AUTOMAKE

AC_LANG(C++)

AC_CONFIG_SRCDIR([src/gpt.cpp])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX_11([noext], [mandatory])
AC_PROG_INSTALL
AC_PROG_RANLIB

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
                 test/Makefile])
AC_OUTPUT

echo "
  ($PACKAGE_NAME) version $PACKAGE_VERSION
  Prefix.........: $prefix
  Debug Build....: $debug
  C++ Compiler...: $CXX $CXXFLAGS $CPPFLAGS
  Linker.........: $LD $LDFLAGS $LIBS
"
