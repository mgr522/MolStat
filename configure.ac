# boilerplate initialization
AC_INIT([cond-hist], [1.0], [mgreuter@u.northwestern.edu])
AC_CONFIG_SRCDIR([configure.ac])

# this is where autoconf-generated files go
AC_CONFIG_AUX_DIR([config])

# local macros
AC_CONFIG_MACRO_DIR([config])

# initialize automake
AM_INIT_AUTOMAKE([foreign -Wall -Werror])
m4_pattern_allow([AM_PROG_AR])
AM_PROG_AR

# specify C++11
AC_LANG([C++])
AX_CXX_COMPILE_STDCXX_11([noext], [mandatory])

# look for GSL
ACX_WITH_GSL

# required options for building static libraries
AC_PROG_RANLIB

# output for the config header file
AC_CONFIG_HEADERS([src/config.h])

# files to be generated
AC_CONFIG_FILES([
	Makefile
		src/Makefile
			src/aux_simulator/Makefile
])

# finish up
AC_OUTPUT
