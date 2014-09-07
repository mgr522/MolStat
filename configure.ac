# boilerplate initialization
AC_INIT([MolStat], [1.0], [mgreuter@u.northwestern.edu])
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

# python is needed for some "make check" scripts
AX_PYTHON
AM_CONDITIONAL([HAVE_PYTHON], [test "$PYTHON" != ":"])

# look for GSL
ACX_WITH_GSL

# required options for building static libraries
AC_PROG_RANLIB

# output for the config header file
AC_CONFIG_HEADERS([src/config.h])

# We need these to build the documentation
AC_CHECK_PROG([HAVE_DOXYGEN],[doxygen],[YES],[NO])
if test $HAVE_DOXYGEN = YES; then
   AC_PATH_PROG(DOXYGEN,doxygen)
fi
AC_CHECK_PROG([HAVE_PDFLATEX],[pdflatex],[YES],[NO])
if test $HAVE_PDFLATEX = YES; then
   AC_PATH_PROG(PDFLATEX,pdflatex)
fi
AM_CONDITIONAL([CAN_BUILD_DOC], [test $HAVE_DOXYGEN = YES && test $HAVE_PDFLATEX = YES])

# files to be generated
AC_CONFIG_FILES([
	Makefile
		doc/Makefile
		doc/fullref.doxygen
		doc/userman.doxygen
		src/Makefile
			src/general/Makefile
				src/general/fitter_tools/Makefile
				src/general/histogram_tools/Makefile
				src/general/random_distributions/Makefile
				src/general/simulator_tools/Makefile
			src/conductance_histogram/Makefile
				src/conductance_histogram/fitter_models/Makefile
				src/conductance_histogram/simulator_models/Makefile
				src/conductance_histogram/tests/Makefile
])

# finish up
AC_OUTPUT
