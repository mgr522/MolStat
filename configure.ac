# boilerplate initialization
AC_INIT([MolStat], [1.2], [mgreuter@u.northwestern.edu])
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



# check the options.
# by default, build both the simulator and the fitter
AC_ARG_ENABLE([fitter],
	[AS_HELP_STRING([--disable-fitter],
		[disable compilation of the fitter @<:@default: no@:>@])],
	[build_fitter=${enableval}], [build_fitter=yes])
AC_ARG_ENABLE([simulator],
	[AS_HELP_STRING([--disable-simulator],
		[disable compilation of the simulator @<:@default: no@:>@])],
	[build_simulator=${enableval}], [build_simulator=yes])

# by default, include the transport modules
AC_ARG_ENABLE([transport-fitter],
	[AS_HELP_STRING([--disable-transport-fitter],
		[disable the transport module for the fitter @<:@default: no@:>@])],
	[with_transport_fit=${enableval}], [with_transport_fit=yes])
# disable the transport module if we're not building the fitter
if test x$build_fitter = xno; then
	with_transport_fit=no
fi

AC_ARG_ENABLE([transport-simulator],
	[AS_HELP_STRING([--disable-transport-simulator],
		[disable the transport module for the simulator @<:@default: no@:>@])],
	[with_transport_sim=${enableval}], [with_transport_sim=yes])
# disable the transport module if we're not building the simulator
if test x$build_simulator = xno; then
	with_transport_sim=no
fi



# figure out what software packages we require
need_gsl=no

# the simulator requires nothing

# the fitter requires GSL
if test x$build_fitter = xyes; then
	need_gsl=yes
fi

# simulator transport module requires nothing

# fitter transport module requires GSL
if test x$with_transport_fit = xyes; then
	need_gsl=yes
fi



# check for the necessary packages

# required for building static libraries
AC_PROG_RANLIB

# python is needed for some "make check" scripts
AX_PYTHON
AM_CONDITIONAL([HAVE_PYTHON], [test "$PYTHON" != ":"])

# look for GSL, if necessary
if test x$need_gsl != xno; then
	ACX_WITH_GSL
	if test x$acx_with_gsl = xno; then
		if test x$need_gsl = xyes; then
			AC_MSG_ERROR([Unable to find GSL.])
		fi
		if test x$need_gsl = xmaybe; then
			AC_MSG_WARN([Unable to find GSL. Requested parts of MolStat that do not require GSL will still be built.])
		fi
	fi
fi



# set flags for all of the options
AM_CONDITIONAL([BUILD_FITTER], [test x$build_fitter = xyes])
AM_CONDITIONAL([TRANSPORT_FITTER], [test x$with_transport_fit = xyes])

AM_CONDITIONAL([BUILD_SIMULATOR], [test x$build_simulator = xyes])
AM_CONDITIONAL([TRANSPORT_SIMULATOR], [test x$with_transport_sim = xyes])



# output for the config header file
AC_CONFIG_HEADERS([src/config.h])

# We need these to build the documentation
AC_CHECK_PROG([HAVE_DOXYGEN], [doxygen], [YES], [NO])
if test $HAVE_DOXYGEN = YES; then
   AC_PATH_PROG(DOXYGEN, doxygen)
fi
AC_CHECK_PROG([HAVE_PDFLATEX], [pdflatex], [YES], [NO])
if test $HAVE_PDFLATEX = YES; then
   AC_PATH_PROG(PDFLATEX, pdflatex)
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
				src/general/tests/Makefile
			src/electron_transport/Makefile
				src/electron_transport/fitter_models/Makefile
				src/electron_transport/simulator_models/Makefile
				src/electron_transport/tests/Makefile
			src/tests/Makefile
])

# finish up
AC_OUTPUT
