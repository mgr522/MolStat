# boilerplate initialization
AC_INIT([MolStat], [1.3.1], [matthew.reuter@stonybrook.edu], , [https://bitbucket.org/mgreuter/molstat])
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

# required for building static libraries
AC_PROG_RANLIB

# python is needed for some "make check" scripts
AX_PYTHON
AM_CONDITIONAL([HAVE_PYTHON], [test "$PYTHON" != ":"])



# check the options.
# by default, build the fitter
AC_ARG_ENABLE([fitter],
	[AS_HELP_STRING([--disable-fitter],
		[disable compilation of the fitter @<:@default: no@:>@])],
	[build_fitter=${enableval}], [build_fitter=yes])

# by default, build the simulator
AC_ARG_ENABLE([simulator],
	[AS_HELP_STRING([--disable-simulator],
		[disable compilation of the simulator @<:@default: no@:>@])],
	[build_simulator=${enableval}], [build_simulator=yes])

# transport module
ACX_MODULE_OPTIONS([transport], [transport])



# figure out what software packages we require
need_gsl=no
need_cvode=no

# the simulator itself requires nothing

# the fitter itself requires GSL
if test x$build_fitter = xyes; then
	need_gsl=yes
fi


# transport simulator module requires nothing


# transport fitter module requires GSL
if test x$with_transport_fit = xyes; then
	need_gsl=yes
fi



# check for the necessary packages

# look for GSL
ACX_WITH_GSL
ACX_SET_PACKAGE([gsl], [GSL])

# look for CVODE, if necessary
ACX_WITH_CVODE
ACX_SET_PACKAGE([cvode], [CVODE])


# set flags for all of the options
# the AM_CONDITIONAL flags are for the makefiles, the AC_DEFINE flags get put
# in config.h for the code to use
AM_CONDITIONAL([BUILD_FITTER], [test x$build_fitter = xyes])
AM_CONDITIONAL([BUILD_SIMULATOR], [test x$build_simulator = xyes])


# transport flags
AM_CONDITIONAL([TRANSPORT_FITTER], [test x$with_transport_fit = xyes])
AM_CONDITIONAL([TRANSPORT_SIMULATOR], [test x$with_transport_sim = xyes])

if test x$with_transport_fit = xyes; then
	AC_DEFINE([BUILD_TRANSPORT_FITTER], [1],
		[Include the transport fitter module.])
else
	AC_DEFINE([BUILD_TRANSPORT_FITTER], [0],
		[Include the transport fitter module.])
fi

if test x$with_transport_sim = xyes; then
	AC_DEFINE([BUILD_TRANSPORT_SIMULATOR], [1],
		[Include the transport simulator module.])
else
	AC_DEFINE([BUILD_TRANSPORT_SIMULATOR], [0],
		[Include the transport simulator module.])
fi



# configure the documentation
# by default, do NOT build the documentation
AC_ARG_ENABLE([documentation],
	[AS_HELP_STRING([--enable-documentation],
		[enable building the documentation (for developers only) @<:@default: no@:>@])],
	[build_documentation=${enableval}], [build_documentation=no])

# check for configuration files, doxygen, and pdflatex, if necessary
if test x$build_documentation = xyes; then
	AC_CHECK_PROG([have_doxygen], [doxygen], [yes], [no])
	if test x$have_doxygen = xyes; then
	   AC_PATH_PROG(DOXYGEN, doxygen)
	else
		AC_MSG_WARN([Unable to build documentation: Cannot find doxygen.])
	   build_documentation=no
	fi

	AC_CHECK_PROG([have_pdflatex], [pdflatex], [yes], [no])
	if test x$have_pdflatex = xyes; then
   	AC_PATH_PROG(PDFLATEX, pdflatex)
   else
   	AC_MSG_WARN([Unable to build documentation: Cannot find pdflatex.])
   	build_documentation=no
   fi
fi
AM_CONDITIONAL([BUILD_DOC], [test x$build_documentation = xyes])

# output for the config header file
AC_CONFIG_HEADERS([src/config.h])

# files to be generated
AC_CONFIG_FILES([
	README
	Makefile
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
					src/electron_transport/tests/fit-asymmetric-resonant.py
					src/electron_transport/tests/fit-options.py
					src/electron_transport/tests/fit-symmetric-nonresonant.py
					src/electron_transport/tests/fit-symmetric-nonresonant-binlog.py
					src/electron_transport/tests/fit-symmetric-resonant.py
			src/tests/Makefile
])

if test x$build_documentation = xyes; then
	AC_CONFIG_FILES([
		doc/Makefile
			doc/fullref.doxygen
			doc/userman.doxygen
	])
fi

# finish up
AC_OUTPUT
