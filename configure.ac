# boilerplate initialization
AC_INIT([MolStat], [1.2], [mgreuter@u.northwestern.edu], , [https://bitbucket.org/mgreuter/molstat])
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
# fitter
AC_ARG_ENABLE([transport-fitter],
	[AS_HELP_STRING([--disable-transport-fitter],
		[disable the transport module for the fitter @<:@default: no@:>@])],
	[with_transport_fit=${enableval}], [with_transport_fit=yes])
# disable the transport module if we're not building the fitter
if test x$build_fitter = xno; then
	with_transport_fit=no
fi

# simulator
AC_ARG_ENABLE([transport-simulator],
	[AS_HELP_STRING([--disable-transport-simulator],
		[disable the transport module for the simulator @<:@default: no@:>@])],
	[with_transport_sim=${enableval}], [with_transport_sim=yes])
# disable the transport module if we're not building the simulator
if test x$build_simulator = xno; then
	with_transport_sim=no
fi


# by default, include the electrochemistry modules
# fitter
AC_ARG_ENABLE([echem-fitter],
	[AS_HELP_STRING([--disable-echem-fitter],
		[disable the electrochemistry module for the fitter @<:@default: no@:>@]
		)],
	[with_echem_fit=${enableval}], [with_echem_fit=yes])
# disable the electrochemistry module if we're not building the fitter
if test x$build_fitter = xno; then
	with_echem_fit=no
fi

# simulator
AC_ARG_ENABLE([echem-simulator],
	[AS_HELP_STRING([--disable-echem-simulator],
		[disable the electrochemistry module for the simulator @<:@default: no@:>@])],
	[with_echem_sim=${enableval}], [with_echem_sim=yes])
# disable the electrochemistry module if we're not building the simulator
if test x$build_simulator = xno; then
	with_echem_sim=no
fi



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


# parts of the electrochemistry simulator module require CVODE
# include them if CVODE is specified, but skip them otherwise
if test x$with_echem_sim = xyes; then
	need_cvode=maybe
fi

# electrochemistry fitter module requires GSL
if test x$with_echem_fit = xyes; then
	need_gsl=yes
fi



# check for the necessary packages

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
AM_CONDITIONAL([HAVE_GSL],
	[test x$need_gsl != no && test x$acx_with_gsl = xyes])
if test x$need_gsl != no && test x$acx_with_gsl = xyes; then
	AC_DEFINE([HAVE_GSL], [1],
		[We can build and link with GSL.])
else
	AC_DEFINE([HAVE_GSL], [0],
		[We can build and link with GSL.])
fi

# look for CVODE, if necessary
if test x$need_cvode != xno; then
	ACX_WITH_CVODE
	if test x$acx_with_cvode = xno; then
		if test x$need_cvode = xyes; then
			AC_MSG_ERROR([Unable to find CVODE.])
		fi
		if test x$need_cvode = xmaybe; then
			AC_MSG_WARN([Unable to find CVODE. Requested parts of MolStat that do not require CVODE will still be built.])
		fi
	fi
fi
AM_CONDITIONAL([HAVE_CVODE],
	[test x$need_cvode != no && test x$acx_with_cvode = xyes])
if test x$need_cvode != no && test x$acx_with_cvode = xyes; then
	AC_DEFINE([HAVE_CVODE], [1],
		[We can build and link with CVODE.])
else
	AC_DEFINE([HAVE_CVODE], [0],
		[We can build and link with CVODE.])
fi


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



# electrochemistry flags
AM_CONDITIONAL([ECHEM_FITTER], [test x$with_echem_fit = xyes])
AM_CONDITIONAL([ECHEM_SIMULATOR], [test x$with_echem_sim = xyes])

if test x$with_echem_fit = xyes; then
	AC_DEFINE([BUILD_ECHEM_FITTER], [1],
		[Include the electrochemistry fitter module.])
else
	AC_DEFINE([BUILD_ECHEM_FITTER], [0],
		[Include the electrochemistry fitter module.])
fi

if test x$with_echem_sim = xyes; then
	AC_DEFINE([BUILD_ECHEM_SIMULATOR], [1],
		[Include the electrochemistry simulator module.])
else
	AC_DEFINE([BUILD_ECHEM_SIMULATOR], [0],
		[Include the electrochemistry simulator module.])
fi


# configure the documentation
# by default, do NOT build the documentation
AC_ARG_ENABLE([documentation],
	[AS_HELP_STRING([--enable-documentation],
		[enable building the documentation (for developers only) @<:@default: no@:>@])],
	[build_documentation=${enableval}], [build_documentation=no])

# check for configuration files, doxygen, and pdflatex, if necessary
if test x$build_documentation = xyes; then
	AC_CHECK_FILE([doc/Makefile.am], ,
		[AC_MSG_WARN([Unable to build documentation: Cannot find doc/Makefile.am])
		build_documentation=no])

	AC_CHECK_FILE([doc/userman.doxygen.in], ,
		[AC_MSG_WARN([Unable to build documentation: Cannot find doc/userman.doxygen.in])
		build_documentation=no])

	AC_CHECK_FILE([doc/fullref.doxygen.in], ,
		[AC_MSG_WARN([Unable to build documentation: Cannot find doc/fullref.doxygen.in])
		build_documentation=no])

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
			src/echem/Makefile
				src/echem/simulator_models/Makefile
				src/echem/tests/Makefile
			src/electron_transport/Makefile
				src/electron_transport/fitter_models/Makefile
				src/electron_transport/simulator_models/Makefile
				src/electron_transport/tests/Makefile
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
