# boilerplate initialization
AC_INIT([cond-hist], [1.0], [mgreuter@u.northwestern.edu])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])

# local macros
AC_CONFIG_MACRO_DIR([config])

m4_pattern_allow([AM_PROG_AR])
AM_PROG_AR

# specify C++11
AC_LANG([C++])
AX_CXX_COMPILE_STDCXX_11([noext], [mandatory])

# look for GSL
AC_CHECK_LIB([m], [cos])
AC_CHECK_LIB([gslcblas], [cblas_dgemm])
AC_CHECK_LIB([gsl], [gsl_blas_dgemm])

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
