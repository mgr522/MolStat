# Check for the GSL.
#
# Sets $acx_with_gsl to yes if GSL is found and in working order; no otherwise
#
# If the GSL is found, sets GSL_INCLUDE, GSL_LDFLAGS, and GSL_LIBS
AC_DEFUN([ACX_WITH_GSL],
[
  acx_with_gsl=no
  AC_ARG_WITH([gsl],
    [AS_HELP_STRING([--with-gsl@<:@=Install DIR@:>@], [Build with GSL.])],
    [
      acx_with_gsl=yes
      GSL_INCLUDE="-I$withval/include"
      GSL_LDFLAGS="-L$withval/lib"
      GSL_LIBS="-lgsl -lgslcblas -lm"
    ]
  )
  
  if test x"$acx_with_gsl" != xno; then
    # store the existing CPPFLAGS and LDFLAGS
    my_CPPFLAGS=$CPPFLAGS
    my_LDFLAGS=$LDFLAGS
    my_LIBS=$LIBS
    
    # append the GSL information and perform the tests
    CPPFLAGS="$GSL_INCLUDE $CPPFLAGS"
    LDFLAGS="$GSL_LDFLAGS $LDFLAGS"

    # Check for the pressence of the necessary GSL headers
    AC_CHECK_HEADER([gsl/gsl_blas.h],
      [],
      [AC_MSG_WARN([Unable to find the gsl/gsl_blas.h header file.])
       acx_with_gsl=no])
    AC_CHECK_HEADER([gsl/gsl_vector.h], [],
      [AC_MSG_WARN([Unable to find the gsl/gsl_vector.h header file.])
       acx_with_gsl=no])
    AC_CHECK_HEADER([gsl/gsl_matrix.h], [],
      [AC_MSG_WARN([Unable to find the gsl/gsl_matrix.h header file.])
       acx_with_gsl=no])
    AC_CHECK_HEADER([gsl/gsl_multifit_nlin.h], [],
      [AC_MSG_WARN([Unable to find the gsl/gsl_multifit_nlin.h header file.])
       acx_with_gsl=no])
    AC_CHECK_HEADER([gsl/gsl_integration.h], [],
      [AC_MSG_WARN([Unable to find the gsl/gsl_integration.h header file.])
       acx_with_gsl=no])

    # Check for the GSL libraries
    AC_CHECK_LIB([m], [cos])
    AC_CHECK_LIB([gslcblas], [cblas_dgemm])
    AC_CHECK_LIB([gsl], [gsl_blas_dgemm])

    # revert the CPPFLAGS and LDFLAGS
    CPPFLAGS=$my_CPPFLAGS
    LDFLAGS=$my_LDFLAGS
    LIBS=$my_LIBS
  fi

  if test x"$acx_with_gsl" != xno; then
    AC_SUBST([GSL_INCLUDE], [$GSL_INCLUDE])
    AC_SUBST([GSL_LDFLAGS], [$GSL_LDFLAGS])
    AC_SUBST([GSL_LIBS], [$GSL_LIBS])
  fi
])
