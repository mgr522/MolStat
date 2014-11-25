AC_DEFUN([ACX_WITH_CVODE],
[
  acx_with_cvode=no
  AC_ARG_WITH([cvode],
    [AS_HELP_STRING([--with-cvode@<:@=Install DIR@:>@], [Build with CVODE.])],
    [
      acx_with_cvode=yes
      CVODE_INCLUDE="-I$withval/include"
      CVODE_LDFLAGS="-L$withval/lib"
      CVODE_LIBS="-lsundials_cvode -lsundials_nvecserial -lm"
    ]
  )
  
  if test "$acx_with_cvode" != no; then
    CPPFLAGS="$CVODE_INCLUDE $CPPFLAGS"
    LDFLAGS="$CVODE_LDFLAGS $LDFLAGS"
    LIBS="$CVODE_LIBS $LIBS"

    # Check for the pressence of the necessary CVODE headers.
    AC_CHECK_HEADER([cvode/cvode.h], [],
      [AC_MSG_ERROR([Unable to find the cvode/cvode.h header file.])])
    AC_CHECK_HEADER([nvector/nvector_serial.h], [],
      [AC_MSG_ERROR([Unable to find the nvector/nvector_serial.h header file.])])
    AC_CHECK_HEADER([cvode/cvode_dense.h], [],
      [AC_MSG_ERROR([Unable to find the cvode/cvode_dense.h header file.])])
    AC_CHECK_HEADER([sundials/sundials_dense.h], [],
      [AC_MSG_ERROR([Unable to find the sundials/sundials_dense.h header file.])])
    AC_CHECK_HEADER([sundials/sundials_types.h], [],
      [AC_MSG_ERROR([Unable to find the sundials/sundials_types.h header file.])])

    # Check for the CVODE libraries
    AC_CHECK_LIB([m], [cos])
    AC_CHECK_LIB([sundials_cvode], [CVodeCreate])
    AC_CHECK_LIB([sundials_nvecserial], [main])
  else
    AC_MSG_ERROR([Unable to find CVODE.])
  fi
])
