# Set the flags for an external package.
#
# Parameter $1: Name of the package (in English, probably lowercase).
# Parameter $2: Name of the package (in uppercase, for use in macro names).
#
# Assumes that the corresponding ACX_WITH_$1 macro has been run such that
# $acx_with_$1 is set.
# 
AC_DEFUN([ACX_SET_PACKAGE],
[
  have_$1=no
  
  if test x$need_$1 != xno; then
    if test x$acx_with_$1 = xno; then
      if test x$need_$1 = xyes; then
        AC_MSG_ERROR([Unable to find $1.])
      fi
      if test x$need_$1 = xmaybe; then
        AC_MSG_WARN([Unable to find $1. Requested parts of MolStat that do not require $1 will still be built.])
      fi
    else
      have_$1=yes
    fi
  fi

  # automake macro
  AM_CONDITIONAL([HAVE_$2], [test x$have_$1 = xyes])

  # C++ define macro
  if test x$have_$1 = xyes; then
    AC_DEFINE([HAVE_$2], [1], [We can build and link with $1.])
  else
    AC_DEFINE([HAVE_$2], [0], [We can build and link with $1.])
  fi
])
