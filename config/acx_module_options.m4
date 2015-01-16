# Given a module name, checks for the disable/enable options and sets the
# appropriate flags.
#
# By default, both the simulator and fitter modules are included.
#
# Parameter $1: The name of the module (in English).
# Parameter $2: The name of the module (directory).
#
# No return value; sets the with_$2_fit and with_$2_sim variables telling
# us whether to build and link the specified module.
#
AC_DEFUN([ACX_MODULE_OPTIONS],
[
	# check for the fitter option
	AC_ARG_ENABLE([$2-fitter],
		[AS_HELP_STRING([--disable-$2-fitter],
			[disable the $1 module for the fitter @<:@default: no@:>@])],
		[with_$2_fit=${enableval}], [with_$2_fit=yes])
	# disable the module if we're not building the fitter
	if test x$build_fitter = xno; then
		with_$2_fit=no
	fi

	# check for the simulator option
	AC_ARG_ENABLE([$2-simulator],
		[AS_HELP_STRING([--disable-$2-simulator],
			[disable the $1 module for the simulator @<:@default: no@:>@])],
		[with_$2_sim=${enableval}], [with_$2_sim=yes])
	# disable the module if we're not building the simulator
	if test x$build_simulator = xno; then
		with_$2_sim=no
	fi
])
