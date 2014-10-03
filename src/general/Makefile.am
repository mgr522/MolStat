AUTOMAKE_OPTIONS = subdir-objects

noinst_LIBRARIES = libmolstat_general.a

SUBDIRS = random_distributions \
	fitter_tools \
	histogram_tools \
	simulator_tools

libmolstat_general_a_SOURCES = string_tools.h \
	string_tools.cc \
	random_distributions/rng.h \
	random_distributions/rng.cc \
	random_distributions/constant.h \
	random_distributions/constant.cc \
	random_distributions/uniform.h \
	random_distributions/uniform.cc \
	random_distributions/normal.h \
	random_distributions/normal.cc \
	random_distributions/lognormal.h \
	random_distributions/lognormal.cc \
	random_distributions/gamma.h \
	random_distributions/gamma.cc \
	fitter_tools/fit_model_interface.h \
	fitter_tools/fit_model_interface.cc \
	histogram_tools/bin_style.h \
	histogram_tools/bin_style.cc \
	histogram_tools/bin_linear.h \
	histogram_tools/bin_linear.cc \
	histogram_tools/bin_log.h \
	histogram_tools/bin_log.cc \
	histogram_tools/histogram_interface.h \
	histogram_tools/histogram1d.h \
	histogram_tools/histogram1d.cc \
	histogram_tools/histogram2d.h \
	histogram_tools/histogram2d.cc \
	simulator_tools/simulate_model_interface.h
