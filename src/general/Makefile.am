AUTOMAKE_OPTIONS = subdir-objects

AM_CPPFLAGS = -I$(top_srcdir)/src

noinst_LIBRARIES = \
	libmolstat_general.a

SUBDIRS = random_distributions \
	fitter_tools \
	histogram_tools \
	simulator_tools \
	tests

libmolstat_general_a_SOURCES = \
	string_tools.h \
	string_tools.cc \
	histogram_tools/counterindex.h \
	histogram_tools/counterindex.cc \
	histogram_tools/bin_style.h \
	histogram_tools/bin_style.cc \
	histogram_tools/bin_linear.h \
	histogram_tools/bin_linear.cc \
	histogram_tools/bin_log.h \
	histogram_tools/bin_log.cc \
	histogram_tools/histogram.h \
	histogram_tools/histogram.cc

if BUILD_SIMULATOR
noinst_LIBRARIES += libmolstat_simulator.a

libmolstat_simulator_a_SOURCES = \
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
	simulator_tools/simulator_exceptions.h \
	simulator_tools/simulator.h \
	simulator_tools/simulator.cc \
	simulator_tools/simulate_model.h \
	simulator_tools/observable.h \
	simulator_tools/simulate_model.cc \
	simulator_tools/composite_simulate_model.cc \
	simulator_tools/simulate_model_factory.cc \
	simulator_tools/identity_tools.h \
	simulator_tools/identity_tools.cc
endif

if BUILD_FITTER
noinst_LIBRARIES += libmolstat_fitter.a

libmolstat_fitter_a_SOURCES = \
	fitter_tools/fit_model_interface.h \
	fitter_tools/fit_model_interface.cc

# fitter requires GSL
libmolstat_fitter_a_CPPFLAGS = $(GSL_INCLUDE) $(AM_CPPFLAGS)
endif
