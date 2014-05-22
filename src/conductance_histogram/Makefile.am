SUBDIRS = simulator_models \
	fitter_models \
	tests

bin_PROGRAMS = simulator fitter

AM_CPPFLAGS = -I../general

simulator_SOURCES = main-simulator.cc

simulator_LDADD = ../general/libaux.a \
	../general/histogram_tools/libhistogram_tools.a \
	../general/random_distributions/libaux_random_distributions.a \
	simulator_models/libsimulator_models.a

fitter_SOURCES = main-fitter.cc

fitter_LDADD = ../general/libaux.a \
	../general/histogram_tools/libhistogram_tools.a \
	fitter_models/libfitter_models.a
