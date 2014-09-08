SUBDIRS = simulator_models \
	fitter_models \
	tests

bin_PROGRAMS = fitter

#simulator_SOURCES = main-simulator.cc

#simulator_LDADD = ../general/libmolstat_general.a \
#	simulator_models/libsimulator_models.a

fitter_SOURCES = main-fitter.cc

fitter_LDADD = ../general/libmolstat_general.a \
	fitter_models/libfitter_models.a
