SUBDIRS = simulator_models \
	fitter_models \
	tests

bin_PROGRAMS = fitter

fitter_SOURCES = main-fitter.cc

fitter_LDADD = ../general/libmolstat_general.a \
	fitter_models/libfitter_models.a
