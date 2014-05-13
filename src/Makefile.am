SUBDIRS = aux_random_distributions \
	aux_simulator \
	aux_binner \
	aux_fitter \
	tests

bin_PROGRAMS = simulator fitter

simulator_SOURCES = main-simulator.cc \
	string_tools.h \
	string_tools.cc

simulator_LDADD = aux_simulator/libaux_simulator.a \
	aux_random_distributions/libaux_random_distributions.a \
	aux_binner/libaux_binner.a

fitter_SOURCES = main-fitter.cc

fitter_LDADD = aux_binner/libaux_binner.a \
	aux_fitter/libaux_fitter.a
