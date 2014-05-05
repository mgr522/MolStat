SUBDIRS = aux_random_distributions aux_simulator aux_binner

bin_PROGRAMS = simulator

simulator_SOURCES = main-simulator.cc \
	string_tools.h \
	string_tools.cc

simulator_LDADD = aux_simulator/libaux_simulator.a \
	aux_random_distributions/libaux_random_distributions.a \
	aux_binner/libaux_binner.a
