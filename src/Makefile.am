SUBDIRS = general \
	electron_transport \
	tests

bin_PROGRAMS = molstat-simulator molstat-fitter

molstat_simulator_SOURCES = main-simulator.h \
	main-simulator.cc \
	main-simulator-inputparse.cc

molstat_simulator_LDADD = \
	electron_transport/simulator_models/libtransport_simulate.a \
	general/libmolstat_general.a 

molstat_fitter_SOURCES = main-fitter.cc

molstat_fitter_LDADD = \
	electron_transport/fitter_models/libtransport_fit.a \
	general/libmolstat_general.a 
