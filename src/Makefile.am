SUBDIRS = general \
	electron_transport

bin_PROGRAMS = molstat-simulator molstat-fitter

molstat_simulator_SOURCES = main-simulator.cc

molstat_simulator_LDADD = general/libmolstat_general.a \
	electron_transport/simulator_models/libtransport_simulate.a

molstat_fitter_SOURCES = main-fitter.cc

molstat_fitter_LDADD = general/libmolstat_general.a \
	electron_transport/fitter_models/libtransport_fit.a
