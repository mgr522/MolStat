SUBDIRS = general \
	electron_transport

bin_PROGRAMS = molstat-simulator

molstat_simulator_SOURCES = main-simulator.cc

molstat_simulator_LDADD = general/libmolstat_general.a \
	electron_transport/simulator_models/libtransport_simulator_models.a
