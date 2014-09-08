SUBDIRS = general \
	electron_transport

bin_PROGRAMS = simulator

simulator_SOURCES = main-simulator.cc

simulator_LDADD = general/libmolstat_general.a \
	electron_transport/simulator_models/libtransport_simulator_models.a
