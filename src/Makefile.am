SUBDIRS = general \
	electron_transport \
	tests

bin_PROGRAMS =

if BUILD_SIMULATOR
bin_PROGRAMS += molstat-simulator

molstat_simulator_SOURCES = main-simulator.h \
	main-simulator.cc \
	main-simulator-inputparse.cc

molstat_simulator_LDADD =

if TRANSPORT_SIMULATOR
molstat_simulator_LDADD += \
	electron_transport/simulator_models/libtransport_simulate.a
endif

molstat_simulator_LDADD += \
	general/libmolstat_simulator.a \
	general/libmolstat_general.a \
	$(AM_LDADD)
endif

if BUILD_FITTER
bin_PROGRAMS += molstat-fitter

molstat_fitter_SOURCES = main-fitter.cc

# fitter requires GSL
molstat_fitter_CPPFLAGS = $(GSL_INCLUDE) $(AM_CPPFLAGS)

molstat_fitter_LDADD =

if TRANSPORT_FITTER
molstat_fitter_LDADD += \
	electron_transport/fitter_models/libtransport_fit.a
endif

molstat_fitter_LDADD += \
	general/libmolstat_fitter.a \
	general/libmolstat_general.a \
	$(GSL_LDFLAGS) $(AM_LDADD) $(GSL_LIBS) $(AM_LIBS)
endif
