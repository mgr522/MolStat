if ECHEM_SIMULATOR
noinst_LIBRARIES = libechem_simulate.a

libechem_simulate_a_SOURCES = \
	echem_simulate_module.h \
	echem_simulate_module.cc \
	observables.h \
	non-nernstian.h \
	non-nernstian.cc

# parts of the electrochemistry simulator module require CVODE
libechem_simulate_a_CPPFLAGS = $(CVODE_INCLUDE) $(AM_CPPFLAGS)
endif
