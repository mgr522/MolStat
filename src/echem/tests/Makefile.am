TESTS =
check_PROGRAMS =

if ECHEM_SIMULATOR
TESTS += \
	simulate-Nernstian 

check_PROGRAMS += \
	simulate-Nernstian

simulate_Nernstian_SOURCES = simulate-Nernstian.cc
simulate_Nernstian_LDADD = ../simulator_models/libechem_simulate.a \
	../../general/libmolstat_simulator.a \
	../../general/libmolstat_general.a

# test that requires CVODE
if HAVE_CVODE

TESTS += \
	simulate-NonNernstian 

check_PROGRAMS += \
	simulate-NonNernstian

simulate_NonNernstian_SOURCES = simulate-NonNernstian.cc 
simulate_NonNernstian_CPPFLAGS = $(CVODE_INCLUDE)
simulate_NonNernstian_LDADD = ../simulator_models/libechem_simulate.a \
	../../general/libmolstat_simulator.a \
	../../general/libmolstat_general.a \
	$(CVODE_LDFLAGS) $(AM_LDADD) $(CVODE_LIBS) $(AM_LIBS)
endif # HAVE_CVODE

endif # ECHEM_SIMULATOR
