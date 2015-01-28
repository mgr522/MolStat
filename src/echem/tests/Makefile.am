TESTS =
check_PROGRAMS =


if ECHEM_SIMULATOR
TESTS += \
	simulate-NonNernstian 

check_PROGRAMS += \
	simulate-NonNernstian 


simulate_NonNernstian_SOURCES = simulate-NonNernstian.cc 

simulate_NonNernstian_LDADD = ../simulator_models/libechem_simulate.a \
	../../general/libmolstat_simulator.a \
	../../general/libmolstat_general.a \
	$(CVODE_LIBS) 
endif
