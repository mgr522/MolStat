TESTS =
check_PROGRAMS =

TEST_EXTENSIONS = .py
PY_LOG_COMPILER = $(PYTHON_BIN)

if TRANSPORT_SIMULATOR
TESTS += \
	simulate-SymOneSite \
	simulate-AsymOneSite \
	simulate-SymTwoSite \
	simulate-AsymTwoSite \
	simulate-CompositeJunction

check_PROGRAMS += \
	simulate-SymOneSite \
	simulate-AsymOneSite \
	simulate-SymTwoSite \
	simulate-AsymTwoSite \
	simulate-CompositeJunction

simulate_SymOneSite_SOURCES = simulate-SymOneSite.cc
simulate_SymOneSite_LDADD = ../simulator_models/libtransport_simulate.a \
	../../general/libmolstat_simulator.a \
	../../general/libmolstat_general.a 

simulate_AsymOneSite_SOURCES = simulate-AsymOneSite.cc
simulate_AsymOneSite_LDADD = ../simulator_models/libtransport_simulate.a \
	../../general/libmolstat_simulator.a \
 	../../general/libmolstat_general.a 

simulate_SymTwoSite_SOURCES = simulate-SymTwoSite.cc
simulate_SymTwoSite_LDADD = ../simulator_models/libtransport_simulate.a \
	../../general/libmolstat_simulator.a \
	../../general/libmolstat_general.a 

simulate_AsymTwoSite_SOURCES = simulate-AsymTwoSite.cc
simulate_AsymTwoSite_LDADD = ../simulator_models/libtransport_simulate.a \
	../../general/libmolstat_simulator.a \
	../../general/libmolstat_general.a 

simulate_CompositeJunction_SOURCES = simulate-CompositeJunction.cc
simulate_CompositeJunction_LDADD = ../simulator_models/libtransport_simulate.a \
	../../general/libmolstat_simulator.a \
	../../general/libmolstat_general.a 
endif

if TRANSPORT_FITTER
if HAVE_PYTHON
TESTS += \
	fit-symmetric-resonant.py \
	fit-symmetric-nonresonant.py \
   fit-symmetric-nonresonant-binlog.py \
	fit-options.py \
	fit-asymmetric-resonant.py
endif
endif

# make sure automake includes the script in a distribution
dist_check_SCRIPTS = \
	fit-symmetric-resonant.py \
	fit-symmetric-nonresonant.py \
	fit-symmetric-nonresonant-binlog.py \
	fit-options.py \
	fit-asymmetric-resonant.py
