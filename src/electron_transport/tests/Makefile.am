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

TESTS += \
	fit-AsymResonant \
	fit-SymNonresonant \
	fit-SymResonant

check_PROGRAMS += \
	fit-AsymResonant \
	fit-SymNonresonant \
	fit-SymResonant

fit_AsymResonant_SOURCES = fit-AsymResonant.cc
fit_AsymResonant_CPPFLAGS = $(GSL_INCLUDE) $(AM_CPPFLAGS)
fit_AsymResonant_LDADD = ../fitter_models/libtransport_fit.a \
	../../general/libmolstat_fitter.a \
	../../general/libmolstat_general.a \
	$(GSL_LDFLAGS) $(AM_LDADD) $(GSL_LIBS) $(AM_LIBS)

fit_SymNonresonant_SOURCES = fit-SymNonresonant.cc
fit_SymNonresonant_CPPFLAGS = $(GSL_INCLUDE) $(AM_CPPFLAGS)
fit_SymNonresonant_LDADD = ../fitter_models/libtransport_fit.a \
	../../general/libmolstat_fitter.a \
	../../general/libmolstat_general.a \
	$(GSL_LDFLAGS) $(AM_LDADD) $(GSL_LIBS) $(AM_LIBS)

fit_SymResonant_SOURCES = fit-SymResonant.cc
fit_SymResonant_CPPFLAGS = $(GSL_INCLUDE) $(AM_CPPFLAGS)
fit_SymResonant_LDADD = ../fitter_models/libtransport_fit.a \
	../../general/libmolstat_fitter.a \
	../../general/libmolstat_general.a \
	$(GSL_LDFLAGS) $(AM_LDADD) $(GSL_LIBS) $(AM_LIBS)
endif

# make sure automake includes the scripts in a distribution
dist_check_SCRIPTS = \
	fit-symmetric-resonant.py \
	fit-symmetric-nonresonant.py \
	fit-symmetric-nonresonant-binlog.py \
	fit-options.py \
	fit-asymmetric-resonant.py

dist_check_DATA = \
	asymmetric-resonant.dat	\
	symmetric-nonresonant.dat \
	symmetric-nonresonant-binlog.dat \
	symmetric-resonant.dat
