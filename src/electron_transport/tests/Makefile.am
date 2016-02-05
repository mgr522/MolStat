AM_CPPFLAGS = -I$(top_srcdir)/src

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
	simulate-SymInterference \
	simulate-RectBarrier \
	simulate-CompositeJunction

check_PROGRAMS += \
	simulate-SymOneSite \
	simulate-AsymOneSite \
	simulate-SymTwoSite \
	simulate-AsymTwoSite \
	simulate-SymInterference \
	simulate-RectBarrier \
	simulate-CompositeJunction

simulate_SymOneSite_SOURCES = simulate-SymOneSite.cc
simulate_SymOneSite_LDADD = ../simulator_models/libtransport_simulate.a \
	../../general/libmolstat_simulator.a \
	../../general/libmolstat_general.a 
if HAVE_GSL
simulate_SymOneSite_LDADD += \
	$(GSL_LDFLAGS) $(GSL_LIBS)
endif # HAVE_GSL

simulate_AsymOneSite_SOURCES = simulate-AsymOneSite.cc
simulate_AsymOneSite_LDADD = ../simulator_models/libtransport_simulate.a \
	../../general/libmolstat_simulator.a \
 	../../general/libmolstat_general.a 
if HAVE_GSL
simulate_AsymOneSite_LDADD += \
	$(GSL_LDFLAGS) $(GSL_LIBS)
endif # HAVE_GSL

simulate_SymTwoSite_SOURCES = simulate-SymTwoSite.cc
simulate_SymTwoSite_LDADD = ../simulator_models/libtransport_simulate.a \
	../../general/libmolstat_simulator.a \
	../../general/libmolstat_general.a
if HAVE_GSL
simulate_SymTwoSite_LDADD += \
	$(GSL_LDFLAGS) $(GSL_LIBS)
endif # HAVE_GSL

simulate_AsymTwoSite_SOURCES = simulate-AsymTwoSite.cc
simulate_AsymTwoSite_LDADD = ../simulator_models/libtransport_simulate.a \
	../../general/libmolstat_simulator.a \
	../../general/libmolstat_general.a
if HAVE_GSL
simulate_AsymTwoSite_LDADD += \
	$(GSL_LDFLAGS) $(GSL_LIBS)
endif # HAVE_GSL

simulate_CompositeJunction_SOURCES = simulate-CompositeJunction.cc
simulate_CompositeJunction_LDADD = ../simulator_models/libtransport_simulate.a \
	../../general/libmolstat_simulator.a \
	../../general/libmolstat_general.a 
if HAVE_GSL
simulate_CompositeJunction_LDADD += \
	$(GSL_LDFLAGS) $(GSL_LIBS)
endif # HAVE_GSL

simulate_SymInterference_SOURCES = simulate-SymInterference.cc
simulate_SymInterference_LDADD = ../simulator_models/libtransport_simulate.a \
	../../general/libmolstat_simulator.a \
	../../general/libmolstat_general.a 
if HAVE_GSL
simulate_SymInterference_LDADD += \
	$(GSL_LDFLAGS) $(GSL_LIBS)
endif # HAVE_GSL

simulate_RectBarrier_SOURCES = simulate-RectBarrier.cc
simulate_RectBarrier_LDADD = ../simulator_models/libtransport_simulate.a \
	../../general/libmolstat_simulator.a \
	../../general/libmolstat_general.a
if HAVE_GSL
simulate_RectBarrier_LDADD += \
	$(GSL_LDFLAGS) $(GSL_LIBS)
endif # HAVE_GSL
endif # TRANSPORT_SIMULATOR

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
	fit-CompositeSymNonresonantBackground \
	fit-ExperimentSymNonresonant \
	fit-SymResonant \
	fit-Interference \
	fit-CompositeInterferenceBackground

check_PROGRAMS += \
	fit-AsymResonant \
	fit-SymNonresonant \
	fit-CompositeSymNonresonantBackground \
	fit-ExperimentSymNonresonant \
	fit-SymResonant \
	fit-Interference \
	fit-CompositeInterferenceBackground

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

fit_CompositeSymNonresonantBackground_SOURCES = fit-CompositeSymNonresonantBackground.cc
fit_CompositeSymNonresonantBackground_CPPFLAGS = $(GSL_INCLUDE) $(AM_CPPFLAGS)
fit_CompositeSymNonresonantBackground_LDADD = ../fitter_models/libtransport_fit.a \
	../../general/libmolstat_fitter.a \
	../../general/libmolstat_general.a \
	$(GSL_LDFLAGS) $(AM_LDADD) $(GSL_LIBS) $(AM_LIBS)

fit_ExperimentSymNonresonant_SOURCES = fit-ExperimentSymNonresonant.cc
fit_ExperimentSymNonresonant_CPPFLAGS = $(GSL_INCLUDE) $(AM_CPPFLAGS)
fit_ExperimentSymNonresonant_LDADD = ../fitter_models/libtransport_fit.a \
	../../general/libmolstat_fitter.a \
	../../general/libmolstat_general.a \
	$(GSL_LDFLAGS) $(AM_LDADD) $(GSL_LIBS) $(AM_LIBS)

fit_CompositeInterferenceBackground_SOURCES = fit-CompositeInterferenceBackground.cc
fit_CompositeInterferenceBackground_CPPFLAGS = $(GSL_INCLUDE) $(AM_CPPFLAGS)
fit_CompositeInterferenceBackground_LDADD = ../fitter_models/libtransport_fit.a \
	../../general/libmolstat_fitter.a \
	../../general/libmolstat_general.a \
	$(GSL_LDFLAGS) $(AM_LDADD) $(GSL_LIBS) $(AM_LIBS)

fit_Interference_SOURCES = fit-Interference.cc
fit_Interference_CPPFLAGS = $(GSL_INCLUDE) $(AM_CPPFLAGS)
fit_Interference_LDADD = ../fitter_models/libtransport_fit.a \
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

# make sure automake includes the script in a distribution
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
