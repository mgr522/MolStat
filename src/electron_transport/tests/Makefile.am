TESTS = simulate-SymOneSite \
	simulate-AsymOneSite \
	simulate-SymTwoSite \
	simulate-AsymTwoSite

if HAVE_PYTHON
TESTS += fit-symmetric-resonant.py \
	fit-symmetric-nonresonant.py \
   fit-symmetric-nonresonant-binlog.py \
	fit-asymmetric-resonant.py \
	fit-options.py
endif

TEST_EXTENSIONS = .py
PY_LOG_COMPILER = $(PYTHON_BIN)

check_PROGRAMS = simulate-SymOneSite \
	simulate-AsymOneSite \
	simulate-SymTwoSite \
	simulate-AsymTwoSite

simulate_SymOneSite_SOURCES = simulate-SymOneSite.cc
simulate_SymOneSite_LDADD = ../../general/libmolstat_general.a \
	../simulator_models/libtransport_simulator_models.a

simulate_AsymOneSite_SOURCES = simulate-AsymOneSite.cc
simulate_AsymOneSite_LDADD = ../../general/libmolstat_general.a \
	../simulator_models/libtransport_simulator_models.a

simulate_SymTwoSite_SOURCES = simulate-SymTwoSite.cc
simulate_SymTwoSite_LDADD = ../../general/libmolstat_general.a \
	../simulator_models/libtransport_simulator_models.a

simulate_AsymTwoSite_SOURCES = simulate-AsymTwoSite.cc
simulate_AsymTwoSite_LDADD = ../../general/libmolstat_general.a \
	../simulator_models/libtransport_simulator_models.a
