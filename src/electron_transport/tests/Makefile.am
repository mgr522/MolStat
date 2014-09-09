TESTS = simulate-SymOneSite \
	simulate-AsymOneSite \
	simulate-SymTwoSite

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
	simulate-SymTwoSite

simulate_SymOneSite_SOURCES = simulate-SymOneSite.cc
simulate_SymOneSite_LDADD = ../../general/libmolstat_general.a \
	../simulator_models/libtransport_simulator_models.a

simulate_AsymOneSite_SOURCES = simulate-AsymOneSite.cc
simulate_AsymOneSite_LDADD = ../../general/libmolstat_general.a \
	../simulator_models/libtransport_simulator_models.a

simulate_SymTwoSite_SOURCES = simulate-SymTwoSite.cc
simulate_SymTwoSite_LDADD = ../../general/libmolstat_general.a \
	../simulator_models/libtransport_simulator_models.a

#test_asymmetric_one_site_SOURCES = asymmetric_one_site.cc
#test_asymmetric_one_site_LDADD = ../simulator_models/libsimulator_models.a

#test_symmetric_two_site_SOURCES = symmetric_two_site.cc
#test_symmetric_two_site_LDADD = ../simulator_models/libsimulator_models.a

#test_asymmetric_two_site_SOURCES = asymmetric_two_site.cc
#test_asymmetric_two_site_LDADD = ../simulator_models/libsimulator_models.a

#test_symmetric_voltage_two_site_SOURCES = symmetric_voltage_two_site.cc
#test_symmetric_voltage_two_site_LDADD = ../simulator_models/libsimulator_models.a
