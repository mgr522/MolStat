TESTS = test_symmetric_one_site

if HAVE_PYTHON
TESTS += fit-symmetric-resonant.py \
	fit-symmetric-nonresonant.py \
   fit-symmetric-nonresonant-binlog.py \
	fit-asymmetric-resonant.py \
	fit-options.py
endif

TEST_EXTENSIONS = .py
PY_LOG_COMPILER = $(PYTHON_BIN)

check_PROGRAMS = test_symmetric_one_site

test_symmetric_one_site_SOURCES = symmetric_one_site.cc
test_symmetric_one_site_LDADD = ../../general/libmolstat_general.a \
	../simulator_models/libtransport_simulator_models.a

#test_asymmetric_one_site_SOURCES = asymmetric_one_site.cc
#test_asymmetric_one_site_LDADD = ../simulator_models/libsimulator_models.a

#test_symmetric_two_site_SOURCES = symmetric_two_site.cc
#test_symmetric_two_site_LDADD = ../simulator_models/libsimulator_models.a

#test_asymmetric_two_site_SOURCES = asymmetric_two_site.cc
#test_asymmetric_two_site_LDADD = ../simulator_models/libsimulator_models.a

#test_symmetric_voltage_two_site_SOURCES = symmetric_voltage_two_site.cc
#test_symmetric_voltage_two_site_LDADD = ../simulator_models/libsimulator_models.a
