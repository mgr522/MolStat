TESTS = test_histogram1d_linear \
	test_histogram1d_log \
	test_histogram2d_linear \
	test_histogram2d_log \
	test_symmetric_one_site \
	test_asymmetric_one_site \
	test_symmetric_voltage_one_site \
	test_symmetric_two_site \
	test_asymmetric_two_site \
	test_symmetric_voltage_two_site

if HAVE_PYTHON
TESTS += fit-symmetric-resonant.py \
	fit-symmetric-nonresonant.py \
	fit-asymmetric-resonant.py
endif

TEST_EXTENSIONS = .py
PY_LOG_COMPILER = $(PYTHON_BIN)

check_PROGRAMS = test_histogram1d_linear \
	test_histogram1d_log \
	test_histogram2d_linear \
	test_histogram2d_log \
	test_symmetric_one_site \
	test_asymmetric_one_site \
	test_symmetric_voltage_one_site \
	test_symmetric_two_site \
	test_asymmetric_two_site \
	test_symmetric_voltage_two_site

test_histogram1d_linear_SOURCES = histogram1d_linear.cc
test_histogram1d_linear_LDADD = ../../general/libmolstat_general.a

test_histogram1d_log_SOURCES = histogram1d_log.cc
test_histogram1d_log_LDADD = ../../general/libmolstat_general.a

test_histogram2d_linear_SOURCES = histogram2d_linear.cc
test_histogram2d_linear_LDADD = ../../general/libmolstat_general.a

test_histogram2d_log_SOURCES = histogram2d_log.cc
test_histogram2d_log_LDADD = ../../general/libmolstat_general.a

test_symmetric_one_site_SOURCES = symmetric_one_site.cc
test_symmetric_one_site_LDADD = ../simulator_models/libsimulator_models.a

test_asymmetric_one_site_SOURCES = asymmetric_one_site.cc
test_asymmetric_one_site_LDADD = ../simulator_models/libsimulator_models.a

test_symmetric_voltage_one_site_SOURCES = symmetric_voltage_one_site.cc
test_symmetric_voltage_one_site_LDADD = ../simulator_models/libsimulator_models.a

test_symmetric_two_site_SOURCES = symmetric_two_site.cc
test_symmetric_two_site_LDADD = ../simulator_models/libsimulator_models.a

test_asymmetric_two_site_SOURCES = asymmetric_two_site.cc
test_asymmetric_two_site_LDADD = ../simulator_models/libsimulator_models.a

test_symmetric_voltage_two_site_SOURCES = symmetric_voltage_two_site.cc
test_symmetric_voltage_two_site_LDADD = ../simulator_models/libsimulator_models.a
