TESTS = test_histogram1d_linear \
	test_histogram1d_log \
	test_histogram2d_linear \
	test_histogram2d_log \
	test_symmetric_one_site \
	test_asymmetric_one_site \
	test_symmetric_voltage_one_site \
	test_symmetric_two_site \
	test_asymmetric_two_site \
	test_symmetric_voltage_two_site \
	test_gsl_std_vector
check_PROGRAMS = test_histogram1d_linear \
	test_histogram1d_log \
	test_histogram2d_linear \
	test_histogram2d_log \
	test_symmetric_one_site \
	test_asymmetric_one_site \
	test_symmetric_voltage_one_site \
	test_symmetric_two_site \
	test_asymmetric_two_site \
	test_symmetric_voltage_two_site \
	test_gsl_std_vector

test_histogram1d_linear_SOURCES = histogram1d_linear.cc
test_histogram1d_linear_LDADD =  ../aux_binner/libaux_binner.a

test_histogram1d_log_SOURCES = histogram1d_log.cc
test_histogram1d_log_LDADD =  ../aux_binner/libaux_binner.a

test_histogram2d_linear_SOURCES = histogram2d_linear.cc
test_histogram2d_linear_LDADD =  ../aux_binner/libaux_binner.a

test_histogram2d_log_SOURCES = histogram2d_log.cc
test_histogram2d_log_LDADD =  ../aux_binner/libaux_binner.a

test_symmetric_one_site_SOURCES = symmetric_one_site.cc
test_symmetric_one_site_LDADD = ../aux_simulator/libaux_simulator.a

test_asymmetric_one_site_SOURCES = asymmetric_one_site.cc
test_asymmetric_one_site_LDADD = ../aux_simulator/libaux_simulator.a

test_symmetric_voltage_one_site_SOURCES = symmetric_voltage_one_site.cc
test_symmetric_voltage_one_site_LDADD = ../aux_simulator/libaux_simulator.a

test_symmetric_two_site_SOURCES = symmetric_two_site.cc
test_symmetric_two_site_LDADD = ../aux_simulator/libaux_simulator.a

test_asymmetric_two_site_SOURCES = asymmetric_two_site.cc
test_asymmetric_two_site_LDADD = ../aux_simulator/libaux_simulator.a

test_symmetric_voltage_two_site_SOURCES = symmetric_voltage_two_site.cc
test_symmetric_voltage_two_site_LDADD = ../aux_simulator/libaux_simulator.a

test_gsl_std_vector_SOURCES = gsl_std_vector.cc
test_gsl_std_vector_LDADD = ../aux_fitter/libaux_fitter.a
