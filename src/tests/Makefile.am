TESTS = test_histogram2d \
	test_symmetric_one_site \
	test_asymmetric_one_site
check_PROGRAMS = test_histogram2d \
	test_symmetric_one_site \
	test_asymmetric_one_site

test_histogram2d_SOURCES = histogram2d.cc
test_histogram2d_LDADD =  ../aux_binner/libaux_binner.a

test_symmetric_one_site_SOURCES = symmetric_one_site.cc
test_symmetric_one_site_LDADD = ../aux_simulator/libaux_simulator.a

test_asymmetric_one_site_SOURCES = asymmetric_one_site.cc
test_asymmetric_one_site_LDADD = ../aux_simulator/libaux_simulator.a
