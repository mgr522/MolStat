TESTS = gsl_std_vector \
	histogram1d_linear \
	histogram1d_log \
	histogram2d_linear \
	histogram2d_log \
	simulate_model_interface_direct \
	simulate_model_interface_indirect

check_PROGRAMS = gsl_std_vector \
	histogram1d_linear \
	histogram1d_log \
	histogram2d_linear \
	histogram2d_log \
	simulate_model_interface_direct \
	simulate_model_interface_indirect

gsl_std_vector_SOURCES = gsl_std_vector.cc
gsl_std_vector_LDADD = ../libmolstat_general.a

histogram1d_linear_SOURCES = histogram1d_linear.cc
histogram1d_linear_LDADD = ../libmolstat_general.a

histogram1d_log_SOURCES = histogram1d_log.cc
histogram1d_log_LDADD = ../libmolstat_general.a

histogram2d_linear_SOURCES = histogram2d_linear.cc
histogram2d_linear_LDADD = ../libmolstat_general.a

histogram2d_log_SOURCES = histogram2d_log.cc
histogram2d_log_LDADD = ../libmolstat_general.a

simulate_model_interface_direct_SOURCES = \
	simulate_model_interface_observables.h \
	simulate_model_interface_models.h \
	simulate_model_interface_direct.cc
simulate_model_interface_direct_LDADD = ../libmolstat_general.a

simulate_model_interface_indirect_SOURCES = \
	simulate_model_interface_observables.h \
	simulate_model_interface_models.h \
	simulate_model_interface_indirect.cc
simulate_model_interface_indirect_LDADD = ../libmolstat_general.a
