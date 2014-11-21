TESTS = string_tools \
	counter_index_functionality \
	histogram1d_linear \
	histogram1d_log \
	histogram2d_linear \
	histogram2d_mixed \
	histogram2d_log

check_PROGRAMS = string_tools \
	counter_index_functionality \
	histogram1d_linear \
	histogram1d_log \
	histogram2d_linear \
	histogram2d_mixed \
	histogram2d_log

string_tools_SOURCES = string_tools.cc
string_tools_LDADD = ../libmolstat_general.a

counter_index_functionality_SOURCES = counter_index_functionality.cc
counter_index_functionality_LDADD = ../libmolstat_general.a

histogram1d_linear_SOURCES = histogram1d_linear.cc
histogram1d_linear_LDADD = ../libmolstat_general.a

histogram1d_log_SOURCES = histogram1d_log.cc
histogram1d_log_LDADD = ../libmolstat_general.a

histogram2d_linear_SOURCES = histogram2d_linear.cc
histogram2d_linear_LDADD = ../libmolstat_general.a

histogram2d_mixed_SOURCES = histogram2d_mixed.cc
histogram2d_mixed_LDADD = ../libmolstat_general.a

histogram2d_log_SOURCES = histogram2d_log.cc
histogram2d_log_LDADD = ../libmolstat_general.a

if BUILD_SIMULATOR
TESTS += \
	simulate_model_interface_direct \
	simulate_model_interface_indirect

check_PROGRAMS += \
	simulate_model_interface_direct \
	simulate_model_interface_indirect

simulate_model_interface_direct_SOURCES = \
	simulate_model_interface_observables.h \
	simulate_model_interface_models.h \
	simulate_model_interface_direct.cc
simulate_model_interface_direct_LDADD = \
	../libmolstat_simulator.a \
	../libmolstat_general.a

simulate_model_interface_indirect_SOURCES = \
	simulate_model_interface_observables.h \
	simulate_model_interface_models.h \
	simulate_model_interface_indirect.cc
simulate_model_interface_indirect_LDADD = \
	../libmolstat_simulator.a \
	../libmolstat_general.a
endif

if BUILD_FITTER
TESTS += gsl_std_vector
check_PROGRAMS += gsl_std_vector

gsl_std_vector_SOURCES = gsl_std_vector.cc
gsl_std_vector_LDADD = \
	../libmolstat_fitter.a \
	../libmolstat_general.a \
	$(GSL_LDFLAGS) $(AM_LDFLAGS) $(GSL_LIBS) $(AM_LIBS)
gsl_std_vector_CPPFLAGS = $(GSL_INCLUDE) $(AM_CPPFLAGS)
endif
