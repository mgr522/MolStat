AM_CPPFLAGS = -I$(top_srcdir)/src

if TRANSPORT_SIMULATOR
noinst_LIBRARIES = libtransport_simulate.a

libtransport_simulate_a_SOURCES = transport_simulate_module.h \
	transport_simulate_module.cc \
	observables.h \
	junction.h \
	junction.cc \
	sym_one_site_channel.h \
	sym_one_site_channel.cc \
	asym_one_site_channel.h \
	asym_one_site_channel.cc \
	sym_two_site_channel.h \
	sym_two_site_channel.cc \
	asym_two_site_channel.h \
	asym_two_site_channel.cc \
	rectangular_barrier.h \
	rectangular_barrier.cc \
	sym_interference.h \
	sym_interference.cc
endif
