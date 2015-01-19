if TRANSPORT_FITTER
noinst_LIBRARIES = libtransport_fit.a

libtransport_fit_a_SOURCES = transport_fit_module.h \
	transport_fit_module.cc \
	symmetric_resonant.h \
	symmetric_resonant.cc \
	symmetric_nonresonant.h \
	symmetric_nonresonant.cc \
	asymmetric_resonant.h \
	asymmetric_resonant.cc \
	symmetric_nonresonant_w_vacuum.h \
	symmetric_nonresonant_w_vacuum.cc

# transport fit module requires GSL
libtransport_fit_a_CPPFLAGS = $(GSL_INCLUDE) $(AM_CPPFLAGS)
endif
