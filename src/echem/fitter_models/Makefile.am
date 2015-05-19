AM_CPPFLAGS = -I$(top_srcdir)/src

if ECHEM_FITTER
# electrochemistry fitter module requires GSL

noinst_LIBRARIES = libechem_fit.a

libechem_fit_a_SOURCES = echem_fit_module.h \
	echem_fit_module.cc

# electrochemistry fit module requires GSL
libechem_fit_a_CPPFLAGS = $(GSL_INCLUDE) $(AM_CPPFLAGS)
endif
