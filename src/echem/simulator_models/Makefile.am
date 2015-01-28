if ECHEM_SIMULATOR
# parts of the electrochemistry simulator module require CVODE
noinst_LIBRARIES = libechem_simulate.a

libechem_simulate_a_SOURCES = \
	echem_simulate_module.h \
	echem_simulate_module.cc \
	observables.h \
	nernstian.h \
	nernstian.cc 

libechem_simulate_a_CPPFLAGS =

# parts that require CVODE
if HAVE_CVODE
libechem_simulate_a_SOURCES += \
	non-nernstian.h \
	non-nernstian.cc

libechem_simulate_a_CPPFLAGS += $(CVODE_INCLUDE)
endif # HAVE_CVODE

libechem_simulate_a_CPPFLAGS += $(AM_CPPFLAGS)

endif # ECHEM_SIMULATOR
