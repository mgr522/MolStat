CPP=clang++

# Compile on MAC
GSL_DIR=/opt/local
GSL_LIB=$(GSL_DIR)/lib
GSL_INCLUDE=$(GSL_DIR)/include

CFLAGS=-O2 -Wall -I$(GSL_INCLUDE)
LIBS=-L$(GSL_LIB) -lgsl -lgslcblas -lm

all: simulator binner simulator-v-2d binner-v-2d
#all: simulator binner fitter simulator-v-2d binner-v-2d
	#cp simulator binner fitter simulator-v-2d binner-v-2d ../bin

simulator: main-simulator.cc
	$(CPP) -o simulator main-simulator.cc \
		$(CFLAGS) $(LIBS)

simulator-v-2d: main-simulator-v-2d.cc
	$(CPP) -o simulator-v-2d main-simulator-v-2d.cc \
		$(CFLAGS) $(LIBS)

binner: main-binner.cc
	$(CPP) -o binner main-binner.cc $(CFLAGS) $(LIBS)

binner-v-2d: main-binner-v-2d.cc
	$(CPP) -o binner-v-2d main-binner-v-2d.cc $(CFLAGS) $(LIBS)

fitter: main-fitter.cc models.h model-asymmetric-resonant.h model-asymmetric-resonant.cc model-symmetric-nonresonant.h model-symmetric-nonresonant.cc model-symmetric-resonant.h model-symmetric-resonant.cc
	$(CPP) -o fitter main-fitter.cc model-asymmetric-resonant.cc model-symmetric-nonresonant.cc model-symmetric-resonant.cc $(CFLAGS) $(LIBS)

clean:
	rm -f *.o simulator binner fitter simulator-v-2d binner-v-2d

distclean:
	rm -f *.o simulator binner fitter simulator-v-2d binner-v-2d ../bin/simulator ../bin/binner ../bin/fitter ../bin/simulator-v-2d ../bin/binner-v-2d
