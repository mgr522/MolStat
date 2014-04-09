CPP=clang++ -std=c++11
CXXFLAGS=-O2 -Wall

# Compile on MAC
GSL_DIR=/opt/local
GSL_LIB=$(GSL_DIR)/lib
GSL_INCLUDE=$(GSL_DIR)/include

INCLUDE=-I$(GSL_INCLUDE)
LIBS=-L$(GSL_LIB) -lgsl -lgslcblas -lm

all: simulator binner simulator-v-2d binner-v-2d
#all: simulator binner fitter simulator-v-2d binner-v-2d
	#cp simulator binner fitter simulator-v-2d binner-v-2d ../bin

simulator: main-simulator.cc
	$(CPP) $(CXXFLAGS) -o simulator main-simulator.cc \
		$(INCLUDE) $(LIBS)

simulator-v-2d: main-simulator-v-2d.cc aux_simulator
	cd aux_simulator; make
	$(CPP) $(CXXFLAGS) -o simulator-v-2d main-simulator-v-2d.cc \
		aux_simulator/libaux_simulator.a $(INCLUDE) $(LIBS)

binner: main-binner.cc
	$(CPP) $(CXXFLAGS) -o binner main-binner.cc $(INCLUDE) $(LIBS)

binner-v-2d: main-binner-v-2d.cc
	$(CPP) $(CXXFLAGS) -o binner-v-2d main-binner-v-2d.cc $(INCLUDE) $(LIBS)

fitter: main-fitter.cc models.h model-asymmetric-resonant.h model-asymmetric-resonant.cc model-symmetric-nonresonant.h model-symmetric-nonresonant.cc model-symmetric-resonant.h model-symmetric-resonant.cc
	$(CPP) $(CXXFLAGS) -o fitter main-fitter.cc model-asymmetric-resonant.cc model-symmetric-nonresonant.cc model-symmetric-resonant.cc $(INCLUDE) $(LIBS)

clean:
	cd aux_simulator; make clean
	rm -f *.o simulator binner fitter simulator-v-2d binner-v-2d

distclean:
	rm -f *.o simulator binner fitter simulator-v-2d binner-v-2d ../bin/simulator ../bin/binner ../bin/fitter ../bin/simulator-v-2d ../bin/binner-v-2d
