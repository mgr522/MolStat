CPP=clang++ -std=c++11
CXXFLAGS=-O2 -Wall

# Compile on MAC
GSL_DIR=/opt/local
GSL_LIB=$(GSL_DIR)/lib
GSL_INCLUDE=$(GSL_DIR)/include

INCLUDE=-I$(GSL_INCLUDE)
LIBS=-L$(GSL_LIB) -lgsl -lgslcblas -lm

all: simulator binner binner-v-2d
#all: simulator binner fitter simulator-v-2d binner-v-2d
	#cp simulator binner fitter simulator-v-2d binner-v-2d ../bin

simulator: main-simulator.cc aux_simulator string_tools.o
	cd aux_simulator; make
	$(CPP) $(CXXFLAGS) -o simulator main-simulator.cc string_tools.o \
		aux_simulator/libaux_simulator.a $(INCLUDE) $(LIBS)

binner: main-binner.cc
	$(CPP) $(CXXFLAGS) -o binner main-binner.cc $(INCLUDE) $(LIBS)

binner-v-2d: main-binner-v-2d.cc
	$(CPP) $(CXXFLAGS) -o binner-v-2d main-binner-v-2d.cc $(INCLUDE) $(LIBS)

fitter: main-fitter.cc models.h model-asymmetric-resonant.h model-asymmetric-resonant.cc model-symmetric-nonresonant.h model-symmetric-nonresonant.cc model-symmetric-resonant.h model-symmetric-resonant.cc
	$(CPP) $(CXXFLAGS) -o fitter main-fitter.cc model-asymmetric-resonant.cc model-symmetric-nonresonant.cc model-symmetric-resonant.cc $(INCLUDE) $(LIBS)

string_tools.o: string_tools.h string_tools.cc
	$(CPP) $(CXXFLAGS) -c -o string_tools.o string_tools.cc $(INCLUDE)

clean:
	cd aux_simulator; make clean
	rm -f *.o simulator binner fitter binner-v-2d

distclean:
	rm -f *.o simulator binner fitter simulator-v-2d binner-v-2d ../bin/simulator ../bin/binner ../bin/fitter ../bin/simulator-v-2d ../bin/binner-v-2d
