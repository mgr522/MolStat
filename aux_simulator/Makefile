CPP=clang++ -std=c++11
CXXFLAGS=-O2 -Wall

# GSL on MAC
GSL_DIR=/opt/local
GSL_LIB=$(GSL_DIR)/lib
GSL_INCLUDE=$(GSL_DIR)/include

INCLUDE=-I$(GSL_INCLUDE)
LIBS=-L$(GSL_LIB) -lgsl -lgslcblas -lm

all: libaux_simulator.a

libaux_simulator.a: rng.o model_interface.o voltage_independent.o
	ar rvs libaux_simulator.a rng.o model_interface.o voltage_independent.o

rng.o: rng.h rng.cc
	$(CPP) $(CXXFLAGS) -c -o rng.o rng.cc $(INCLUDE)

model_interface.o: model_interface.h model_interface.cc rng.h
	$(CPP) $(CXXFLAGS) -c -o model_interface.o model_interface.cc $(INCLUDE)

voltage_independent.o: voltage_independent.h voltage_independent.cc model_interface.h rng.h
	$(CPP) $(CXXFLAGS) -c -o voltage_independent.o voltage_independent.cc $(INCLUDE)

clean:
	rm -f libaux_simulator.a *.o
