CPP=clang++ -std=c++11
CXXFLAGS=-O2 -Wall

# GSL on MAC
GSL_DIR=/opt/local
GSL_LIB=$(GSL_DIR)/lib
GSL_INCLUDE=$(GSL_DIR)/include

INCLUDE=-I$(GSL_INCLUDE)
LIBS=-L$(GSL_LIB) -lgsl -lgslcblas -lm

all: libaux_simulator.a

libaux_simulator.a: rng.o symmetric_voltage_independent.o symmetric_voltage_one_site.o
	ar rvs libaux_simulator.a rng.o symmetric_voltage_independent.o symmetric_voltage_one_site.o

rng.o: rng.h rng.cc
	$(CPP) $(CXXFLAGS) -c -o rng.o rng.cc $(INCLUDE)

symmetric_voltage_independent.o: symmetric_voltage_independent.h symmetric_voltage_independent.cc model_interface.h rng.h
	$(CPP) $(CXXFLAGS) -c -o symmetric_voltage_independent.o symmetric_voltage_independent.cc $(INCLUDE)

symmetric_voltage_one_site.o: symmetric_voltage_one_site.h symmetric_voltage_one_site.cc symmetric_voltage_independent.h model_interface.h rng.h
	$(CPP) $(CXXFLAGS) -c -o symmetric_voltage_one_site.o symmetric_voltage_one_site.cc $(INCLUDE)

clean:
	rm -f libaux_simulator.a *.o
