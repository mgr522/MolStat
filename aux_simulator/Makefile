CPP=clang++ -std=c++11
CXXFLAGS=-O2 -Wall

# GSL on MAC
GSL_DIR=/opt/local
GSL_LIB=$(GSL_DIR)/lib
GSL_INCLUDE=$(GSL_DIR)/include

INCLUDE=-I$(GSL_INCLUDE)
LIBS=-L$(GSL_LIB) -lgsl -lgslcblas -lm

all: libaux_simulator.a

libaux_simulator.a: rng.o
	ar rvs libaux_simulator.a rng.o

rng.o: rng.h rng.cc
	$(CPP) $(CXXFLAGS) -c -o rng.o rng.cc $(INCLUDE)

clean:
	rm -f libaux_simulator.a *.o
