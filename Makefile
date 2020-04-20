all: lac-test

LIBGEOM=/home/salvi/project/libgeom
INCLUDES=-I$(LIBGEOM)
LIBS=-L$(LIBGEOM)/release -lgeom

CXXFLAGS=-Wall -std=c++17 -pedantic -g $(INCLUDES)

lac-test: lac-test.o discrete-lac.o
	$(CXX) -o $@ $^ $(LIBS)
