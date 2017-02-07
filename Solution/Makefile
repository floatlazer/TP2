# Makefile TD n°2
CXX = g++ -std=c++14
CXXFLAGS = -O2 -Wall -pedantic -march=native
#CXXFLAGS = -g -Wall -pedantic -march=native
LIBS = -lm

all:    Mandelbrot.exe matvec.exe

Mandelbrot.exe: Mandelbrot.cpp lodepng/lodepng.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

matvec.exe:	matvec.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

cleanall:
	@rm -rf *.o *~ *.exe
